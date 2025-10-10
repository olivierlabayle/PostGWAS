module PostGWAS

using CSV
using DataFrames
using Parquet
using TMLE
using Combinatorics
using RCall
using Downloads
using JSON

function extract_variant_info(variant_id)
    chr, pos, ref, alt, assembly = split(variant_id, "_")
    return (CHR=replace(chr, "chr" => ""), POS=parse(Int, pos), REF=ref, ALT=alt)
end

join_variant_info(chr_col, pos_col, ref_col, alt_col) = 
    string.(string.(chr_col), ":", string.(pos_col), ":", string.(ref_col), ":", string.(alt_col))

function add_universal_ID!(df; cols=[:CHROM, :POS, :REF, :ALT])
    transform!(df, 
        cols => join_variant_info => :UNIV_ID)
    return df
end

function get_harmonized_qtl_data!(parquet_file)
    eqtl_data = read_parquet(parquet_file) |> DataFrame
    tissue = split(basename(parquet_file), ".")[1]
    transform!(eqtl_data, :variant_id => ByRow(extract_variant_info) => AsTable)
    add_universal_ID!(eqtl_data)
    eqtl_data.QTL_TISSUE .= tissue
    return eqtl_data
end

function get_harmonized_finemapping_results(finemapping_file)
    finemapping_results = CSV.read(finemapping_file, DataFrame, delim="\t")
    rename!(finemapping_results, "#CHROM" => "CHR")
    add_universal_ID!(finemapping_results)
    return finemapping_results
end

function get_finemapped_loci(finemapping_file)
    finemapping_df = CSV.read(finemapping_file, DataFrame)
    add_universal_ID!(finemapping_df)
    return groupby(finemapping_df, :LOCUS_ID)
end

function GTEx_columns_selection(;qtl_type=:eQTL)
    if qtl_type == :eQTL
        return [:UNIV_ID,
                :QTL_TISSUE,
                :gene_id => :QTL_GENE_ID, 
                :pval_nominal => :QTL_PVAL, 
                :slope => :QTL_BETA, 
                :slope_se => :QTL_BETA_SE]
    elseif qtl_type == :sQTL
        return [:UNIV_ID,
                :QTL_TISSUE,
                :phenotype_id => :QTL_PHENOTYPE_ID,
                :gene_id => :QTL_GENE_ID, 
                :pval_nominal => :QTL_PVAL, 
                :slope => :QTL_BETA, 
                :slope_se => :QTL_BETA_SE]
    else
        error("Unsupported QTL type: $qtl_type")
    end
end

function get_GTEx_genes(variants, parquet_files; qtl_type=:eQTL)
    qtls_columns = GTEx_columns_selection(;qtl_type=qtl_type)
    variant_target_pairs = []
    for parquet_file in parquet_files
        sig_cis_eqtls = get_harmonized_qtl_data!(parquet_file)
        merged_results = innerjoin(
            select(variants, [:UNIV_ID]),
            select(sig_cis_eqtls, qtls_columns...), 
            on =:UNIV_ID
        )
        push!(variant_target_pairs, merged_results)
    end
    qtl_data = vcat(variant_target_pairs...)
    qtl_data.QTL_TYPE .= qtl_type
    return qtl_data
end

function map_variants_to_targets(GTEx_dirs, variants)
    eQTLs_GTEx_dir, sQTLs_GTEx_dir = GTEx_dirs
    GTEx_eQTLs = PostGWAS.get_GTEx_genes(locus, filter(endswith(".parquet"), readdir(eQTLs_GTEx_dir, join=true)); qtl_type=:eQTL)
    GTEx_sQTLs = get_GTEx_genes(variants, filter(endswith(".parquet"), readdir(sQTLs_GTEx_dir, join=true)); qtl_type=:sQTL)

end

function download_GTex_QTLs(;file="GTEx_Analysis_v10_eQTL.tar", gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data"))
    tar_file = joinpath(gtex_cache_dir, file)
    outputdir = joinpath(gtex_cache_dir, replace(file, ".tar" => "_updated"))
    if !isdir(outputdir)
        @info "Downloading GETx $file, this may take a while..."
        Downloads.download("https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/$file", tar_file)
        run(`tar -xf $(tar_file) -C $(gtex_cache_dir)`)
        rm(tar_file)
    end
    return outputdir
end

function download_all_GTEx_QTLs_v10(;gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data"))
    if !isdir(gtex_cache_dir)
        mkpath(gtex_cache_dir)
    end
    eQTLs_GTEx_dir = download_GTex_QTLs(file="GTEx_Analysis_v10_eQTL.tar", gtex_cache_dir=gtex_cache_dir)
    sQTLs_GTEx_dir = download_GTex_QTLs(file="GTEx_Analysis_v10_sQTL.tar", gtex_cache_dir=gtex_cache_dir)
    return (eQTLs_GTEx_dir, sQTLs_GTEx_dir)
end

function get_filtered_gwas_results(gwas_file)
    gwas_results = CSV.read(gwas_file, DataFrame; select=["CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"])
    add_universal_ID!(gwas_results, cols=[:CHROM, :GENPOS, :ALLELE0, :ALLELE1])
    return select!(gwas_results, :UNIV_ID, :BETA, :SE, :LOG10P, :A1FREQ => :ALT_FREQ)
end

function download_ensembl_vep_cache(;vep_cache_dir=joinpath(ENV["HOME"], "vep_data"))
    @info "Installing Ensembl VEP and cache in $vep_cache_dir, this may take a while..."
    run(`docker run --rm -it --platform linux/amd64 -v $vep_cache_dir:/data ensemblorg/ensembl-vep INSTALL.pl -a cf -s homo_sapiens -y GRCh38`)
end

function recover_universal_IDs_from_ensembl_ann(ensembl_ann)
    chr, pos, varid, ref, alt, _ = split(ensembl_ann["input"], "\t")
    return string(chr, ":", pos, ":", ref, ":", alt)
end

function add_ensembl_annotations!(locus; vep_cache_dir=joinpath(ENV["HOME"], "vep_data"))
    tmpdir = mktempdir()
    ensembl_input = select(locus,
        "CHROM" => "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        )
    for col in ["QUAL", "FILTER", "INFO", "FORMAT"]
        ensembl_input[!, col] .= missing
    end
    CSV.write(joinpath(tmpdir, "input.vcf"), ensembl_input; delim="\t", missingstring=".")

    run(
        `docker run -v $vep_cache_dir:/data -v $tmpdir:/mnt/in_out ensemblorg/ensembl-vep \
            vep --cache --offline --format vcf --json --force_overwrite \
            --input_file /mnt/in_out/input.vcf \
            --output_file /mnt/in_out/output.json \
            --regulatory`
    )
    output_df = DataFrame(
        FULL_ENSEMBL_ANNOTATIONS = JSON.parse.(readlines(joinpath(tmpdir, "output.json")))
    )
    transform!(output_df, 
        "FULL_ENSEMBL_ANNOTATIONS" => (x -> recover_universal_IDs_from_ensembl_ann.(x)) => :UNIV_ID,
        "FULL_ENSEMBL_ANNOTATIONS" => (col -> get.(col, "most_severe_consequence", missing)) => "MOST_SEVERE_CONSEQUENCE",
    )
    leftjoin!(locus, output_df, on = :UNIV_ID)

    rm(tmpdir; force=true, recursive=true)
    return locus
end

function julia_main(gwas_file, finemapping_file, gtex_dir; 
    output_prefix="post_gwas",
    vep_cache_dir=joinpath(ENV["HOME"], "vep_data"),
    gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data")
    )
    # Download VEP cache data if not already done
    PostGWAS.download_ensembl_vep_cache(vep_cache_dir=vep_cache_dir)
    # Download GTEx QTL data if not already done
    GTEx_dirs = PostGWAS.download_all_GTEx_QTLs_v10(;gtex_cache_dir=gtex_cache_dir)
    # Load GWAS results
    gwas_results = PostGWAS.get_filtered_gwas_results(gwas_file)
    # For each finemapped locus:
    # - add GWAS results
    # - add Ensembl annotations
    # - map to GTEx eQTLs and sQTLs
    loci = PostGWAS.get_finemapped_loci(finemapping_file)
    for locus in loci
        leftjoin!(locus, gwas_results, on = :UNIV_ID)
        add_ensembl_annotations!(locus; vep_cache_dir=vep_cache_dir)
        map_variants_to_GTEx_genes(locus, parquet_files)
    end

    # Map finemapped variants to the genes they are eQTLs for
    
    CSV.write(string(output_prefix, ".finemapped_variant_egene_pairs.tsv"), variant_target_pairs; delim="\t", header=true)

    # Get interaction candidates
    egenes_with_more_than_one_variant = filter(
        :UNIV_IDs => x -> length(split(x, ",")) > 1,
        combine(
            groupby(variant_target_pairs, [:GENE_ID, :TISSUE]), 
            :UNIV_ID => (x -> join(x, ",")) => :UNIV_IDs
        )
    )
    egenes_with_more_than_one_variant.ESTIMANDS = map(egenes_with_more_than_one_variant.UNIV_IDs) do variant_str
        variants = split(variant_str, ",")
        map(combinations(variants, 2)) do variants_comb
            factorialEstimand(
                AIE,
                variants_comb,
                phenotype,
                dataset=dataset,
                confounders=confounders,
                outcome_extra_covariates=outcome_extra_covariates,
                positivity_constraint=positivity_constraint
            )
        end
    end
end

end
