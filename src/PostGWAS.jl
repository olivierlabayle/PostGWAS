module PostGWAS

using CSV
using DataFrames
using Parquet
using TMLE
using Combinatorics
using RCall
using Downloads

function extract_variant_info(variant_id)
    chr, pos, ref, alt, assembly = split(variant_id, "_")
    return (CHR=replace(chr, "chr" => ""), POS=parse(Int, pos), REF=ref, ALT=alt)
end

join_variant_info(chr_col, pos_col, ref_col, alt_col) = 
    string.(string.(chr_col), ":", string.(pos_col), ":", string.(ref_col), ":", string.(alt_col))

function add_universal_ID!(df; cols=[:CHR, :POS, :REF, :ALT])
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

    GTEx_eQTLs = get_GTEx_genes(variants, filter(endswith(".parquet"), readdir(eQTLs_GTEx_dir, join=true)); qtl_type=:eQTL)
    GTEx_sQTLs = get_GTEx_genes(variants, filter(endswith(".parquet"), readdir(sQTLs_GTEx_dir, join=true)); qtl_type=:sQTL)

end

function download_GTex_QTLs(;file="GTEx_Analysis_v10_eQTL.tar", output_dir="assets/")
    tar_file = joinpath(output_dir, file)
    outputdir = joinpath(output_dir, replace(file, ".tar" => "_updated"))
    if !isdir(outputdir)
        Downloads.download("https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/$file", tar_file)
        run(`tar -xf $(tar_file) -C $(output_dir)`)
        rm(tar_file)
    end
    return outputdir
end

function download_all_GTEx_QTLs_v10(;output_dir="assets/")
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    eQTLs_GTEx_dir = download_GTex_QTLs(file="GTEx_Analysis_v10_eQTL.tar", output_dir=output_dir)
    sQTLs_GTEx_dir = download_GTex_QTLs(file="GTEx_Analysis_v10_sQTL.tar", output_dir=output_dir)
    return (eQTLs_GTEx_dir, sQTLs_GTEx_dir)
end

function get_filtered_gwas_results(gwas_file; maf_threshold=0.005)
    gwas_results = CSV.read(gwas_file, DataFrame, delim="\t")
    return subset(gwas_results, :A1FREQ => x -> x .>= maf_threshold, skipmissing=true)
end

function get_filtered_finemapping_results(finemapping_file; pip_threshold=0.95)
    finemapping_results = get_harmonized_finemapping_results(finemapping_file)
    return subset(finemapping_results,
        :CS => x -> x .!== missing,
        :PIP => x -> x .> pip_threshold
    )
end

function julia_main(gwas_file, finemapping_file, gtex_dir; output_prefix="post_gwas")
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "post_gwas")
    phenotype = "SEVERE_COVID_19"
    gwas_file = "test/assets/results.all_chr.EUR.SEVERE_COVID_19.gwas.tsv"
    finemapping_file = "test/assets/results.all_chr.EUR.SEVERE_COVID_19.finemapping.tsv"
    gtex_dir = "../assets/GTEx_Analysis_v10_eQTL_updated"
    confounders = ["PC1", "PC2", "PC3", "PC4", "PC5"]
    outcome_extra_covariates = ["AGE", "SEX"]
    positivity_constraint = 0.01
    maf_threshold = 0.005
    pip_threshold = 0.
    workdir = "../assets/"


    gwas_results = get_filtered_gwas_results(gwas_file; maf_threshold=maf_threshold)
    finemapping_results = get_filtered_finemapping_results(finemapping_file; pip_threshold=pip_threshold)
    variants = innerjoin(
        gwas_results,
        select(finemapping_results, [:ID, :UNIV_ID, :PIP, :CS, :LOCUS_ID]),
        on = :ID
    )
    # Map finemapped variants to the genes they are eQTLs for
    GTEx_dirs = download_all_GTEx_QTLs_v10(;output_dir=workdir)
    variant_target_pairs = map_variants_to_GTEx_genes(finemapping_results, parquet_files)
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

end # module PostGWAS
