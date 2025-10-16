function variant_info_from_gtex_variant_id(variant_id)
    chr, pos, ref, alt, assembly = split(variant_id, "_")
    return (CHR=replace(chr, "chr" => ""), POS=parse(Int, pos), REF=ref, ALT=alt)
end

function get_harmonized_qtl_data(parquet_file)
    eqtl_data = read_parquet(parquet_file) |> DataFrame
    tissue = split(basename(parquet_file), ".")[1]
    transform!(eqtl_data, :variant_id => ByRow(variant_info_from_gtex_variant_id) => AsTable)
    add_universal_ID!(eqtl_data; cols=[:CHR, :POS, :REF, :ALT])
    eqtl_data.QTL_TISSUE .= tissue
    return eqtl_data
end

function GTEx_columns_selection(;qtl_type=:EQTL)
    if qtl_type == :EQTL
        return [:UNIV_ID => :UNIV_ID,
                :QTL_TISSUE => :QTL_TISSUE,
                :gene_id => :QTL_GENE_ID, 
                :pval_nominal => :QTL_PVAL, 
                :slope => :QTL_BETA, 
                :slope_se => :QTL_BETA_SE]
    elseif qtl_type == :SQTL
        return [:UNIV_ID => :UNIV_ID,
                :QTL_TISSUE => :QTL_TISSUE,
                :phenotype_id => :QTL_PHENOTYPE_ID,
                :gene_id => :QTL_GENE_ID, 
                :pval_nominal => :QTL_PVAL, 
                :slope => :QTL_BETA, 
                :slope_se => :QTL_BETA_SE]
    else
        error("Unsupported QTL type: $qtl_type")
    end
end

function get_qtl_data_matching_locus(locus, parquet_file, qtls_columns)
    sig_cis_eqtls = PostGWAS.get_harmonized_qtl_data(parquet_file)
    return innerjoin(
        select(locus, [:UNIV_ID]),
        select(sig_cis_eqtls, qtls_columns...), 
        on =:UNIV_ID
    )
end

function add_GTEx_qtls_info!(locus, parquet_files; qtl_type=:EQTL)
    qtls_columns = PostGWAS.GTEx_columns_selection(;qtl_type=qtl_type)
    # Get all TISSUE/TARGET/VARIANT association matching any of the variant in the locus
    tasks = map(parquet_files) do parquet_file
        Threads.@spawn PostGWAS.get_qtl_data_matching_locus(locus, parquet_file, qtls_columns)
    end
    variant_tissue_target_qtl_data = fetch.(tasks)
    variant_tissue_target_qtl_data = vcat(variant_tissue_target_qtl_data...)
    # Aggregate the infomation at the variant level
    info_cols = filter(!=(:UNIV_ID), last.(qtls_columns))
    variants_qtl_info = combine(groupby(variant_tissue_target_qtl_data, :UNIV_ID),
        AsTable(info_cols) =>
        (x -> Ref(Tables.columntable(x))) => Symbol(string("GTEX_", qtl_type, "_INFO"))
    )
    # Enrich the locus with the data
    return leftjoin!(locus, variants_qtl_info, on=:UNIV_ID)
end

function add_GTEX_info!(locus; gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data"))
    # Add eQTL info
    eqtls_dir = joinpath(gtex_cache_dir, "GTEx_Analysis_v10_eQTL_updated")
    parquet_files = filter(endswith(".parquet"), readdir(eqtls_dir, join=true))
    PostGWAS.add_GTEx_qtls_info!(locus, parquet_files; qtl_type=:EQTL)
    # Add sQTL info
    sqtls_dir = joinpath(gtex_cache_dir, "/Users/olabayle/gtex_data/GTEx_Analysis_v10_sQTL_updated")
    parquet_files = filter(endswith(".parquet"), readdir(sqtls_dir, join=true))
    PostGWAS.add_GTEx_qtls_info!(locus, parquet_files; qtl_type=:SQTL)
    return locus
end

function download_GTex_QTLs(;file="GTEx_Analysis_v10_eQTL.tar", gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data"))
    tar_file = joinpath(gtex_cache_dir, file)
    outputdir = joinpath(gtex_cache_dir, replace(file, ".tar" => "_updated"))
    if !isdir(outputdir)
        @info "Downloading GETx $file, this may take a while..."
        Downloads.download("https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/$file", tar_file)
        run(`tar -xf $(tar_file) -C $(gtex_cache_dir)`)
        rm(tar_file)
    else
        @info "GTEx $outputdir already exists, skipping."
    end
    return outputdir
end

function download_all_GTEx_QTLs_v10(;gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data"))
    @info "Downloading GTEx QTL data."
    if !isdir(gtex_cache_dir)
        mkpath(gtex_cache_dir)
    end
    eQTLs_GTEx_dir = download_GTex_QTLs(file="GTEx_Analysis_v10_eQTL.tar", gtex_cache_dir=gtex_cache_dir)
    sQTLs_GTEx_dir = download_GTex_QTLs(file="GTEx_Analysis_v10_sQTL.tar", gtex_cache_dir=gtex_cache_dir)
    return (eQTLs_GTEx_dir, sQTLs_GTEx_dir)
end
