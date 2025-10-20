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

function get_qtl_data_matching_locus(variant_ids, parquet_file, qtls_columns)
    sig_cis_eqtls = PostGWAS.get_harmonized_qtl_data(parquet_file)
    return innerjoin(
        variant_ids,
        select(sig_cis_eqtls, qtls_columns...), 
        on = :UNIV_ID
    )
end

function get_GTEx_table_matching_variants(variant_ids, parquet_files; qtl_type=:EQTL)
    qtls_columns = PostGWAS.GTEx_columns_selection(;qtl_type=qtl_type)
    # Get all TISSUE/TARGET/VARIANT association matching any of the variant in the locus
    tasks = map(parquet_files) do parquet_file
        Threads.@spawn PostGWAS.get_qtl_data_matching_locus(variant_ids, parquet_file, qtls_columns)
    end
    variant_tissue_target_qtl_data = fetch.(tasks)
    return vcat(variant_tissue_target_qtl_data...)
end

function make_GTEx_table!(db, variant_ids, qtl_type; gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data"))
    qtl_file = if qtl_type == :EQTL
        "GTEx_Analysis_v10_eQTL_updated"
    elseif qtl_type ==:SQTL
        "GTEx_Analysis_v10_sQTL_updated"
    else
        throw(ArgumentError(string("Unknown qtl_type:", qtl_type)))
    end
    # Get QTL info
    qtl_dir = joinpath(gtex_cache_dir, qtl_file)
    parquet_files = filter(endswith(".parquet"), readdir(qtl_dir, join=true))
    QTL_table = PostGWAS.get_GTEx_table_matching_variants(variant_ids, parquet_files; qtl_type=qtl_type)
    # Make Table
    table_name = string("GTEX_", qtl_type)
    SQLite.load!(QTL_table, db, table_name)
    index_name = string(table_name, "_UNIV_ID_INDEX")
    SQLite.createindex!(db, table_name, index_name, "UNIV_ID", unique=false)
end

function make_gtex_tables!(db; gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data"))
    loci_variants = all_loci_variant_ids(db)
    @info "Creating GTEx eQTL table..."
    make_GTEx_table!(db, loci_variants, :EQTL; gtex_cache_dir=gtex_cache_dir)
    @info "Creating GTEx sQTL table..."
    make_GTEx_table!(db, loci_variants, :SQTL; gtex_cache_dir=gtex_cache_dir)
    return 0
end