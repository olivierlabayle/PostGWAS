join_variant_info(chr_col, pos_col, ref_col, alt_col) = 
    string.(string.(chr_col), ":", string.(pos_col), ":", string.(ref_col), ":", string.(alt_col))

function add_universal_ID!(df; cols=[:CHROM, :POS, :REF, :ALT])
    transform!(df, 
        cols => join_variant_info => :UNIV_ID)
    return df
end

function make_gwas_table!(db, gwas_file)
    gwas_results = CSV.read(gwas_file, DataFrame; select=["CHROM", "ID", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"])
    rename!(gwas_results, "GENPOS" => "POS", "ALLELE0" => "REF", "ALLELE1" => "ALT", "A1FREQ" => "ALT_FREQ")
    PostGWAS.add_universal_ID!(gwas_results, cols=[:CHROM, :POS, :REF, :ALT])
    select!(gwas_results, :UNIV_ID, :CHROM, :POS, :REF, :ALT, :BETA, :SE, :LOG10P, :ALT_FREQ, :ID)
    SQLite.load!(gwas_results, db, "GWAS_RESULTS")
    SQLite.createindex!(db, "GWAS_RESULTS", "GWAS_RESULTS_UNIV_ID_INDEX", "UNIV_ID")
end

function get_harmonized_finemapping_results(finemapping_file)
    finemapping_results = CSV.read(finemapping_file, DataFrame, delim="\t")
    rename!(finemapping_results, "#CHROM" => "CHR")
    add_universal_ID!(finemapping_results)
    return finemapping_results
end

function make_finemapping_table!(db, finemapping_file)
    finemapping_df = CSV.read(finemapping_file, DataFrame)
    PostGWAS.add_universal_ID!(finemapping_df)
    select!(finemapping_df, :LOCUS_ID, :UNIV_ID, :CHROM, :POS, :REF, :ALT, :ID, :PIP, :CS, :PHASED_R2)
    SQLite.load!(finemapping_df, db, "FINEMAPPING_RESULTS")
    SQLite.createindex!(db, "FINEMAPPING_RESULTS", "FINEMAPPING_RESULTS_LOCUS_ID_INDEX", "LOCUS_ID", unique=false)
    SQLite.createindex!(db, "FINEMAPPING_RESULTS", "FINEMAPPING_RESULTS_LOCUS_ID_UNIV_ID_INDEX", ["LOCUS_ID", "UNIV_ID"])
end

function get_loci_ids(db)
    return DataFrame(DBInterface.execute(db, "SELECT DISTINCT LOCUS_ID FROM FINEMAPPING_RESULTS")).LOCUS_ID
end

function get_locus_from_id(db, locus_id)
    DBInterface.execute(
        db,
        "SELECT CHROM, POS, UNIV_ID AS ID, REF, ALT FROM FINEMAPPING_RESULTS WHERE LOCUS_ID = '$locus_id'"
    ) |> DataFrame
end

add_gwas_info!(locus, gwas_results) = leftjoin!(locus, gwas_results, on = :UNIV_ID)