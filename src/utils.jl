join_variant_info(chr_col, pos_col, ref_col, alt_col) = 
    string.(string.(chr_col), ":", string.(pos_col), ":", string.(ref_col), ":", string.(alt_col))

function add_universal_ID!(df; cols=[:CHROM, :POS, :REF, :ALT])
    transform!(df, 
        cols => join_variant_info => :UNIV_ID)
    return df
end

function get_filtered_gwas_results(gwas_file)
    gwas_results = CSV.read(gwas_file, DataFrame; select=["CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"])
    add_universal_ID!(gwas_results, cols=[:CHROM, :GENPOS, :ALLELE0, :ALLELE1])
    return select!(gwas_results, :UNIV_ID, :BETA, :SE, :LOG10P, :A1FREQ => :ALT_FREQ)
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

add_gwas_info!(locus, gwas_results) = leftjoin!(locus, gwas_results, on = :UNIV_ID)