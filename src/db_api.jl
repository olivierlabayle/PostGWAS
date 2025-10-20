function get_loci_ids(db)
    return DataFrame(DBInterface.execute(db, "SELECT DISTINCT LOCUS_ID FROM FINEMAPPING_RESULTS")).LOCUS_ID
end

function get_locus_from_id(db, locus_id)
    DBInterface.execute(
        db,
        "SELECT CHROM, POS, UNIV_ID, REF, ALT FROM FINEMAPPING_RESULTS WHERE LOCUS_ID = '$locus_id'"
    ) |> DataFrame
end


function all_loci_variant_ids(db)
    DBInterface.execute(
            db,
            "SELECT DISTINCT UNIV_ID FROM FINEMAPPING_RESULTS"
    ) |> DataFrame
end