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


function get_locus_table(db, locus_id)
    DBInterface.execute(
            db,
            """SELECT 
                    GWAS_RESULTS.CHROM, 
                    GWAS_RESULTS.POS, 
                    GWAS_RESULTS.UNIV_ID AS ID, 
                    GWAS_RESULTS.REF, 
                    GWAS_RESULTS.ALT, 
                    GWAS_RESULTS.ALT_FREQ,
                    GWAS_RESULTS.LOG10P,
                    GWAS_RESULTS.BETA,
                    GWAS_RESULTS.SE,
                    MOST_SEVERE_CONSEQUENCES.MOST_SEVERE_CONSEQUENCE,
                    FINEMAPPING_RESULTS.LOCUS_ID,
                    FINEMAPPING_RESULTS.PIP, 
                    FINEMAPPING_RESULTS.CS,
                    FINEMAPPING_RESULTS.PHASED_R2
                FROM GWAS_RESULTS 
                INNER JOIN FINEMAPPING_RESULTS ON GWAS_RESULTS.UNIV_ID = FINEMAPPING_RESULTS.UNIV_ID
                INNER JOIN MOST_SEVERE_CONSEQUENCES ON GWAS_RESULTS.UNIV_ID = MOST_SEVERE_CONSEQUENCES.UNIV_ID
                WHERE FINEMAPPING_RESULTS.LOCUS_ID = '$locus_id'
            """
    ) |> DataFrame
end

function get_ensembl_transcript_consequences(db, variant_id)
    DBInterface.execute(
                db,
                """SELECT 
                        *
                    FROM TRANSCRIPT_CONSEQUENCES 
                    WHERE TRANSCRIPT_CONSEQUENCES.UNIV_ID = '$variant_id'
                """
    ) |> DataFrame
end

function get_GTEx_qtl_info(db, variant_id; table="GTEX_EQTL")
    DBInterface.execute(
                db,
                """SELECT 
                        *
                    FROM $table
                    WHERE $table.UNIV_ID = '$variant_id'
                """
    ) |> DataFrame
end