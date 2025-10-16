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