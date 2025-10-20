function download_ensembl_vep_cache(;vep_cache_dir=joinpath(ENV["HOME"], "vep_data"))
    @info "Installing Ensembl VEP and cache in $vep_cache_dir, this may take a while..."
    run(`docker run --rm -it --platform linux/amd64 -v $vep_cache_dir:/data ensemblorg/ensembl-vep INSTALL.pl -a cf -s homo_sapiens -y GRCh38`)
end

function recover_universal_IDs_from_ensembl_ann(ensembl_ann)
    chr, pos, varid, ref, alt, _ = split(ensembl_ann["input"], "\t")
    return string(chr, ":", pos, ":", ref, ":", alt)
end

function add_ensembl_annotations!(
    ensembl_annotation_tables,
    locus; 
    vep_cache_dir=joinpath(ENV["HOME"], "vep_data")
    )
    # Make ENSEMBL VEP Inputs
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
    # Run ENSEMBL VEP
    run(
        `docker run -v $vep_cache_dir:/data -v $tmpdir:/mnt/in_out ensemblorg/ensembl-vep \
            vep --cache --offline --format vcf --json --force_overwrite \
            --input_file /mnt/in_out/input.vcf \
            --output_file /mnt/in_out/output.json \
            --regulatory`
    )
    # Parse ENSEMBL VEP Results into DataFrames
    ensembl_anns = JSON.parse.(readlines(joinpath(tmpdir, "output.json")))
    colnames_dict = Dict(
        "transcript_consequences" => Dict(key => missing for key in ["amino_acids",
            "variant_allele",
            "gene_id",
            "consequence_terms",
            "distance",
            "codons",
            "protein_start",
            "cds_start",
            "protein_end",
            "cdna_start",
            "strand",
            "transcript_id",
            "cds_end",
            "impact",
            "flags",
            "cdna_end"
        ]),
        "regulatory_feature_consequences" => Dict(key => missing for key in [
            "impact",
            "biotype",
            "variant_allele",
            "regulatory_feature_id",
            "consequence_terms"
        ]),
        "intergenic_consequences" => Dict(key => missing for key in [
            "impact",
            "variant_allele",
            "consequence_terms"
        ])
    )
    ann_keys = keys(colnames_dict)
    locus_annotation_dict = Dict(key => [] for key in ann_keys)
    for ensembl_ann in ensembl_anns
        for ann_key in ann_keys
            if haskey(ensembl_ann, ann_key)
                homogeneous_dicts = map(ensembl_ann[ann_key]) do ann_dict
                    ann_dict["univ_id"] = ensembl_ann["id"]
                    merge(colnames_dict[ann_key], ann_dict)
                end
                push!(
                    locus_annotation_dict[ann_key],
                    DataFrame(homogeneous_dicts)
                )
            end
        end
    end
    locus_annotation_dict = Dict(ann_key => vcat(vals...) for (ann_key, vals) in locus_annotation_dict)
    locus_annotation_dict["most_severe_consequences"] = DataFrame([(univ_id = row["id"], most_severe_consequence = row["most_severe_consequence"]) for row in ensembl_anns])
    # Update ensembl_annotation_tables with the locus specific info
    for (ann_key, ann_table) in ensembl_annotation_tables
        append!(ann_table, locus_annotation_dict[ann_key])
    end
    # CLeanup
    rm(tmpdir; force=true, recursive=true)
end


function make_ensembl_tables!(db; vep_cache_dir=joinpath(ENV["HOME"], "vep_data"))
    @info "Creating ENSEMBL VEP annotation tables..."
    # Initialise tables
    ensembl_annotation_tables = Dict(key => DataFrame() 
        for key in (
            "transcript_consequences", 
            "regulatory_feature_consequences", 
            "intergenic_consequences", 
            "most_severe_consequences"
        )
    )
    # Update tables with each locus using the ENSEMBL VEP
    locus_ids = get_loci_ids(db)
    for (locus_index, locus_id) in enumerate(locus_ids)
        @info "Annotating locus: $locus_id ($locus_index/$(length(locus_ids)))"
        locus = PostGWAS.get_locus_from_id(db, locus_id)
        PostGWAS.add_ensembl_annotations!(ensembl_annotation_tables, locus; vep_cache_dir=vep_cache_dir)
    end

    # Finalising the tables
    for (ann_key, ann_table) in ensembl_annotation_tables
        ## Make uppercase column names
        rename!(uppercase, ann_table)
        # Replace vectors by strings
        for colname in names(ann_table)
            if eltype(ann_table[!, colname]) <: AbstractVector
                ann_table[!, colname] = join.(ann_table[!, colname], ";")
            end
        end
        # Make table and index
        table_name = uppercase(ann_key)
        SQLite.load!(ann_table, db, table_name)
        index_name = string(table_name, "UNIV_ID_INDEX")
        SQLite.createindex!(db, table_name, index_name, "UNIV_ID", unique=false)
    end
    return 0
end