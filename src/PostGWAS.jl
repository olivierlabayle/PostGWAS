module PostGWAS

using CSV
using DataFrames
using Parquet
using TMLE
using Combinatorics
using RCall
using Downloads
using JSON
using Base.Threads
using JLD2
using ArgParse

include("utils.jl")
include("ensembl.jl")
include("gtex.jl")
include("cli.jl")

export main

function main(ARGS)
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    @info "Running PostGWAS CLI: $cmd"
    cmd_settings = settings[cmd]
    if cmd == "enrich-loci"
        make_rich_loci_dataset(
            cmd_settings["gwas-file"], 
            cmd_settings["fine-mapping-file"];
            output_file=cmd_settings["output-file"],
            vep_cache_dir=cmd_settings["vep-cache-dir"],
            gtex_cache_dir=cmd_settings["gtex-cache-dir"]
        )
    else
        throw(ArgumentError("Unknown command: $cmd"))
    end

    return 0
end

function make_rich_loci_dataset(gwas_file, finemapping_file; 
    output_file="rich_loci.jld2",
    vep_cache_dir=joinpath(ENV["HOME"], "vep_data"),
    gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data")
    )
    # Download VEP cache data if not already done
    download_ensembl_vep_cache(vep_cache_dir=vep_cache_dir)
    # Download GTEx QTL data if not already done
    PostGWAS.download_all_GTEx_QTLs_v10(;gtex_cache_dir=gtex_cache_dir)
    # Load GWAS results
    gwas_results = PostGWAS.get_filtered_gwas_results(gwas_file)
    # For each finemapped locus:
    # - add GWAS results
    # - add Ensembl annotations
    # - map to GTEx eQTLs and sQTLs
    loci = PostGWAS.get_finemapped_loci(finemapping_file)
    @info "$(length(loci)) loci to be enriched."
    for (locus_key, locus) in pairs(loci) 
        # locus_key, locus = first(pairs(loci))
        @info "Processing locus: $(locus_key.LOCUS_ID)"
        locus = DataFrame(locus)
        PostGWAS.add_gwas_info!(locus, gwas_results)
        PostGWAS.add_ensembl_annotations!(locus; vep_cache_dir=vep_cache_dir)
        PostGWAS.add_GTEX_info!(locus; gtex_cache_dir=gtex_cache_dir)
        jldopen(io -> io[locus_key.LOCUS_ID] = locus, output_file, "a+")
    end
    @info "Done."
    return 0
end

end
