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
using SQLite

include("utils.jl")
include("ensembl.jl")
include("gtex.jl")
include("cli.jl")
include("db_api.jl")

export main

function main(ARGS)
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    @info "Running PostGWAS CLI: $cmd"
    cmd_settings = settings[cmd]
    if cmd == "enrich-loci"
        make_db(
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

function make_db(gwas_file, finemapping_file; 
    output_file="postgwas.s3db",
    vep_cache_dir=joinpath(ENV["HOME"], "vep_data"),
    gtex_cache_dir=joinpath(ENV["HOME"], "gtex_data")
    )
    # Download VEP cache data if not already done
    download_ensembl_vep_cache(vep_cache_dir=vep_cache_dir)
    # Download GTEx QTL data if not already done
    PostGWAS.download_all_GTEx_QTLs_v10(;gtex_cache_dir=gtex_cache_dir)
    # Make DB
    isfile(output_file) && rm(output_file)
    db = SQLite.DB(output_file)

    # Make GWAS results table
    PostGWAS.make_gwas_table!(db, gwas_file)
    # Make Fine-Mapping results table
    PostGWAS.make_finemapping_table!(db, finemapping_file)
    # Make ENSEMBL annottations tables
    PostGWAS.make_ensembl_tables!(db; vep_cache_dir=vep_cache_dir)
    # Make GTEx annottations tables
    PostGWAS.make_gtex_tables!(db; gtex_cache_dir=gtex_cache_dir)

    @info "Done."
    return 0
end

end
