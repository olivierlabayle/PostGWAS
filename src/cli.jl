function cli_settings()
    s = ArgParseSettings(
        description="PostGWAS CLI",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(PostGWAS))
    )

    @add_arg_table! s begin
        "enrich-loci"
            action = :command
            help = "Runs meta-analysis across GWAS results."
    end

    @add_arg_table! s["enrich-loci"] begin
        "gwas-file"
            arg_type = String
            required = true
            help = "List of REGENIE results files to be meta-analysed."
    
        "fine-mapping-file"
            arg_type = String
            help = "List of groups to exclude from the meta-analysis, comma separated."

        "--output-file"
            arg_type = String
            help = "rich_loci.jld2."
            default = "STDERR"

        "--vep-cache-dir"
            arg_type = String
            help = "Prefix to output files."
            default = joinpath(ENV["HOME"], "vep_data")
        
        "--gtex-cache-dir"
            arg_type = String
            help = "Prefix to output files."
            default = joinpath(ENV["HOME"], "gtex_data")
    end

    return s
end