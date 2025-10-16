using Test
using PostGWAS

vep_cache_dir = joinpath(ENV["HOME"], "vep_data")
gtex_cache_dir = joinpath(ENV["HOME"], "gtex_data")
gwas_file = "test/assets/META_ANALYSIS.SEVERE_COVID_19.gwas.tsv"
finemapping_file = "test/assets/META_ANALYSIS.SEVERE_COVID_19.finemapping.tsv"
output_file="rich_loci.jld2"



PostGWAS.make_rich_loci_dataset(gwas_file, finemapping_file)