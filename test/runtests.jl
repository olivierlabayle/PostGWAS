using Test
using PostGWAS

vep_cache_dir = joinpath(ENV["HOME"], "vep_data")
gtex_cache_dir = joinpath(ENV["HOME"], "gtex_data")
phenotype = "SEVERE_COVID_19"
gwas_file = "test/assets/META_ANALYSIS.SEVERE_COVID_19.gwas.tsv"
finemapping_file = "test/assets/META_ANALYSIS.SEVERE_COVID_19.finemapping.tsv"
gtex_dir = "../assets/GTEx_Analysis_v10_eQTL_updated"
confounders = ["PC1", "PC2", "PC3", "PC4", "PC5"]
outcome_extra_covariates = ["AGE", "SEX"]
positivity_constraint = 0.01
maf_threshold = 0.005
pip_threshold = 0.
workdir = "../assets/"

