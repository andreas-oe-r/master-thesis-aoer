
# Load all required modules for the job
module load perl/5.24.0
module load htslib/1.10.2
module load bcftools/1.9
module load bedtools/2.27.1
module load vt/0.5772
module load ensembl-vep/100.1

# Script
# ------

vep \
--input_file path/to/input/input_vcf_file.vcf.gz \
--fork 14 \
--offline \
--cache \
--config path/to/config_file.txt \
--format vcf \
--everything \
--no_escape \
--failed 1 \
--distance 500 \
--stats_text \
--stats_file STATS.tsv \
--output_file path/to/output/output_vcf_file.vcf.gz \
--vcf \
--compress_output bgzip \
--force

