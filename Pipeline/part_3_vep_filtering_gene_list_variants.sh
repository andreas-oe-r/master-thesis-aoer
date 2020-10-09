
# Load all required modules for the job
module load perl/5.24.0
module load htslib/1.10.2
module load bcftools/1.9
module load bedtools/2.27.1
module load vt/0.5772
module load ensembl-vep/100.1

# Script
# ------


# Version 5 – IN USE
filter_vep \
--input_file path/to/input/input_vcf_file.vcf.gz \
--gz \
--format vcf \
--filter \
	"Consequence is not synonymous_variant \
	and Consequence is not intron_variant \
	and Consequence match not intergenic \
	and Consequence match not intragenic \
	and Consequence is not upstream_gene_variant \
	and Consequence is not downstream_gene_variant \
	and BIOTYPE match not pseudogene" \
--only_matched | 
filter_vep \
--output_file path/to/output/output_vcf_file.vcf.gz \
--filter \
	"(gnomADg_AF < 0.01 or not gnomADg_AF) \
	and (gnomAD_AF < 0.01 or not gnomAD_AF)" \
--only_matched \
--force_overwrite
