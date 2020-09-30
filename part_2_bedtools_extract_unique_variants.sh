
# Load all required modules for the job
module load bedtools/2.28.0

# Script
# ------

bedtools intersect -wa -header \
-a path/to/input/input_vcf_file.vcf.gz \
-b path/to/bed_extraction_file.bed \
> path/to/output/output_vcf_file.vcf

gzip path/to/output_vcf_file.vcf.gz
