
# Load all required modules for the job
module load bcftools/1.10

# Script
# ------

bcftools +split-vep path/to/input/input_vcf_file.vcf.gz -f \
'%CHROM:%POS %ID %REF %ALT %QUAL %TYPE [%GT] [%AD] [%DP] [%PL] [%GQ] %CSQ\n' -A tab \
> path/to/output/output_tab–delimited_text_file.txt
