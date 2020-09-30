
## Load all required modules for the job
module load java/1.7.0
module load igv/2.4.9

# Script
# ------

xvfb-run -d path/to/igv_executable.sh -b path/to/patient/igv_batch_file.bat


