#!/bin/bash
#SBATCH -p short
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir=./
#SBATCH --job-name=Enrichment
#SBATCH --output=Enrichment-%A_%a.out

umask 002
module purge
module add R/4.1.0-foss-2021a
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

# ## Some Paths
SCRIPTS_PATH="/gpfs3/well/singlecell/scripts/"

# Job Arguments
SAMPLES_FILE=$1

# Task Arguments
PROJECT_ID=$(awk -F '\t' -v "line=$SLURM_ARRAY_TASK_ID" 'NR==line {print $1}' $SAMPLES_FILE)
GEX_LIBRARY_ID=$(awk -F '\t' -v "line=$SLURM_ARRAY_TASK_ID" 'NR==line {print $2}' $SAMPLES_FILE)
GEX_DATA_DIR=$(awk -F '\t' -v "line=$SLURM_ARRAY_TASK_ID" 'NR==line {print $3}' $SAMPLES_FILE)
BACKEND_DIR=$(awk -F '\t' -v "line=$SLURM_ARRAY_TASK_ID" 'NR==line {print $4}' $SAMPLES_FILE)
ORGANISM=$(awk -F '\t' -v "line=$SLURM_ARRAY_TASK_ID" 'NR==line {print $5}' $SAMPLES_FILE)

echo "********************************************************"
echo "* Job Details"
echo "********************************************************"
echo "SGE job ID       : "$SLURM_ARRAY_JOB_ID
echo "SGE task ID      : "$SLURM_ARRAY_TASK_ID
echo "Run on host      : "`hostname`
echo "Operating system : "`uname -s`
echo "Username         : "`whoami`
echo "Started at       : "`date`
echo
echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "Samples File     : ${SAMPLES_FILE}"
echo
echo "********************************************************"
echo "* Task Parameters"
echo "********************************************************"
echo "Project ID          : ${PROJECT_ID}"
echo "GEX Library ID      : ${GEX_LIBRARY_ID}"
echo "Input path         : ${GEX_DATA_DIR}"
echo "Backend data path         : ${BACKEND_DIR}"
echo "Organism        : ${ORGANISM}"
echo
echo "********************************************************"
echo "* Software Paths"
echo "********************************************************"
echo "Scripts Path        : ${SCRIPTS_PATH}"
echo

# opt$project_id <- "P210341"
# opt$gex_library_id <- "AYO9038A1"
# opt$input.dir <- "/well/singlecell/P210341/10X-gex-grouped/AYO9038A1"
# opt$backend.data.dir <- "/well/singlecell/references/reference_gene_sets/human"
# opt$organism <- "hsapiens" #"mmusculus"

[[ ! -z "${PROJECT_ID}" ]] && PROJECT_ID="--project_id ${PROJECT_ID}"
[[ ! -z "${GEX_LIBRARY_ID}" ]] && GEX_LIBRARY_ID="--gex_library_id ${GEX_LIBRARY_ID}"
[[ ! -z "${GEX_DATA_DIR}" ]] && GEX_DATA_DIR="--input.dir ${GEX_DATA_DIR}"
[[ ! -z "${BACKEND_DIR}" ]] && BACKEND_DIR="--backend.data.dir ${BACKEND_DIR}"
[[ ! -z "${ORGANISM}" ]] && ORGANISM="--organism ${ORGANISM}"

echo "********************************************************"
echo "["`date`"] Running tool"
echo "********************************************************"
echo
(
set -x
Rscript ${SCRIPTS_PATH}/run.qsub.scRNAseq_enrichment_analysis.R \
${PROJECT_ID} \
${GEX_LIBRARY_ID} \
${GEX_DATA_DIR} \
${BACKEND_DIR} \
${ORGANISM}
)

echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
