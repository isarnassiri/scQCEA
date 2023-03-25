#!/bin/bash

#$ -cwd
#$ -N seurat -j y
#$ -P htseq.prjc.low
#$ -q htseq.qd@@htseq.hgd
#$ -pe shmem 4

module purge
module add R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

# ## Some Paths
SCRIPTS_PATH="/well/singlecell/scripts"
CONDA_PATH="/well/singlecell/scripts/"
ENV_NAME="seurat312"

# Job Arguments
SAMPLES_FILE=$1

# Task Arguments
PROJECT_ID=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $1}' $SAMPLES_FILE)
GEX_LIBRARY_ID=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $2}' $SAMPLES_FILE)
GEX_SEQUENCING_ID=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $3}' $SAMPLES_FILE)
GEX_DATA_DIR=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $4}' $SAMPLES_FILE)
HTO_LIBRARY_ID=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $5}' $SAMPLES_FILE)
HTO_SEQUENCING_ID=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $6}' $SAMPLES_FILE)
HTO_DATA_DIR=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $7}' $SAMPLES_FILE)
ADT_LIBRARY_ID=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $8}' $SAMPLES_FILE)
ADT_SEQUENCING_ID=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $9}' $SAMPLES_FILE)
ADT_DATA_DIR=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $10}' $SAMPLES_FILE)
OUTPUT_DIR=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $11}' $SAMPLES_FILE)
DEMUX_METHOD=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $12}' $SAMPLES_FILE)
DEMUX_MS_AT=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $13}' $SAMPLES_FILE)
SAVE_RAW=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $14}' $SAMPLES_FILE)
SAVE_FILTERED=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $15}' $SAMPLES_FILE)
SAVE_SPLITTED=$(awk -F '\t' -v "line=$SGE_TASK_ID" 'NR==line {print $16}' $SAMPLES_FILE)

echo "********************************************************"
echo "* Job Details"
echo "********************************************************"
echo "SGE job ID       : "$JOB_ID
echo "SGE task ID      : "$SGE_TASK_ID
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
echo "GEX Sequencing ID   : ${GEX_SEQUENCING_ID}"
echo "GEX filt. mtx path  : ${GEX_DATA_DIR}"
echo "HTO Library ID      : ${HTO_LIBRARY_ID}"
echo "HTO Sequencing ID   : ${HTO_SEQUENCING_ID}"
echo "HTO raw mtx path    : ${HTO_DATA_DIR}"
echo "ADT Library ID      : ${ADT_LIBRARY_ID}"
echo "ADT Sequencing ID   : ${ADT_SEQUENCING_ID}"
echo "ADT raw mtx path    : ${ADT_DATA_DIR}"
echo "Output path         : ${OUTPUT_DIR}"
echo "Demux method        : ${DEMUX_METHOD}"
echo "MULTIseq AutoThresh : ${DEMUX_MS_AT}"
echo "Save raw mtx        : ${SAVE_RAW}"
echo "Save filtered mtx   : ${SAVE_FILTERED}"
echo "Save splitted mtx   : ${SAVE_SPLITTED}"
echo
echo "********************************************************"
echo "* Software Paths"
echo "********************************************************"
echo "Scripts Path        : ${SCRIPTS_PATH}"
echo "Conda Path          : ${CONDA_PATH}"
echo

[[ ! -z "${OUTPUT_DIR}" ]] && OUTPUT_DIR="--output_dir ${OUTPUT_DIR}"
[[ ! -z "${PROJECT_ID}" ]] && PROJECT_ID="--project_id ${PROJECT_ID}"
[[ ! -z "${GEX_LIBRARY_ID}" ]] && GEX_LIBRARY_ID="--gex_library_id ${GEX_LIBRARY_ID}"
[[ ! -z "${HTO_LIBRARY_ID}" ]] && HTO_LIBRARY_ID="--hto_library_id ${HTO_LIBRARY_ID}"
[[ ! -z "${ADT_LIBRARY_ID}" ]] && ADT_LIBRARY_ID="--adt_library_id ${ADT_LIBRARY_ID}"
[[ ! -z "${GEX_SEQUENCING_ID}" ]] && GEX_SEQUENCING_ID="--gex_sequencing_id ${GEX_SEQUENCING_ID}"
[[ ! -z "${HTO_SEQUENCING_ID}" ]] && HTO_SEQUENCING_ID="--hto_sequencing_id ${HTO_SEQUENCING_ID}"
[[ ! -z "${ADT_SEQUENCING_ID}" ]] && ADT_SEQUENCING_ID="--adt_sequencing_id ${ADT_SEQUENCING_ID}"
[[ ! -z "${GEX_DATA_DIR}" ]] && GEX_DATA_DIR="--gex_data_dir ${GEX_DATA_DIR}"
[[ ! -z "${HTO_DATA_DIR}" ]] && HTO_DATA_DIR="--hto_data_dir ${HTO_DATA_DIR}"
[[ ! -z "${ADT_DATA_DIR}" ]] && ADT_DATA_DIR="--adt_data_dir ${ADT_DATA_DIR}"
[[ ! -z "${DEMUX_METHOD}" ]] && DEMUX_METHOD="--demux_method ${DEMUX_METHOD}"
[[ ! -z "${DEMUX_MS_AT}" && "${DEMUX_MS_AT^^}" =~ ^(0|F|FALSE)$ ]] && DEMUX_MS_AT="--demux_ms_autothresh FALSE" || DEMUX_MS_AT="--demux_ms_autothresh TRUE"
[[ ! -z "${SAVE_RAW}" && "${SAVE_RAW^^}" =~ ^(0|F|FALSE)$ ]] && SAVE_RAW="--save_full_matrix FALSE" || SAVE_RAW="--save_full_matrix TRUE"
[[ ! -z "${SAVE_FILTERED}" && "${SAVE_FILTERED^^}" =~ ^(0|F|FALSE)$ ]] && SAVE_FILTERED="--save_filtered_matrix FALSE" || SAVE_FILTERED="--save_filtered_matrix TRUE"
[[ ! -z "${SAVE_SPLITTED}" && ! "${SAVE_SPLITTED^^}" =~ ^(0|F|FALSE)$ ]] && SAVE_SPLITTED="--save_splitted_matrices TRUE" || SAVE_SPLITTED="--save_splitted_matrices FALSE"

echo "********************************************************"
echo "["`date`"] Activating CONDA environment: ${ENV_NAME}"
echo "********************************************************"
echo

echo "********************************************************"
echo "["`date`"] Running tool"
echo "********************************************************"
echo
(
set -x
Rscript ${SCRIPTS_PATH}/run.qsub.process_feature_sample.R \
${SAVE_RAW} ${SAVE_FILTERED} ${SAVE_SPLITTED} \
${OUTPUT_DIR} ${PROJECT_ID} ${DEMUX_METHOD} ${DEMUX_MS_AT} \
${GEX_LIBRARY_ID} ${GEX_SEQUENCING_ID} ${GEX_DATA_DIR} \
${HTO_LIBRARY_ID} ${HTO_SEQUENCING_ID} ${HTO_DATA_DIR} \
${ADT_LIBRARY_ID} ${ADT_SEQUENCING_ID} ${ADT_DATA_DIR}
)

echo "********************************************************"
echo "["`date`"] De-activating CONDA environment"
echo "********************************************************"

echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0

