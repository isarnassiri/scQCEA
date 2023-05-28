#!/bin/bash
#$ -cwd
umask 002

SCRIPTS_PATH=$(pwd)

#----- Font Color
RED='\033[1;31m' # 1 means bold
NC='\033[0m'     # No Color

function red {
    printf "${bold}${RED}$@${NC}${normal}\n" 
}

#--------------------------------------------

while getopts ":i:c:g:v:e:p:o:f:" opt; do
  case $opt in
    i) PID="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1;;
  esac
done

#---
types=('vdj' 'gex' 'feat' 'arc' 'atac' 'vdj-grouped' 'gex-grouped' 'feat-grouped' 'arc-grouped' 'atac-grouped')

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME
cd ${SCRIPTS_PATH}

for TYPE in "${types[@]}"
do

FOLDER='10X-'$TYPE

if [ -d "$FOLDER" ]; then

echo "$FOLDER:"

mkdir -p ${SCRIPTS_PATH}/Inputs/10X-$TYPE/

cd 10X-$TYPE/
tar czf ${SCRIPTS_PATH}/Inputs/10X-$TYPE/metrics_summary_$TYPE'_'$NAMEFOLDER.tar.gz */outs/*summary.csv

cd ${SCRIPTS_PATH}/Inputs/10X-$TYPE/
tar -xvzf metrics_summary_$TYPE'_'$NAMEFOLDER.tar.gz
rm 'metrics_summary_'$TYPE'_'$NAMEFOLDER'.tar.gz'

cd ${SCRIPTS_PATH}/

fi

done

#--- 
if [ -d "10X-arc" ]; then
cp ${SCRIPTS_PATH}/samples.arc ${SCRIPTS_PATH}/Inputs/10X-arc/
fi

#--- copy samples.metadata
cp ${SCRIPTS_PATH}/samples.metadata ${SCRIPTS_PATH}/Inputs/

done

#---  
for NAME in "${PID[@]}"
do
NAMEFOLDER=$NAME
cd ${SCRIPTS_PATH}
if [ -d "10X-aggregate" ]; then
cd '10X-aggregate'

mkdir ${SCRIPTS_PATH}/Inputs/10X-aggregate/

tar czf ${SCRIPTS_PATH}/Inputs/10X-aggregate/metadata_$NAMEFOLDER.tar.gz */objects/raw_matrix/metadata.tsv.gz

cd ${SCRIPTS_PATH}/Inputs/10X-aggregate/

tar -xvzf metadata_$NAMEFOLDER.tar.gz
rm metadata_$NAMEFOLDER.tar.gz

fi
done

echo $(red "Ready for next steps.")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fi
