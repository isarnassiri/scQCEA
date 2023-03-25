#!/bin/bash
#$ -cwd
umask 002

#----- Font Color
RED='\033[1;31m' # 1 means bold
NC='\033[0m'     # No Color

function red {
    printf "${bold}${RED}$@${NC}${normal}\n" 
}

#--------------------------------------------
echo $(red "If you want to run a specific module, rerun the script and set the following parameters after project ID:")
echo "-i Project ID (e.g., P200231)"
echo "-c TRUE to Check Log Files"
echo "-g TRUE to Generate Aggregate File"
echo "-v TRUE to Generate CSV"
echo "-e TRUE to run Enrichment"
echo "-p TRUE to Convert PDF"
echo "-o TRUE to Orgenize Files Folders"
echo "-f TRUE to Generate an FTP link"

while getopts ":i:c:g:v:e:p:o:f:" opt; do
  case $opt in
    i) PID="$OPTARG";;
    c) CheckLogFiles="$OPTARG";;
    g) GenerateAggregateFile="$OPTARG";;
    v) GenerateCSV="$OPTARG";;
    e) Enrichment="$OPTARG";;
    p) ConvertPDF="$OPTARG";;
    o) OrgenizeFilesFolders="$OPTARG";;
    f) FTPLinks="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1;;
  esac
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running Selected Modules
if
[ ! -z $CheckLogFiles ] || 
[ ! -z $GenerateAggregateFile ] ||
[ ! -z $GenerateCSV ] ||
[ ! -z $Enrichment ] ||
[ ! -z $ConvertPDF ] ||
[ ! -z $OrgenizeFilesFolders ] ||
[ ! -z $FTPLinks ]
then
echo $(red "Running Selected Modules")

#-------------------------------------------- CheckLogFiles module
if [ ! -z $CheckLogFiles ] 
then 

printf "Argument CheckLogFiles is %s\n" "$CheckLogFiles"

file_name=${PID[@]}
current_time=$(date "+%Y.%m.%d-%H.%M")
log_fileName=/well/singlecell/_logs_/$file_name.$current_time.log
echo "Log file name: " "$log_fileName"

types=('vdj' 'gex' 'feat' 'arc' 'vdj-grouped' 'gex-grouped' 'feat-grouped')

for NAME in "${PID[@]}"
do

echo $(red "${PID[@]}")
echo ${PID[@]} >> $log_fileName

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER

if [ -d "Inputs" ]; then
echo $(red "Inputs/ folder exists; renamed.")
mv Inputs Inputs.$current_time
fi

for TYPE in "${types[@]}"
do

FOLDER='10X-'$TYPE

if [ -d "$FOLDER" ]; then

if echo $TYPE | grep --quiet 'grouped'; then
     TYPE_revised=${TYPE%"-grouped"}.grouped
     NumberExpectedFolders=$(cat samples.$TYPE_revised | wc -l)
     else
     NumberExpectedFolders=$(cat samples.$TYPE | wc -l)
fi

cd $FOLDER

# count "find" results
NumberLogFiles=$(find -maxdepth 1 -name '*.log' -exec printf %c {} + | wc -c)
NumberLogFolders=$(ls -l | grep ^d | wc -l)

count=0

echo $FOLDER >> $log_fileName

for f in $(find -maxdepth 1 -name '*.log')
do
     
     #grep -rnw $f -e 'Pipestance completed successfully!'
     
     if grep --quiet $f -e 'Pipestance completed successfully!'; then
     ((count++))
     else
      echo $(red "$f")
      echo $f >> $log_fileName

     fi
done
 
if [ "$NumberExpectedFolders" == "$NumberLogFolders" ]; then
    echo $(red "$FOLDER: MetaData and Number of folders matched.");
    else
    echo $(red "$FOLDER: MetaData and Number of folders did NOT match.");
fi

if [ "$NumberLogFiles" == "$NumberLogFolders" ] && [ "$NumberLogFiles" != "$count" ]; then
    echo $(red "$FOLDER: $((NumberLogFiles-count)) failed sample(s).")
fi

if [ "$NumberLogFiles" != "$NumberLogFolders" ] && [ "$NumberLogFiles" != "$count" ]; then
     echo $(red "$FOLDER: $((NumberLogFiles-count)) failed sample(s) and $((NumberLogFolders-NumberLogFiles)) missed log files." )
fi

if [ "$NumberLogFiles" != "$NumberLogFolders" ] && [ "$NumberLogFiles" == "$count" ]; then
    echo $(red "$FOLDER: $((NumberLogFolders-NumberLogFiles)) missed log files." )
fi

if [ "$NumberLogFiles" == "$NumberLogFolders" ] && [ "$NumberLogFiles" == "$count" ]; then
    echo $(red "$FOLDER: Passed.")
fi

cd ..
fi

done
done

fi

#-------------------------------------------- GenerateAggregateFile module

if [ ! -z $GenerateAggregateFile ] 
then 
   printf "Argument GenerateAggregateFile is %s\n" "$GenerateAggregateFile"
   
#-------------------------------------------- generate the aggregation files

cd '/well/singlecell/'${PID[@]}
 
module add R-bundle-Bioconductor/3.12-foss-2020b-R-4.0.3
Rscript /gpfs3/well/singlecell/scripts/create_input_file_run.qsub.scRNAseq_enrichment_analysis.R --project_id ${PID[@]} --backend.data.dir /well/singlecell/references/reference_gene_sets/

find . -name "*_aggregation*" | while read line; do 
    echo  $(red "$line file has been generated.");  done

fi

#-------------------------------------------- GenerateCSV module
if [ ! -z $GenerateCSV ] 
then 
   printf "Argument GenerateCSV is %s\n" "$GenerateCSV"

#---- generate csv files per group

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME

DIR='/well/singlecell/'$NAMEFOLDER'/'

cd $DIR

if [ -d "scripts_mat2csv" ]; then
echo $(red "scripts_mat2csv/ folder exists; renamed.")
mv scripts_mat2csv scripts_mat2csv.$current_time
fi

mkdir $DIR'/scripts_mat2csv/'

declare x=0
for ADDRESS in $(find . -name 'filtered_feature_bc_matrix.h5' | sed 's/filtered_feature_bc_matrix.h5//'); do 

x=$((x + 1));
 
NAME=$(echo ${ADDRESS#*/} | cut -f1 -d"/")
 
## Input files
CELL_DIR=$(echo ${ADDRESS#*/});
 
echo -e '#!/bin/bash
#$ -cwd
#$ -N mat2csv -j y
#$ -P htseq.prjc.low
#$ -q htseq.qd@@htseq.hgd
#$ -pe shmem 2

module add CellRanger/7.0.0

cd '$DIR'
cellranger mat2csv '$ADDRESS'filtered_feature_bc_matrix.h5 '$ADDRESS'read_count.csv' > $DIR'/scripts_mat2csv/mat2csv_'$x'.sh'

cd $DIR'/scripts_mat2csv/'
qsub -q short.qc -P buck.prjc $DIR'/scripts_mat2csv/mat2csv_'$x'.sh';

done
cd '/well/singlecell/'

done

 
#-------------------------------------------- wait for scripts
y=$(qstat | grep mat2csv | wc -l)
echo $(red "$y running mat2csv scripts.")

while [ $y -ne 0 ]
do  

y=$(qstat | grep mat2csv | wc -l)

sleep 2  
done

echo $(red "mat2csv step finished.")

fi

#-------------------------------------------- Enrichment module
if [ ! -z $Enrichment ] 
then 
   printf "Argument Enrichment is %s\n" "$Enrichment"
#-------------------------------------------- generating input files for report 

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER

TYPE='grouped'
FILE='gex_'$TYPE'_aggregation'

if [ -f "$FILE" ]; then
 
    t=$(wc -l 'gex_'$TYPE'_aggregation')
    qsub -t 1-${t%% *} /well/singlecell/scripts/run.qsub.scRNAseq_enrichment_analysis_qsub.sh $PWD'/gex_'$TYPE'_aggregation'

fi

TYPE='ungrouped'
FILE='gex_'$TYPE'_aggregation'

if [ -f "$FILE" ]; then
 
    t=$(wc -l 'gex_'$TYPE'_aggregation')
    qsub -t 1-${t%% *} /well/singlecell/scripts/run.qsub.scRNAseq_enrichment_analysis_qsub.sh $PWD'/gex_'$TYPE'_aggregation'
fi

done

#-------------------------------------------- wait for scripts
y=$(qstat | grep Enrichment | wc -l)
echo $(red "$y running Enrichment scripts.")

while [ $y -ne 0 ]
do  

y=$(qstat | grep Enrichment | wc -l)

sleep 2  
done

echo $(red "Enrichment step finished.")

#-------------------------------------------- check Enrichment log files
cd '/well/singlecell/'${PID[@]}'/'
count=0
all=0
echo 'Enrichment Analysis:' >> $log_fileName

for f in $(find -maxdepth 1 -name 'Enrichment.*')
do
     if grep --quiet $f -e 'Done'; then
     ((count++))
     else
     echo $f >> $log_fileName
     fi
     ((all++))
done
 
if [ "$all" == "$count" ]; then
    echo $(red "All enrichment analysis script(s) succeed.");
    else
    echo $(red "$((all-count)) enrichment analysis script(s) failed." )
fi

fi

#-------------------------------------------- ConvertPDF module
if [ ! -z $ConvertPDF ] 
then 
   printf "Argument ConvertPDF is %s\n" "$ConvertPDF"
#-------------------------------------------- convert pdf to png

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER'/Inputs/'

for NAME in $(find . -name '*.pdf' -printf "%f\n" | sed 's/.pdf//'); do
for ADDRESS in $(find . -name $NAME'.pdf' | sed 's/.pdf//'); do 

echo "$NAME"
echo "$ADDRESS" 

echo -e '#!/bin/bash
#$ -cwd
#$ -N convertPDGtPNG -j y
#$ -P htseq.prjc.low
#$ -q htseq.qd@@htseq.hgd
#$ -pe shmem 1

cd /well/singlecell/'$NAMEFOLDER'/Inputs/

convert -density 300 -trim '$ADDRESS'.pdf' -quality 100 $ADDRESS'.png;
rm '$ADDRESS'.pdf;' > '/well/singlecell/'$NAMEFOLDER'/scripts_mat2csv/convertPDGtPNG_'$NAME'.sh'
 
qsub -q short.qc -P buck.prjc '/well/singlecell/'$NAMEFOLDER'/scripts_mat2csv/convertPDGtPNG_'$NAME'.sh';

done
done
done

#-------------------------------------------- wait for scripts
y=$(qstat | grep convertPDGtPNG | wc -l)
echo $(red "$y running convertPDGtPNG scripts.")

while [ $y -ne 0 ]
do  

y=$(qstat | grep convertPDGtPNG | wc -l)

sleep 2  
done

echo $(red "convertPDGtPNG step finished.")

cd '/well/singlecell/'$NAMEFOLDER'/Inputs/'
rm *.sh.* convertPDGtPNG.* -f
echo $(red "Log files removed.")

fi

#-------------------------------------------- OrgenizeFilesFolders module
if [ ! -z $OrgenizeFilesFolders ] 
then 
   printf "Argument OrgenizeFilesFolders is %s\n" "$OrgenizeFilesFolders"
#-------------------------------------------- add the metrics_summary.csv to the Inputs folder
types=('vdj' 'gex' 'feat' 'arc' 'atac' 'vdj-grouped' 'gex-grouped' 'feat-grouped' 'arc-grouped' 'atac-grouped')

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER

for TYPE in "${types[@]}"
do

FOLDER='10X-'$TYPE

if [ -d "$FOLDER" ]; then

echo "$FOLDER:"

mkdir -p /well/singlecell/$NAMEFOLDER/Inputs/10X-$TYPE/

cd 10X-$TYPE/
tar czf /well/singlecell/$NAMEFOLDER/Inputs/10X-$TYPE/metrics_summary_$TYPE'_'$NAMEFOLDER.tar.gz */outs/*summary.csv

cd /well/singlecell/$NAMEFOLDER/Inputs/10X-$TYPE/
tar -xvzf metrics_summary_$TYPE'_'$NAMEFOLDER.tar.gz
rm 'metrics_summary_'$TYPE'_'$NAMEFOLDER'.tar.gz'

cd /well/singlecell/$NAMEFOLDER/

fi

done

#--- samples.arc
if [ -d "10X-arc" ]; then
cp /well/singlecell/$NAMEFOLDER/samples.arc /well/singlecell/$NAMEFOLDER/Inputs/10X-arc/
fi

#--- copy samples.metadata
cp /well/singlecell/$NAMEFOLDER/samples.metadata /well/singlecell/$NAMEFOLDER/Inputs/

done


#--- for HTO projects
for NAME in "${PID[@]}"
do
NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-aggregate" ]; then
cd '10X-aggregate'

mkdir /well/singlecell/$NAMEFOLDER/Inputs/10X-aggregate/

tar czf /well/singlecell/$NAMEFOLDER/Inputs/10X-aggregate/metadata_$NAMEFOLDER.tar.gz */objects/raw_matrix/metadata.tsv.gz

cd /well/singlecell/$NAMEFOLDER/Inputs/10X-aggregate/

tar -xvzf metadata_$NAMEFOLDER.tar.gz
rm metadata_$NAMEFOLDER.tar.gz

fi
done

fi

#-------------------------------------------- FTPLinks module
if [ ! -z $FTPLinks ] 
then 
   printf "Argument FTPLinks is %s\n" "$FTPLinks"
NAMEFOLDER=${PID[@]}
DIR='/well/singlecell/'$NAMEFOLDER'/'

cd $DIR

/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER
/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER > /well/singlecell/$NAMEFOLDER/Inputs/PInf.txt
/apps/htseq/bin/getProjectDetails $NAMEFOLDER >> /well/singlecell/$NAMEFOLDER/Inputs/PInf.txt

chmod -R u+wr,g+wr ./

echo $(red "Ready for next steps.")
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running All Modules 
else
echo $(red "Running All Modules")

#-------------------------------------------- check log files
file_name=${PID[@]}
current_time=$(date "+%Y.%m.%d-%H.%M")
log_fileName=/well/singlecell/_logs_/$file_name.$current_time.log
echo "Log file name: " "$log_fileName"

types=('vdj' 'gex' 'feat' 'arc' 'vdj-grouped' 'gex-grouped' 'feat-grouped')

for NAME in "${PID[@]}"
do

echo $(red "${PID[@]}")
echo ${PID[@]} >> $log_fileName

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER

if [ -d "Inputs" ]; then
echo $(red "Inputs/ folder exists; renamed.")
mv Inputs Inputs.$current_time
fi

for TYPE in "${types[@]}"
do

FOLDER='10X-'$TYPE

if [ -d "$FOLDER" ]; then

if echo $TYPE | grep --quiet 'grouped'; then
     TYPE_revised=${TYPE%"-grouped"}.grouped
     NumberExpectedFolders=$(cat samples.$TYPE_revised | wc -l)
     else
     NumberExpectedFolders=$(cat samples.$TYPE | wc -l)
fi

cd $FOLDER

# count "find" results
NumberLogFiles=$(find -maxdepth 1 -name '*.log' -exec printf %c {} + | wc -c)
NumberLogFolders=$(ls -l | grep ^d | wc -l)

count=0

echo $FOLDER >> $log_fileName

for f in $(find -maxdepth 1 -name '*.log')
do
     
     #grep -rnw $f -e 'Pipestance completed successfully!'
     
     if grep --quiet $f -e 'Pipestance completed successfully!'; then
     ((count++))
     else
      echo $(red "$f")
      echo $f >> $log_fileName

     fi
done
 
if [ "$NumberExpectedFolders" == "$NumberLogFolders" ]; then
    echo $(red "$FOLDER: MetaData and Number of folders matched.");
    else
    echo $(red "$FOLDER: MetaData and Number of folders did NOT match.");
fi

if [ "$NumberLogFiles" == "$NumberLogFolders" ] && [ "$NumberLogFiles" != "$count" ]; then
    echo $(red "$FOLDER: $((NumberLogFiles-count)) failed sample(s).")
fi

if [ "$NumberLogFiles" != "$NumberLogFolders" ] && [ "$NumberLogFiles" != "$count" ]; then
     echo $(red "$FOLDER: $((NumberLogFiles-count)) failed sample(s) and $((NumberLogFolders-NumberLogFiles)) missed log files." )
fi

if [ "$NumberLogFiles" != "$NumberLogFolders" ] && [ "$NumberLogFiles" == "$count" ]; then
    echo $(red "$FOLDER: $((NumberLogFolders-NumberLogFiles)) missed log files." )
fi

if [ "$NumberLogFiles" == "$NumberLogFolders" ] && [ "$NumberLogFiles" == "$count" ]; then
    echo $(red "$FOLDER: Passed.")
fi

cd ..
fi

done
done

#-------------------------------------------- generate the aggregation files

cd '/well/singlecell/'${PID[@]}
 
module add R-bundle-Bioconductor/3.12-foss-2020b-R-4.0.3
Rscript /gpfs3/well/singlecell/scripts/create_input_file_run.qsub.scRNAseq_enrichment_analysis.R --project_id ${PID[@]} --backend.data.dir /well/singlecell/references/reference_gene_sets/

find . -name "*_aggregation*" | while read line; do 
    echo  $(red "$line file has been generated.");  done

#~~~~ if gex exists and genome type is human/mouse

FILE='gex_others_aggregation'

if [ ! -f "$FILE" ]; then

#---- generate csv files per group

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME

DIR='/well/singlecell/'$NAMEFOLDER'/'

cd $DIR

if [ -d "scripts_mat2csv" ]; then
echo $(red "scripts_mat2csv/ folder exists; renamed.")
mv scripts_mat2csv scripts_mat2csv.$current_time
fi

mkdir $DIR'/scripts_mat2csv/'

declare x=0
for ADDRESS in $(find . -name 'filtered_feature_bc_matrix.h5' | sed 's/filtered_feature_bc_matrix.h5//'); do 

x=$((x + 1));
 
NAME=$(echo ${ADDRESS#*/} | cut -f1 -d"/")
 
## Input files
CELL_DIR=$(echo ${ADDRESS#*/});
 
echo -e '#!/bin/bash
#$ -cwd
#$ -N mat2csv -j y
#$ -P htseq.prjc.low
#$ -q htseq.qd@@htseq.hgd
#$ -pe shmem 2

module add CellRanger/7.0.0

cd '$DIR'
cellranger mat2csv '$ADDRESS'filtered_feature_bc_matrix.h5 '$ADDRESS'read_count.csv' > $DIR'/scripts_mat2csv/mat2csv_'$x'.sh'

cd $DIR'/scripts_mat2csv/'
qsub -q short.qc -P buck.prjc $DIR'/scripts_mat2csv/mat2csv_'$x'.sh';

done
cd '/well/singlecell/'

done

 
#-------------------------------------------- wait for scripts
y=$(qstat | grep mat2csv | wc -l)
echo $(red "$y running mat2csv scripts.")

while [ $y -ne 0 ]
do  

y=$(qstat | grep mat2csv | wc -l)

sleep 2  
done

echo $(red "mat2csv step finished.")


#-------------------------------------------- generating input files for report 

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER

TYPE='grouped'
FILE='gex_'$TYPE'_aggregation'

if [ -f "$FILE" ]; then
 
    t=$(wc -l 'gex_'$TYPE'_aggregation')
    qsub -t 1-${t%% *} /well/singlecell/scripts/run.qsub.scRNAseq_enrichment_analysis_qsub.sh $PWD'/gex_'$TYPE'_aggregation'

fi

TYPE='ungrouped'
FILE='gex_'$TYPE'_aggregation'

if [ -f "$FILE" ]; then
 
    t=$(wc -l 'gex_'$TYPE'_aggregation')
    qsub -t 1-${t%% *} /well/singlecell/scripts/run.qsub.scRNAseq_enrichment_analysis_qsub.sh $PWD'/gex_'$TYPE'_aggregation'
fi

done

#-------------------------------------------- wait for scripts
y=$(qstat | grep Enrichment | wc -l)
echo $(red "$y running Enrichment scripts.")

while [ $y -ne 0 ]
do  

y=$(qstat | grep Enrichment | wc -l)

sleep 2  
done

echo $(red "Enrichment step finished.")

#-------------------------------------------- check Enrichment log files
cd '/well/singlecell/'${PID[@]}'/'
count=0
all=0
echo 'Enrichment Analysis:' >> $log_fileName

for f in $(find -maxdepth 1 -name 'Enrichment.*')
do
     if grep --quiet $f -e 'Done'; then
     ((count++))
     else
     echo $f >> $log_fileName
     fi
     ((all++))
done
 
if [ "$all" == "$count" ]; then
    echo $(red "All enrichment analysis script(s) succeed.");
    else
    echo $(red "$((all-count)) enrichment analysis script(s) failed." )
fi

#-------------------------------------------- convert pdf to png

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER'/Inputs/'

for NAME in $(find . -name '*.pdf' -printf "%f\n" | sed 's/.pdf//'); do
for ADDRESS in $(find . -name $NAME'.pdf' | sed 's/.pdf//'); do 

echo "$NAME"
echo "$ADDRESS" 

echo -e '#!/bin/bash
#$ -cwd
#$ -N convertPDGtPNG -j y
#$ -P htseq.prjc.low
#$ -q htseq.qd@@htseq.hgd
#$ -pe shmem 1

cd /well/singlecell/'$NAMEFOLDER'/Inputs/

convert -density 300 -trim '$ADDRESS'.pdf' -quality 100 $ADDRESS'.png;
rm '$ADDRESS'.pdf;' > '/well/singlecell/'$NAMEFOLDER'/scripts_mat2csv/convertPDGtPNG_'$NAME'.sh'
 
qsub -q short.qc -P buck.prjc '/well/singlecell/'$NAMEFOLDER'/scripts_mat2csv/convertPDGtPNG_'$NAME'.sh';

done
done
done

#-------------------------------------------- wait for scripts
y=$(qstat | grep convertPDGtPNG | wc -l)
echo $(red "$y running convertPDGtPNG scripts.")

while [ $y -ne 0 ]
do  

y=$(qstat | grep convertPDGtPNG | wc -l)

sleep 2  
done

echo $(red "convertPDGtPNG step finished.")

cd '/well/singlecell/'$NAMEFOLDER'/Inputs/'
rm *.sh.* convertPDGtPNG.* -f
echo $(red "Log files removed.")

fi
#~~~~ fi gex exists and genome type is human/mouse

#-------------------------------------------- add the metrics_summary.csv to the Inputs folder
types=('vdj' 'gex' 'feat' 'arc' 'atac' 'vdj-grouped' 'gex-grouped' 'feat-grouped' 'arc-grouped' 'atac-grouped')

for NAME in "${PID[@]}"
do

NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER

for TYPE in "${types[@]}"
do

FOLDER='10X-'$TYPE

if [ -d "$FOLDER" ]; then

echo "$FOLDER:"

mkdir -p /well/singlecell/$NAMEFOLDER/Inputs/10X-$TYPE/

cd 10X-$TYPE/
tar czf /well/singlecell/$NAMEFOLDER/Inputs/10X-$TYPE/metrics_summary_$TYPE'_'$NAMEFOLDER.tar.gz */outs/*summary.csv

cd /well/singlecell/$NAMEFOLDER/Inputs/10X-$TYPE/
tar -xvzf metrics_summary_$TYPE'_'$NAMEFOLDER.tar.gz
rm 'metrics_summary_'$TYPE'_'$NAMEFOLDER'.tar.gz'

cd /well/singlecell/$NAMEFOLDER/

fi

done

#--- samples.arc
if [ -d "10X-arc" ]; then
cp /well/singlecell/$NAMEFOLDER/samples.arc /well/singlecell/$NAMEFOLDER/Inputs/10X-arc/
fi

#--- copy samples.metadata
cp /well/singlecell/$NAMEFOLDER/samples.metadata /well/singlecell/$NAMEFOLDER/Inputs/

done


#--- for HTO projects
for NAME in "${PID[@]}"
do
NAMEFOLDER=$NAME
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-aggregate" ]; then
cd '10X-aggregate'

mkdir /well/singlecell/$NAMEFOLDER/Inputs/10X-aggregate/

tar czf /well/singlecell/$NAMEFOLDER/Inputs/10X-aggregate/metadata_$NAMEFOLDER.tar.gz */objects/raw_matrix/metadata.tsv.gz

cd /well/singlecell/$NAMEFOLDER/Inputs/10X-aggregate/

tar -xvzf metadata_$NAMEFOLDER.tar.gz
rm metadata_$NAMEFOLDER.tar.gz

fi
done

#-------- double check
cd '/well/singlecell/'$NAMEFOLDER'/Inputs/'
rm *.sh.* *convert* -f
echo $(red "Log files removed.")

#-------------------------------------------- generate ftp link
NAMEFOLDER=${PID[@]}
DIR='/well/singlecell/'$NAMEFOLDER'/'

cd $DIR

/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER
/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER > /well/singlecell/$NAMEFOLDER/Inputs/PInf.txt
/apps/htseq/bin/getProjectDetails $NAMEFOLDER >> /well/singlecell/$NAMEFOLDER/Inputs/PInf.txt

username=$(/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER)
ftp=${username:30:8}
mkdir /data/vftp/$ftp/additional_analyses/

cp /well/singlecell/$NAMEFOLDER/samples.metadata /data/vftp/$ftp/additional_analyses/

chmod -R u+wr,g+wr ./

echo $(red "Ready for next steps.")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fi
