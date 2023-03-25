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

declare -A PID
PID=$1

NAMEFOLDER=${PID[@]}
DIR='/well/singlecell/'$NAMEFOLDER'/'

cd $DIR
chmod -R u+wr,g+wr ./

/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER
/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER > /well/singlecell/$NAMEFOLDER/Inputs/PInf.txt
/apps/htseq/bin/getProjectDetails $NAMEFOLDER >> /well/singlecell/$NAMEFOLDER/Inputs/PInf.txt

username=$(/apps/htseq/bin/solexaApiComplete get_project_details $NAMEFOLDER)
ftp=${username:30:8}
#--------------------------------------------

types=('vdj' 'gex' 'feat' 'arc' 'atac' 'vdj-grouped' 'gex-grouped' 'feat-grouped' 'arc-grouped' 'atac-grouped')

mkdir /data/vftp/$ftp/additional_analyses/

for TYPE in "${types[@]}"
do
cd '/well/singlecell/'$NAMEFOLDER

if [ -d "10X-"$TYPE ]; then

echo $(red "Web_summary for $TYPE exists." )
cd '10X-'$TYPE/
tar czf 10x_QC_report_$TYPE'_'$NAMEFOLDER.tar.gz */outs/web_summary.html
cp 10x_QC_report_$TYPE'_'$NAMEFOLDER.tar.gz /data/vftp/$ftp/additional_analyses/

fi
done

#-------------
cd '/well/singlecell/'$NAMEFOLDER
if [ -f "OGC_Interactive_QC_Report_"$NAMEFOLDER".zip" ]; then
echo $(red "OGC_Interactive_QC_Report exists." )
cp OGC_Interactive_QC_Report_$NAMEFOLDER.zip /data/vftp/$ftp/additional_analyses/
else
echo $(red "OGC_Interactive_QC_Report does NOT exist." )
fi

#-------------
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-gex" ]; then
echo $(red "Read_count for 10X-gex exists." )
cd '10X-gex'
tar czf ReadCount_10X-gex'_'$NAMEFOLDER.tar.gz */outs/read_count.csv
cp ReadCount_10X-gex'_'$NAMEFOLDER.tar.gz /data/vftp/$ftp/additional_analyses/
fi

#-------------
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-gex-grouped" ]; then
echo $(red "Read_count for 10X-gex-grouped exists." )
cd '10X-gex-grouped'
tar czf ReadCount_10X-gex-grouped'_'$NAMEFOLDER.tar.gz */outs/read_count.csv
cp ReadCount_10X-gex-grouped'_'$NAMEFOLDER.tar.gz /data/vftp/$ftp/additional_analyses/
 
fi

#-------------
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-feat-grouped" ]; then
echo $(red "Read_count for 10X-feat-grouped exists." )
cd '10X-feat-grouped'
tar czf ReadCount_10X-feat-grouped'_'$NAMEFOLDER.tar.gz */outs/read_count.csv
cp ReadCount_10X-feat-grouped'_'$NAMEFOLDER.tar.gz /data/vftp/$ftp/additional_analyses/
 
fi

#-------------
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-feat" ]; then
echo $(red "Read_count for 10X-feat exists." )
cd '10X-feat'
tar czf ReadCount_10X-feat'_'$NAMEFOLDER.tar.gz */outs/read_count.csv
cp ReadCount_10X-feat'_'$NAMEFOLDER.tar.gz /data/vftp/$ftp/additional_analyses/
 
fi

#-------------
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-aggregate" ]; then
echo $(red "10X-aggregate exists." )
cd '10X-aggregate'
tar czf 10X-aggregate'_'$NAMEFOLDER.tar.gz *
cp 10X-aggregate'_'$NAMEFOLDER.tar.gz /data/vftp/$ftp/additional_analyses/
 
fi

#-------------
cd '/well/singlecell/'$NAMEFOLDER
if [ -d "10X-multi" ]; then
echo $(red "10X-multi exists." )
cd '10X-multi'
tar czf '10x_QC_report_multi_'$NAMEFOLDER'.tar.gz' */outs/per_sample_outs/*/web_summary.html
cp '10x_QC_report_multi_'$NAMEFOLDER'.tar.gz' /data/vftp/$ftp/additional_analyses/
fi


#========================== check files

#--- 10X reports
types=('vdj' 'gex' 'feat' 'arc' 'atac' 'vdj-grouped' 'gex-grouped' 'feat-grouped' 'arc-grouped' 'atac-grouped')

for TYPE in "${types[@]}"
do
cd '/well/singlecell/'$NAMEFOLDER

if [ -d "10X-"$TYPE ]; then

cd '10X-'$TYPE/

#---
file1=10x_QC_report_$TYPE'_'$NAMEFOLDER.tar.gz
file2=/data/vftp/$ftp/additional_analyses/10x_QC_report_$TYPE'_'$NAMEFOLDER.tar.gz

if((`stat -c%s "$file1"`==`stat -c%s "$file2"`));then
  echo $(red "-- $TYPE: Same Size -- " )
else
  echo $(red "-- $TYPE: Different Size -- " )
#---
fi
fi
done

#--- OGC interactive reports
cd '/well/singlecell/'$NAMEFOLDER

if [ -f "OGC_Interactive_QC_Report_"$NAMEFOLDER".zip" ]; then

#---
file1='OGC_Interactive_QC_Report_'$NAMEFOLDER'.zip'
file2=/data/vftp/$ftp/additional_analyses/'OGC_Interactive_QC_Report_'$NAMEFOLDER'.zip'

if((`stat -c%s "$file1"`==`stat -c%s "$file2"`));then
  echo $(red "-- OGC_Interactive_QC_Report: Same Size -- " )
else
  echo $(red "-- OGC_Interactive_QC_Report: Different Size -- " )
#---
fi
fi

#==========================
cd /data/vftp/$ftp/additional_analyses/
find . -maxdepth 1 -type f -name '*' -type f -exec md5sum "{}" + >> checksum_other.md5
chmod -R u+wr,g+wr ./
echo $(red "Content of FTP link: " )
echo "$(ls)"

CREDENTIALS=$(/apps/htseq/bin/solexaApiComplete get_or_create_project ${PID[@]})
vars=( $CREDENTIALS )
echo 'ftp://'${vars[0]}':'${vars[1]}'@bsg-ftp.well.ox.ac.uk/additional_analyses/'
echo '/data/vftp/'${vars[0]}'/additional_analyses/'

