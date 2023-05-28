
#### Transfer QC metrics summary to the Inputs folder

Before running the `GenerateInteractiveQCReport()` function, you need to transfer the `metrics_summary.csv` and `samples.metadata` files to the Inputs/ folder. You can use the `MetricsSummary.sh` script to automatise transferring files for a scRNA-Seq project (e.g., P230078) as follows:

You need to copy MetricsSummary.sh script to the folder of inputs files. 

```{r,eval=FALSE}
#- copy MetricsSummary.sh script to the project folder:
cp ~/MetricsSummary.sh ~/P230078/

#- change your current directory to the project folder:
cd ~/P230078/

#- make MetricsSummary.sh script excutable
chmod u+x MetricsSummary.sh

#- run the  MetricsSummary.sh script
./MetricsSummary.sh -i P230078 
