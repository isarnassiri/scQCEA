### Run Cell Type Enrichment Analysis for Multiple Samples
**On an HPC Cluster:**
Create a tab-separated file (e.g., gex_aggregation) with input data set name (e.g., P220386), sample name (e.g., HAN9935A100), a gene-cell count matrix path (e.g., ~/P220386/10X-gex/HAN9935A100), a repository of reference gene sets path (e.g., ~/references/reference_gene_sets/human), and reference genome name (e.g., hsapiens) as follows:

```{r,eval=FALSE}
P220386 HAN9935A100     ~/P220386/10X-gex/HAN9935A100    ~/references/reference_gene_sets/human   hsapiens
P220386 HAN9935A101     ~/P220386/10X-gex/HAN9935A101    ~/references/reference_gene_sets/human   hsapiens
P220386 HAN9935A102     ~/P220386/10X-gex/HAN9935A102    ~/references/reference_gene_sets/human   hsapiens
```

Use this file to submit enrichment analysis scripts in parallel on the HPC cluster using `CellTypeEnrichment_Multiple_Samples.sh` script as follows:

```{r,eval=FALSE}
NAMEFOLDER=P220386
cd '~/'$NAMEFOLDER
FILE='gex_aggregation'
t=$(wc -l 'gex_aggregation')
qsub -t 1-${t%% *} CellTypeEnrichment_Multiple_Samples.sh $PWD'/gex_aggregation'
```

You need to copy the files in this directory to the folder of inputs files and modify line 5-6, and 10 in `CellTypeEnrichment_Multiple_Samples.sh` depending on the set-up of the HPC.
