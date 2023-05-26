
#############################################################################################
###################################### Cell Type Enrichment #################################
#############################################################################################

#'@import DropletUtils
#'@import reshape2
#'@import dplyr
#'@import Matrix
#'@import DT
#'@import NMF
#'@import plotly
#'@import stringr
#'@import AUCell
#'@import GSEABase
#'@import rstudioapi
#'@import GEOquery
#'@importFrom data.table fread
#'@importFrom ggplot2 ggplot ggsave
#'@export
#'@name CellTypeEnrichment
#'@title Cell Type Enrichment Analysis for Gene Expression scRNA-seq data
#'@description Cell type annotation on scRNA-Seq data is a pre-step for generating an interactive QC report with scQCEA. scQCEA provides a function called CellTypeEnrichment, for automatic cell type identification and visualization on the gene-by-cell count matrix.
#'@author {Isar Nassiri}
#'@param SamplesMetadata
#'Metadata of samples including the following headers: Project Number,	LIMS ID,	Sample Name,	Index	Library Type,	Genome,	Flowcell ID,	Lane Number,	Sequencing ID
#'@param ReadCount
#'Gene-cell count matrix from 10X CellRanger count
#'@param GTF
#'We convert Ensembl ids to gene names/symbols by parsing this GTF (General Transfer Format) file
#'@param BackendDataDir
#' We used Human Protein Atlas database (version 21.0) to generate a repository of reference gene sets that are exclusively expressed in each cell type
#'@param tSNECellranger 
#'tSNE projections from 10X CellRanger count
#'@param UMAPCellranger 
#'UMAP projections from 10X CellRanger count
#'@param RawFeatureDir
#'A folder including raw feature-barcode matrices from 10X CellRanger count (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
#'@param FilteredFeatureBarcodes
#'A folder including filtered feature-barcode matrices from 10X CellRanger count (barcodes.tsv.gz)
#'@param SampleName
#'Name of an indicated sample
#'@param nCores
#'Number of cores to use for computation.
#'@param aucMaxRank
#'Number of expressed/detected genes that are going to be used to calculate the AUC. As a default, the first 250 genes are used to calculate the AUC.
#'@return You can find the results in R object under title of 'RESULTsGenomicFeatures' and 'RESULTsChromatinState'.
#'@examples
# library("scQCEA")
# csQCEAdir <- system.file("extdata", package = "scQCEA")
# DataTyep <- '10X-gex'
# SampleName <- '481207_03'
# SamplesMetadata = paste(csQCEAdir, 'Inputs/samples.metadata', sep = '/' )
# ReadCount = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'outs', sep = '/')
# GTF = paste(csQCEAdir, 'ensembl_human.txt', sep = '/')
# BackendDataDir = paste(csQCEAdir, 'ReferenceGeneSets/', sep = '/')
# tSNECellranger = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, '/outs/analysis/tsne/gene_expression_2_components', sep = '/')
# UMAPCellranger =  paste(csQCEAdir, 'Inputs', DataTyep, SampleName, '/outs/analysis/umap/gene_expression_2_components', sep = '/')
# RawFeatureDir = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'outs/raw_feature_bc_matrix', sep = '/')
# FilteredFeatureBarcodes = paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'outs/filtered_feature_bc_matrix', sep = '/')
# aucMaxRank = 250;
# nCores = 1;
# CellTypeEnrichment(SampleName, SamplesMetadata, ReadCount, GTF, BackendDataDir, tSNECellranger, UMAPCellranger, RawFeatureDir, FilteredFeatureBarcodes, aucMaxRank = 250, nCores = 1 )
#'@export
CellTypeEnrichment <- NULL
CellTypeEnrichment <- function(SampleName, SamplesMetadata, ReadCount, GTF, BackendDataDir, tSNECellranger, UMAPCellranger, RawFeatureDir, FilteredFeatureBarcodes, aucMaxRank = 250, nCores = 1 )
{
  # ---------------------------------- Create output directory
  output.dir_perSample <- paste(csQCEAdir, 'Inputs', DataTyep, SampleName, sep = '/')

  # ---------------------------------- read the scRNAseq profile
  TP_profile = fread(paste0(ReadCount, '/read_count.csv'), stringsAsFactors = FALSE)
  TP_profile = as.data.frame(TP_profile)
  row.names(TP_profile) = TP_profile$V1
  TP_profile_sub = TP_profile[,-1]
  colnames(TP_profile_sub) = gsub('-.*','',colnames(TP_profile_sub))
  
  # ---------------------------------- read repository of gene symbols and Ensemble transcript IDs
  GTF = fread(GTF, stringsAsFactors = FALSE, header = TRUE)
  GTF = as.data.frame(GTF)
  
  # ---------------------------------- add gene symbols to the gene expression profile and remove the gene IDs
  TP_profile_sub = data.frame(gene_id = row.names(TP_profile_sub), TP_profile_sub)
  TP_profile_sub = merge(GTF, TP_profile_sub, by = 'gene_id')
  
  TP_profile_sub$gene_name = make.names(TP_profile_sub$gene_name,unique=TRUE)
  row.names(TP_profile_sub) = TP_profile_sub$gene_name
  TP_profile_sub = TP_profile_sub[,-c(1,2)]
  
  # ---------------------------------- enrichment input
  list = list.files(BackendDataDir, pattern = '.tsv')
  
  my_read_csv <- function(x) {
    out <- fread(x, stringsAsFactors = FALSE, select = 'Gene')
    site <- gsub('blood_cell_category_rna_|blood_cell_category_rna_|.tsv.gz|_lineage|_Lineage|_Group','', x)
    site <- gsub('blood_cell_category_rna_|_Cell','', site)
    cbind(cellType=site, out)
  }
  
  repository <- lapply(paste(BackendDataDir,list,sep = '/'), my_read_csv)
  repository = do.call(rbind.data.frame, repository)
  
  # ---------------------------------- create a list of gene sets
  celltype = unique(repository$cellType)
  geneSets <- list()
  for(i in 1:length(celltype))
  {
    temp = repository[which(repository$cellType == celltype[i]),]
    geneSets[[i]] <- temp$Gene
  }
  
  names(geneSets) = gsub('.tsv', '', gsub(paste0(BackendDataDir, '/'), '', celltype) )
  
  # ---------------------------------- AUCell
  exprMatrix <- as.matrix(TP_profile_sub)
  cells_rankings <- AUCell_buildRankings(exprMatrix)

  # ---------------------------------- reports Genes in the gene sets NOT available in the dataset
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores=nCores, aucMaxRank = aucMaxRank)
  
  # ---------------------------------- cells assignment
  selectedThresholds <- NULL
  par(mfrow=c(3,3))  
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)
  
  ## ------ export cell Assignment as text
  cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
  assignmentTable <- melt(cellsAssigned, value.name="cell")
  colnames(assignmentTable)[2] <- "geneSet"
  
  # ---------------------------------- save cell assignments as a heatmap
  assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
  
  assignmentMat = assignmentMat[!duplicated(rowSums(assignmentMat)),]
  dim(assignmentMat)
  
  aheatmap(assignmentMat, scale="none", color="-RdYlBu2:100", legend=FALSE, filename = paste0(output.dir_perSample, '/', 'Celltype_assignment_HeatMap.png'), height = 10, width = 8)
  
  # ---------------------------------- check if any cells assigned in more than one cell type
  Freq_geneSet <- data.frame(table(assignmentTable[,"geneSet"]), stringsAsFactors = FALSE)
  Freq_geneSet$Var1 = as.character(Freq_geneSet$Var1)
  colnames(Freq_geneSet)[1] = 'geneSet'
  
  # ---------------------------------- [Keep one assignment per cell - keep most frequently enriched one] Order data frame rows according to vector with specific order
  assignmentTable_rearranged = left_join(data.frame(geneSet=Freq_geneSet$geneSet),assignmentTable,by="geneSet")
  assignmentTable_dedup = assignmentTable_rearranged[!duplicated(assignmentTable_rearranged$cell),]
  
  # ---------------------------------- add number per enrichment cluster
  Freq_assignmentTable_dedup = data.frame(table(assignmentTable_dedup[,"geneSet"]), stringsAsFactors = FALSE)
  Freq_assignmentTable_dedup = Freq_assignmentTable_dedup[order(Freq_assignmentTable_dedup$Freq, decreasing = TRUE),]
  Freq_assignmentTable_dedup$cluster = 1:dim(Freq_assignmentTable_dedup)[1]
  colnames(Freq_assignmentTable_dedup)[1] = 'geneSet'
  
  assignmentTable_dedup_ClusterNumber = left_join(assignmentTable_dedup, Freq_assignmentTable_dedup, by = 'geneSet')
  assignmentTable_dedup_ClusterNumber = assignmentTable_dedup_ClusterNumber[order(assignmentTable_dedup_ClusterNumber$Freq,decreasing = TRUE),]
  
  # ---------------------------------- save cell assignments as a text file
  write.table(assignmentTable_dedup, paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'Celltype_assignment_toCells.txt', sep = '/'), quote = FALSE, row.names = FALSE, sep = '\t')
  write.table(assignmentTable_dedup_ClusterNumber, paste(csQCEAdir, 'Inputs', DataTyep, SampleName, 'Celltype_assignment_toCells_PlusClusters.txt', sep = '/'), quote = FALSE, row.names = FALSE, sep = '\t')
  
  # ---------------------------------- superimpose the cell type on cell Range clusters [tSNE]
  tSNE_Cellranger = fread(paste0(tSNECellranger,'/projection.csv'), stringsAsFactors = FALSE)
  tSNE_Cellranger = as.data.frame(tSNE_Cellranger)
  row.names(tSNE_Cellranger) = gsub('-.*', '',  tSNE_Cellranger$Barcode)
  
  # ---------------------------------- keep the enriched cells [?]
  tSNE_Cellranger = tSNE_Cellranger[which(row.names(tSNE_Cellranger) %in% assignmentTable_dedup_ClusterNumber$cell),]
  dim(tSNE_Cellranger)
  
  # ----------------------------------
  cellsTsne = as.matrix(tSNE_Cellranger[,c(2,3)])
  
  # ---------------------------------- select 6 top dominant enriched cell types
  selectedThresholds <- getThresholdSelected(cells_assignment)
  selectedCellTypes = unique(assignmentTable_dedup_ClusterNumber$geneSet)[1:6]
  selectedThresholds = selectedThresholds[which(names(selectedThresholds)%in%selectedCellTypes)]
  
  # ---------------------------------- visualization of tSNE plot + functional annotation
  png(paste0(output.dir_perSample, '/', 'tSNE_Plot.png'))
  par(mar=c(2,2,2,2))
  for(geneSetName in selectedCellTypes[1])
  {
    nBreaks <- 5 # Number of levels in the color palettes
    # Color palette for the cells that do not pass the threshold
    colorPal_Neg <- colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
    # Color palette for the cells that pass the threshold
    colorPal_Pos <- colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
    
    passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
    
    if(sum(passThreshold) > 0 )
    {
      aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
      
      # Assign cell color
      cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])),
                     setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
      
      plot(cellsTsne, main=geneSetName,
           sub="Pink/red cells pass the threshold",
           col=cellColor[rownames(cellsTsne)], pch=16, cex = 0.8)
    }
  }
  dev.off()  
  
  # ---------------------------------- superimpose the cell type on cell Range clusters [UMAP]
  UMAP_Cellranger = fread(paste0(UMAPCellranger, '/projection.csv'), stringsAsFactors = FALSE)
  UMAP_Cellranger = as.data.frame(UMAP_Cellranger)
  row.names(UMAP_Cellranger) = gsub('-.*', '',  UMAP_Cellranger$Barcode)
  
  # ---------------------------------- keep the enriched cells [?]
  UMAP_Cellranger = UMAP_Cellranger[which(row.names(UMAP_Cellranger) %in% assignmentTable_dedup_ClusterNumber$cell),]
  dim(UMAP_Cellranger)
  
  # ----------------------------------
  cellsUMAP = as.matrix(UMAP_Cellranger[,c(2,3)])
  
  # ---------------------------------- select 6 top dominant enriched cell types
  selectedThresholds <- getThresholdSelected(cells_assignment)
  selectedCellTypes = unique(assignmentTable_dedup_ClusterNumber$geneSet)[1:6]
  selectedThresholds = selectedThresholds[which(names(selectedThresholds)%in%selectedCellTypes)]
  
  # ---------------------------------- visualization of UMAP plot + functional annotation
  png(paste0(output.dir_perSample, '/', 'UMAP_Plot.png'))
  par(mar=c(2,2,2,2))
  
  for(geneSetName in selectedCellTypes[1])
  {
    nBreaks <- 5 # Number of levels in the color palettes
    # Color palette for the cells that do not pass the threshold
    colorPal_Neg <- colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
    # Color palette for the cells that pass the threshold
    colorPal_Pos <- colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
    
    passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
    
    if(sum(passThreshold) > 0 )
    {
      aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
      
      # Assign cell color
      cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])),
                     setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
      
      plot(cellsUMAP, main=geneSetName,
           sub="Pink/red cells pass the threshold",
           col=cellColor[rownames(cellsUMAP)], pch=16, cex = 0.8)
    }
  }
  dev.off()
  
  # ---------------------------------- knee plot [calculation]
  RawFeatureDir = RawFeatureDir
  barcode.path <- paste(RawFeatureDir, "barcodes.tsv.gz", sep = '/')
  features.path <- paste(RawFeatureDir, "features.tsv.gz", sep = '/')
  matrix.path <- paste(RawFeatureDir, "matrix.mtx.gz", sep = '/')
  
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  br.out <- barcodeRanks(mat)
  
  # ---------------------------------- knee plot - Refined - visualization
  Input = data.frame(cell = row.names(br.out), ranking = br.out$rank, originalFreq = br.out$total, whiteList = rep(NA, length( br.out$total)))
  
  EA_WhiteList = assignmentTable_dedup_ClusterNumber
  EA_WhiteList$cell = paste0(EA_WhiteList$cell, '-1')
  
  Input$whiteList[which(Input$cell %in% EA_WhiteList$cell)] = TRUE
  Input$whiteList[-which(Input$cell %in% EA_WhiteList$cell)] = FALSE
  EstimatedNumberCells = as.numeric(table(Input$whiteList)[2])
  
  g1 <- ggplot(Input, aes(x = .data$ranking, y = .data$originalFreq)) +
    geom_line(alpha = 0.5, size = 2, aes(color = whiteList)) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Cell Barcode Rank") +
    ylab("Cell Barcode Frequency (UMI counts)") +
    labs(title="Barcode Rank Plot - Refined") +
    theme_bw() +
    theme(text = element_text(size=15), legend.text = element_text(size = 15), axis.title = element_text(size = 15), plot.title = element_text(size = 15, face = "bold", color = "black")) +
    scale_color_manual('', labels = c("Background", "Cells"), values = c(`TRUE` = "darkblue", `FALSE` = "gray")) +
    geom_label(data = Input %>% filter(whiteList) %>% filter(ranking == max(ranking)), aes(label = paste0("( Number of Cells: ", EstimatedNumberCells, ")")), hjust = 0, nudge_x = -0.5, nudge_y = -0.2, size = 4) +
    geom_hline(yintercept = metadata(br.out)$knee, linetype='dotted', col = 'darkgreen') +
    annotate("text", x = 20, y = metadata(br.out)$knee, label = "Knee", vjust = -0.5, size = 5) +
    geom_hline(yintercept = metadata(br.out)$inflection, linetype='dotted', col = 'brown') +
    annotate("text", x = 20, y = metadata(br.out)$inflection, label = "Inflection", vjust = -0.5, size = 5)
  
  ggsave(
    paste0(output.dir_perSample, '/', 'BarcodeRankPlot_EB.png'),
    plot = g1,
    device = "png",
    dpi = 600
  )
  
  # ---------------------------------- knee plot - Standard - visualization
  Input = data.frame(cell = row.names(br.out), ranking = br.out$rank, originalFreq = br.out$total, whiteList = rep(NA, length( br.out$total)))
  
  FilteredFeatureBarcodes = fread(paste0(FilteredFeatureBarcodes, '/barcodes.tsv.gz'), stringsAsFactors = FALSE, header = FALSE)
  length(which(Input$cell %in% FilteredFeatureBarcodes$V1))
  
  Input$whiteList[which(Input$cell %in% FilteredFeatureBarcodes$V1)] = TRUE
  Input$whiteList[-which(Input$cell %in% FilteredFeatureBarcodes$V1)] = FALSE
  
  EstimatedNumberCells = as.numeric(table(Input$whiteList)[2])
  
  g2 <- ggplot(Input, aes(x = .data$ranking, y = .data$originalFreq)) +
    geom_line(alpha = 0.5, size = 2, aes(color = whiteList)) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Cell Barcode Rank") +
    ylab("Cell Barcode Frequency (UMI counts)") +
    labs(title="Barcode Rank Plot - Standard") +
    theme_bw() +
    theme(text = element_text(size=15), legend.text = element_text(size = 15), axis.title = element_text(size = 15), plot.title = element_text(size = 15, face = "bold", color = "black")) +
    scale_color_manual('', labels = c("Background", "Cells"), values = c(`TRUE` = "darkblue", `FALSE` = "gray")) +
    geom_label(data = Input %>% filter(whiteList) %>% filter(ranking == max(ranking)), aes(label = paste0("( Number of Cells: ", EstimatedNumberCells, ")")), hjust = 0, nudge_x = -0.5, nudge_y = -0.2, size = 4) +
    geom_hline(yintercept = metadata(br.out)$knee, linetype='dotted', col = 'darkgreen') +
    annotate("text", x = 20, y = metadata(br.out)$knee, label = "Knee", vjust = -0.5, size = 5) +
    geom_hline(yintercept = metadata(br.out)$inflection, linetype='dotted', col = 'brown') +
    annotate("text", x = 20, y = metadata(br.out)$inflection, label = "Inflection", vjust = -0.5, size = 5)
  
  ggsave(
    paste0(output.dir_perSample, '/', 'BarcodeRankPlot_10X.png'),
    plot = g2,
    device = "png",
    dpi = 600
  )
  
  # ---------------------------------- Total UMI vs Number of Detected Genes - visualization
  colnames(TP_profile_sub) = paste0(colnames(TP_profile_sub), '-1')
  
  read_count = TP_profile_sub[apply(TP_profile_sub, 1, function(x) !all(x==0)),]
  dim(read_count)
  
  TotalUMIcount = data.frame(Barcode = colnames(read_count), totalUMICount = colSums(read_count), whiteList = rep(NA, length( colnames(read_count) )))
  
  EA_WhiteList = assignmentTable_dedup_ClusterNumber
  EA_WhiteList$cell = paste0(EA_WhiteList$cell, '-1')
  
  TotalUMIcount$whiteList[which(TotalUMIcount$Barcode %in% EA_WhiteList$cell)] = TRUE
  TotalUMIcount$whiteList[-which(TotalUMIcount$Barcode %in% EA_WhiteList$cell)] = FALSE
  
  read_count2 = read_count %>% mutate_if(is.numeric, ~1 * (. != 0))
  
  TotalUMIcount$nbrGenesAboveZero = colSums(read_count2)
  
  TotalUMIcount = TotalUMIcount[order(TotalUMIcount$nbrGenesAboveZero, decreasing = TRUE),]
  TotalUMIcount$rank = 1:dim(TotalUMIcount)[1]
  
  TotalUMIcount$whiteList[TotalUMIcount$whiteList == 'TRUE'] = "Cells"
  TotalUMIcount$whiteList[TotalUMIcount$whiteList == 'FALSE'] = "Background"
  table(TotalUMIcount$whiteList)

  g3 <- ggplot(TotalUMIcount, aes(x = .data$nbrGenesAboveZero, y = .data$totalUMICount)) +
    geom_point(alpha = 0.5, size = 3, aes(color = whiteList)) +
    geom_point(data = TotalUMIcount[which(TotalUMIcount$whiteList == "Background"),], color = "red") + 
    xlab("Total UMI Count") +
    ylab("Number of Detected Genes") +
    labs(title="Quantification Plot") +
    theme_bw() +
    theme(text = element_text(size=15), legend.text = element_text(size = 15), axis.title = element_text(size = 15), plot.title = element_text(size = 15, face = "bold", color = "black"))  

  ggsave(
    paste0(output.dir_perSample, '/', 'TotalUMIvsDetectedGenes.png'),
    plot = g3,
    device = "png",
    dpi = 600
  )
  
  # ----------------------------------
  cat(paste0("\033[0;", 47, "m", "You can find the results in: ", "\033[0m","\n", output.dir_perSample))

}

