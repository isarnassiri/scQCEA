

# ---------------------------------- essential libraries
suppressPackageStartupMessages({
  library(AUCell)
  library(GSEABase)
  library(data.table)
  library(DT)
  library(NMF)
  library(plotly)
  library(GEOquery)
  library(data.table) ## 1.9.2 or 1.9.3
  library(stringr)
  library(argparser)
  library(optparse)
  library(Matrix)
  library(DropletUtils)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# ---------------------------------- Inputs Paths
INPUT_PATH = "/well/singlecell/"
BACK_END_DATA = "/well/singlecell/references/ensembl_gene_name_id/"

# ---------------------------------- Parse arguments
option_list <- list(
  make_option(c("--project_id"), default = NULL,
              help = "Project ID to be used to name the main Seurat object."),
  make_option(c("--gex_library_id"), default = NULL,
              help = "The Library ID of the Gene Expression sequenced sample (re-sequenced samples will have the same Library ID)."),
  make_option(c("--input.dir"), default = NULL,
              help = "Path to the input files directory."),
  make_option(c("--backend.data.dir"), default = NULL,
              help = "Path to the backend directory."),
  make_option(c("--organism"), default = NULL,
              help = "Type of genome")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------------
project.id <- opt$project_id
gex.library.id <- opt$gex_library_id
input.dir <- opt$input.dir
backend.data.dir <- opt$backend.data.dir
organism <- opt$organism

if (length(grep('grouped', input.dir)) == 1)
  # catch integer(0)
{
  output_dir = paste0(INPUT_PATH, opt$project_id, "/Inputs/10X-gex-grouped/")
} else{
  output_dir = paste0(INPUT_PATH, opt$project_id, "/Inputs/10X-gex/")
}

# ---------------------------------- Validating arguments
if (is.null(opt$project_id)) {
  stop(paste0("Argument not provided. --project_id"))
}
if (!is.null(opt$gex_library_id) &&
    length(grep("[^a-zA-Z0-9_-]", opt$gex_library_id)) > 0) {
  stop(
    paste0(
      "Invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --gex_library_id ",
      opt$gex_library_id
    )
  )
}
if (is.null(opt$backend.data.dir)) {
  stop(
    paste0(
      "Argument not provided or directory provided does not exist. --backend.data.dir ",
      opt$backend.data.dir
    )
  )
}

# ---------------------------------- Create output directory
output.dir_perSample <- paste(output_dir, gex.library.id, sep = "/")
dir.create(output.dir_perSample,
           showWarnings = FALSE,
           recursive = TRUE) # recursive=TRUE will create parent directory first, if it doesn't exist.

# ---------------------------------- read the scRNAseq profile
TP_profile = fread(paste0(input.dir, '/outs/read_count.csv'), stringsAsFactors = F)
TP_profile = as.data.frame(TP_profile)
row.names(TP_profile) = TP_profile$V1
TP_profile_sub = TP_profile[, -1]
colnames(TP_profile_sub) = gsub('-.*', '', colnames(TP_profile_sub))

# ---------------------------------- read repository of gene symbols and Ensemble transcript IDs
if (organism == "hsapiens") {
  GTF = fread(
    paste0(BACK_END_DATA, 'ensembl_human.txt'),
    stringsAsFactors = F,
    header = T
  )
  GTF = as.data.frame(GTF)
}

if (organism == "mmusculus") {
  GTF = fread(
    paste0(BACK_END_DATA, 'ensembl_moueSymID_humanSym.txt'),
    stringsAsFactors = F,
    header = T
  )
  GTF = as.data.frame(GTF)
  GTF = GTF[, which(colnames(GTF) %in% c('gene_name', 'gene_id'))]
}

# ---------------------------------- add gene symbols to the gene expression profile and remove the gene IDs
TP_profile_sub = data.frame(gene_id = row.names(TP_profile_sub), TP_profile_sub)
TP_profile_sub = merge(GTF, TP_profile_sub, by = 'gene_id')

TP_profile_sub$gene_name = make.names(TP_profile_sub$gene_name, unique =
                                        T)
row.names(TP_profile_sub) = TP_profile_sub$gene_name
TP_profile_sub = TP_profile_sub[, -c(1, 2)]

# ---------------------------------- enrichment input
list = list.files(backend.data.dir, pattern = '.tsv')

my_read_csv <- function(x) {
  out <- fread(x, stringsAsFactors = F, select = 'Gene')
  site <-
    gsub(
      'blood_cell_category_rna_|blood_cell_category_rna_|.tsv.gz|_lineage|_Lineage|_Group',
      '',
      x
    )
  site <- gsub('blood_cell_category_rna_|_Cell', '', site)
  cbind(cellType = site, out)
}

repository <-
  lapply(paste(backend.data.dir, list, sep = '/'), my_read_csv)
repository <- do.call(rbind.data.frame, repository)

# ---------------------------------- create a list of gene sets
celltype = unique(repository$cellType)
geneSets <- list()
for (i in 1:length(celltype))
{
  temp = repository[which(repository$cellType == celltype[i]), ]
  geneSets[[i]] <- temp$Gene
}

names(geneSets) = gsub('.tsv', '', gsub(paste0(backend.data.dir, '/'), '', celltype))

# ---------------------------------- AUCell
set.seed(123)

exprMatrix <- as.matrix(TP_profile_sub)
cells_rankings <- AUCell_buildRankings(exprMatrix)

# ---------------------------------- reports Genes in the gene sets NOT available in the dataset
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores = 1)

# ---------------------------------- cells assignment

selectedThresholds <- NULL

set.seed(123)
par(mfrow = c(3, 3)) # PROVIDE ENOUGH SPACE IN PLOTS ENVIRONMENT OF THE RSTUDIO
cells_assignment <-
  AUCell_exploreThresholds(cells_AUC, plotHist = FALSE, assign = TRUE)

## ------ export cell Assignment as text
cellsAssigned <- lapply(cells_assignment, function(x)
  x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name = "cell")
colnames(assignmentTable)[2] <- "geneSet"

# ---------------------------------- save cell assignments as a heatmap
assignmentMat <-
  table(assignmentTable[, "geneSet"], assignmentTable[, "cell"])
assignmentMat = assignmentMat[!duplicated(rowSums(assignmentMat)), ]
dim(assignmentMat)

set.seed(123)
aheatmap(
  assignmentMat,
  scale = "none",
  color = "-RdYlBu2:100",
  legend = FALSE,
  filename = paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'Celltype_assignment_HeatMap.pdf'
  ),
  height = 10,
  width = 8
)

# ---------------------------------- check if any cells assigned in more than one cell type
Freq_geneSet <-
  data.frame(table(assignmentTable[, "geneSet"]), stringsAsFactors = F)
Freq_geneSet$Var1 = as.character(Freq_geneSet$Var1)
colnames(Freq_geneSet)[1] = 'geneSet'

# ---------------------------------- [Keep one assignment per cell - keep most frequently enriched one] Order data frame rows according to vector with specific order
library(dplyr)
assignmentTable_rearranged = left_join(data.frame(geneSet = Freq_geneSet$geneSet),
                                       assignmentTable,
                                       by = "geneSet")
assignmentTable_dedup = assignmentTable_rearranged[!duplicated(assignmentTable_rearranged$cell), ]

# ---------------------------------- add number per enrichment cluster
dim(assignmentTable_dedup)[1] > dim(exprMatrix)[2] # if it's TRUE means if have duplications in enrichment results

Freq_assignmentTable_dedup = data.frame(table(assignmentTable_dedup[, "geneSet"]), stringsAsFactors = F)
Freq_assignmentTable_dedup = Freq_assignmentTable_dedup[order(Freq_assignmentTable_dedup$Freq, decreasing = T), ]
Freq_assignmentTable_dedup$cluster = 1:dim(Freq_assignmentTable_dedup)[1]
colnames(Freq_assignmentTable_dedup)[1] = 'geneSet'

assignmentTable_dedup_ClusterNumber = left_join(assignmentTable_dedup, Freq_assignmentTable_dedup, by = 'geneSet')
assignmentTable_dedup_ClusterNumber = assignmentTable_dedup_ClusterNumber[order(assignmentTable_dedup_ClusterNumber$Freq, decreasing = T), ]

# ---------------------------------- save cell assignments as a text file
write.table(
  assignmentTable_dedup,
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'Celltype_assignment_toCells.txt'
  ),
  quote = F,
  row.names = F,
  sep = '\t'
)
write.table(
  assignmentTable_dedup_ClusterNumber,
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'Celltype_assignment_toCells_PlusClusters.txt'
  ),
  quote = F,
  row.names = F,
  sep = '\t'
)

# ---------------------------------- superimpose the cell type on cell Range clusters [tSNE]

if (length(grep('grouped', input.dir)) == 1)
  # catch integer(0)
{
  tSNE_Cellranger = fread(
    paste0(
      INPUT_PATH,
      project.id,
      '/10X-gex-grouped/',
      gex.library.id,
      '/outs/analysis/tsne/gene_expression_2_components/projection.csv'
    ),
    stringsAsFactors = F
  )
  
} else{
  tSNE_Cellranger = fread(
    paste0(
      INPUT_PATH,
      project.id,
      '/10X-gex/',
      gex.library.id,
      '/outs/analysis/tsne/gene_expression_2_components/projection.csv'
    ),
    stringsAsFactors = F
  )
}

tSNE_Cellranger = as.data.frame(tSNE_Cellranger)
row.names(tSNE_Cellranger) = gsub('-.*', '',  tSNE_Cellranger$Barcode)

# ---------------------------------- keep the enriched cells [?]
tSNE_Cellranger = tSNE_Cellranger[which(row.names(tSNE_Cellranger) %in% assignmentTable_dedup_ClusterNumber$cell), ]
dim(tSNE_Cellranger)

# ----------------------------------
cellsTsne = as.matrix(tSNE_Cellranger[, c(2, 3)])

# ---------------------------------- select 6 top dominant enriched cell types
selectedThresholds <- getThresholdSelected(cells_assignment)
selectedCellTypes = unique(assignmentTable_dedup_ClusterNumber$geneSet)[1:6]
selectedThresholds = selectedThresholds[which(names(selectedThresholds) %in%
                                                selectedCellTypes)]

# ---------------------------------- visualization of tSNE plot + functional annotation
set.seed(123)
pdf(
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'tSNE_Plot.pdf'
  ),
  height = 8,
  width = 10,
  useDingbats = F
)
par(mfrow = c(2, 3)) # Splits the plot into two rows and three columns

for (geneSetName in selectedCellTypes)
{
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <-
    grDevices::colorRampPalette(c("black", "blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <-
    grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  passThreshold <-
    getAUC(cells_AUC)[geneSetName, ] >  selectedThresholds[geneSetName]
  
  if (sum(passThreshold) > 0)
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName, ], passThreshold)
    
    # Assign cell color
    cellColor <-
      c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks = nBreaks)], names(aucSplit[[1]])),
        setNames(colorPal_Pos[cut(aucSplit[[2]], breaks = nBreaks)], names(aucSplit[[2]])))
    
    plot(
      cellsTsne,
      main = geneSetName,
      sub = "Pink/red cells pass the threshold",
      col = cellColor[rownames(cellsTsne)],
      pch = 16,
      cex = 0.8
    )
  }
}
dev.off()

# ---------------------------------- superimpose the cell type on cell Range clusters [UMAP]
if (length(grep('grouped', input.dir)) == 1)
  # catch integer(0)
{
  UMAP_Cellranger = fread(
    paste0(
      INPUT_PATH,
      project.id,
      '/10X-gex-grouped/',
      gex.library.id,
      '/outs/analysis/umap/gene_expression_2_components/projection.csv'
    ),
    stringsAsFactors = F
  )
  
} else{
  UMAP_Cellranger = fread(
    paste0(
      INPUT_PATH,
      project.id,
      '/10X-gex/',
      gex.library.id,
      '/outs/analysis/umap/gene_expression_2_components/projection.csv'
    ),
    stringsAsFactors = F
  )
}

UMAP_Cellranger = as.data.frame(UMAP_Cellranger)
row.names(UMAP_Cellranger) = gsub('-.*', '',  UMAP_Cellranger$Barcode)

# ---------------------------------- keep the enriched cells [?]
UMAP_Cellranger = UMAP_Cellranger[which(row.names(UMAP_Cellranger) %in% assignmentTable_dedup_ClusterNumber$cell), ]
dim(UMAP_Cellranger)

# ----------------------------------
cellsUMAP = as.matrix(UMAP_Cellranger[, c(2, 3)])

# ---------------------------------- select 6 top dominant enriched cell types
selectedThresholds <- getThresholdSelected(cells_assignment)
selectedCellTypes = unique(assignmentTable_dedup_ClusterNumber$geneSet)[1:6]
selectedThresholds = selectedThresholds[which(names(selectedThresholds) %in%
                                                selectedCellTypes)]

# ---------------------------------- visualization of UMAP plot + functional annotation
set.seed(123)
pdf(
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'UMAP_Plot.pdf'
  ),
  height = 8,
  width = 10,
  useDingbats = F
)
par(mfrow = c(2, 3)) # Splits the plot into two rows and three columns

for (geneSetName in selectedCellTypes)
{
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <-
    grDevices::colorRampPalette(c("black", "blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <-
    grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  passThreshold <-
    getAUC(cells_AUC)[geneSetName, ] >  selectedThresholds[geneSetName]
  
  if (sum(passThreshold) > 0)
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName, ], passThreshold)
    
    # Assign cell color
    cellColor <-
      c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks = nBreaks)], names(aucSplit[[1]])),
        setNames(colorPal_Pos[cut(aucSplit[[2]], breaks = nBreaks)], names(aucSplit[[2]])))
    
    plot(
      cellsUMAP,
      main = geneSetName,
      sub = "Pink/red cells pass the threshold",
      col = cellColor[rownames(cellsUMAP)],
      pch = 16,
      cex = 0.8
    )
  }
}
dev.off()

## -------------------------------------------------------------------------- knee plot

matrix_dir = paste0(input.dir, '/outs/raw_feature_bc_matrix/')
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# Computing barcode rank statistics:
br.out <- barcodeRanks(mat)

## -------------------------------------------------------------------------- knee plot - 2

#------- Enrichnment based [Cell barcode error correction]

Input = data.frame(
  cell = row.names(br.out),
  ranking = br.out$rank,
  originalFreq = br.out$total,
  whiteList = rep(NA, length(br.out$total))
)
dim(Input)

EA_WhiteList = assignmentTable_dedup_ClusterNumber
EA_WhiteList$cell = paste0(EA_WhiteList$cell, '-1')
dim(EA_WhiteList)
length(unique(EA_WhiteList$cell))

Input$whiteList[which(Input$cell %in% EA_WhiteList$cell)] = TRUE
Input$whiteList[-which(Input$cell %in% EA_WhiteList$cell)] = FALSE
table((Input$whiteList))

EstimatedNumberCells = as.numeric(table(Input$whiteList)[2])

g1 <-
  ggplot(Input, aes(x = .data$ranking, y = .data$originalFreq)) +
  geom_line(alpha = 0.5, size = 2, aes(color = whiteList)) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Cell Barcode Rank") +
  ylab("Cell Barcode Frequency (UMI counts)") +
  labs(title = "Barcode Rank Plot - EB") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(
      size = 15,
      face = "bold",
      color = "black"
    )
  ) +
  scale_color_manual(
    '',
    labels = c("Cells", "Background"),
    values = c(`TRUE` = "#4DB6D0", `FALSE` = "gray")
  ) +
  geom_label(
    data = Input %>% dplyr::filter(whiteList) %>% dplyr::filter(ranking == max(ranking)),
    aes(label = paste0(
      "( Number of Cells: ", EstimatedNumberCells, ")"
    )),
    hjust = 0,
    nudge_x = 0.1,
    size = 4
  ) +
  geom_hline(
    yintercept = metadata(br.out)$knee,
    linetype = 'dotted',
    col = 'darkgreen'
  ) +
  annotate(
    "text",
    x = 20,
    y = metadata(br.out)$knee,
    label = "Knee",
    vjust = -0.5,
    size = 5
  ) +
  geom_hline(
    yintercept = metadata(br.out)$inflection,
    linetype = 'dotted',
    col = 'brown'
  ) +
  annotate(
    "text",
    x = 20,
    y = metadata(br.out)$inflection,
    label = "Inflection",
    vjust = -0.5,
    size = 5
  )

pdf(
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'BarcodeRankPlot_EB.pdf'
  ),
  height = 8,
  width = 8,
  useDingbats = F
)
g1
dev.off()

#--------- Knee plot - 10X based
Input = data.frame(
  cell = row.names(br.out),
  ranking = br.out$rank,
  originalFreq = br.out$total,
  whiteList = rep(NA, length(br.out$total))
)

EA_WhiteList2 = fread(
  paste0(
    input.dir,
    '/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
  ),
  stringsAsFactors = F,
  header = F
)
dim(EA_WhiteList2)
length(which(Input$cell %in% EA_WhiteList2$V1))
dim(Input)

Input$whiteList[which(Input$cell %in% EA_WhiteList2$V1)] = TRUE
Input$whiteList[-which(Input$cell %in% EA_WhiteList2$V1)] = FALSE
table(is.na(Input$whiteList))

EstimatedNumberCells = as.numeric(table(Input$whiteList)[2])

g2 <-
  ggplot(Input, aes(x = .data$ranking, y = .data$originalFreq)) +
  geom_line(alpha = 0.5, size = 2, aes(color = whiteList)) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Cell Barcode Rank") +
  ylab("Cell Barcode Frequency (UMI counts)") +
  labs(title = "Barcode Rank Plot - SKP") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(
      size = 15,
      face = "bold",
      color = "black"
    )
  ) +
  scale_color_manual(
    '',
    labels = c("Cells", "Background"),
    values = c(`TRUE` = "#4DB6D0", `FALSE` = "gray")
  ) +
  geom_label(
    data = Input %>% dplyr::filter(whiteList) %>% dplyr::filter(ranking == max(ranking)),
    aes(label = paste0(
      "( Number of Cells: ", EstimatedNumberCells, ")"
    )),
    hjust = 0,
    nudge_x = 0.1,
    size = 4
  ) +
  geom_hline(
    yintercept = metadata(br.out)$knee,
    linetype = 'dotted',
    col = 'darkgreen'
  ) +
  annotate(
    "text",
    x = 20,
    y = metadata(br.out)$knee,
    label = "Knee",
    vjust = -0.5,
    size = 5
  ) +
  geom_hline(
    yintercept = metadata(br.out)$inflection,
    linetype = 'dotted',
    col = 'brown'
  ) +
  annotate(
    "text",
    x = 20,
    y = metadata(br.out)$inflection,
    label = "Inflection",
    vjust = -0.5,
    size = 5
  )

pdf(
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'BarcodeRankPlot_10X.pdf'
  ),
  height = 8,
  width = 8,
  useDingbats = F
)
g2
dev.off()

#--------------

#-- TotalUMIcount
colnames(TP_profile_sub) = paste0(colnames(TP_profile_sub), '-1')

read_count = TP_profile_sub[apply(TP_profile_sub, 1, function(x)
  ! all(x == 0)), ]
dim(read_count)

TotalUMIcount = data.frame(
  Barcode = colnames(read_count),
  totalUMICount = colSums(read_count),
  whiteList = rep(NA, length(colnames(read_count)))
)

setwd(output.dir_perSample)
EA_Cluster = fread(paste0(
  output.dir_perSample,
  '/',
  list.files(pattern = '_Celltype_assignment_toCells_PlusClusters.txt')
),
stringsAsFactors = F)
EA_Cluster$cell = paste0(EA_Cluster$cell, '-1')
dim(EA_Cluster)
length(unique(EA_Cluster$cell))

TotalUMIcount$Cluster[which(TotalUMIcount$Barcode %in% EA_Cluster$cell)] = TRUE
TotalUMIcount$Cluster[-which(TotalUMIcount$Barcode %in% EA_Cluster$cell)] = FALSE
table((TotalUMIcount$Cluster))

read_count2 = read_count %>% mutate_if(is.numeric, ~ 1 * (. != 0))
TotalUMIcount$nbrGenesAboveZero = colSums(read_count2)

TotalUMIcount = TotalUMIcount[order(TotalUMIcount$nbrGenesAboveZero, decreasing = T), ]
TotalUMIcount$rank = 1:dim(TotalUMIcount)[1]

TotalUMIcount$Cluster[TotalUMIcount$Cluster == 'TRUE'] = "Cells"
TotalUMIcount$Cluster[TotalUMIcount$Cluster == 'FALSE'] = "Background"
table(TotalUMIcount$Cluster)

library(ggrepel)
g3 <-
  ggplot(TotalUMIcount,
         aes(x = .data$nbrGenesAboveZero, y = .data$totalUMICount)) +
  geom_point(alpha = 0.5, size = 3, aes(color = Cluster)) +
  geom_point(data = TotalUMIcount[which(TotalUMIcount$Cluster == "Background"), ], color = "red") +
  xlab("Total UMI Count") +
  ylab("Number of Detected Genes") +
  labs(title = "Quantification Plot") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(
      size = 15,
      face = "bold",
      color = "black"
    )
  )

pdf(
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'TotalUMIvsDetectedGenes.pdf'
  ),
  height = 8,
  width = 8,
  useDingbats = F
)
g3
dev.off()

#------- Just Enrichnment based [Cell barcode error correction]
Input = data.frame(
  cell = row.names(br.out),
  ranking = br.out$rank,
  originalFreq = br.out$total,
  whiteList = rep(NA, length(br.out$total))
)
dim(Input)

# EB
EA_WhiteList = fread(
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'Celltype_assignment_toCells_PlusClusters.txt'
  ),
  stringsAsFactors = F,
  header = T
)
EA_WhiteList$cell = paste0(EA_WhiteList$cell, '-1')
dim(EA_WhiteList)
length(unique(EA_WhiteList$cell))

# 10X
WhiteList_10X = fread(
  paste0(
    input.dir,
    '/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
  ),
  stringsAsFactors = F,
  header = F
)
dim(WhiteList_10X)
length(which(Input$cell %in% WhiteList_10X$V1))

setdiff(WhiteList_10X$V1, EA_WhiteList$cell)

Input$whiteList[which(Input$cell %in% setdiff(WhiteList_10X$V1, EA_WhiteList$cell))] = TRUE
Input$whiteList[-which(Input$cell %in% setdiff(WhiteList_10X$V1, EA_WhiteList$cell))] = FALSE
table((Input$whiteList))

Input = Input[(Input$whiteList), ]
dim(Input)
EstimatedNumberCells = as.numeric(table(Input$whiteList)[1])

library(ggplot2)
g4 <-
  ggplot(Input, aes(x = .data$ranking, y = .data$originalFreq)) +
  geom_point(alpha = 0.5,
             size = 3,
             color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Cell Barcode Rank") +
  ylab("Cell Barcode Frequency (UMI counts)") +
  labs(title = "Barcode Rank Plot (SKP - EB)") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(
      size = 15,
      face = "bold",
      color = "black"
    )
  ) +
  geom_hline(
    yintercept = metadata(br.out)$knee,
    linetype = 'dotted',
    col = 'darkgreen'
  ) +
  annotate(
    "text",
    x = 20,
    y = metadata(br.out)$knee,
    label = "Knee",
    vjust = -0.5,
    size = 5
  ) +
  geom_hline(
    yintercept = metadata(br.out)$inflection,
    linetype = 'dotted',
    col = 'brown'
  ) +
  annotate(
    "text",
    x = 20,
    y = metadata(br.out)$inflection,
    label = "Inflection",
    vjust = -0.5,
    size = 5
  )

pdf(
  paste0(
    output.dir_perSample,
    '/',
    project.id,
    '_',
    gex.library.id,
    '_',
    'BarcodeRankPlot_EB_FilterOut.pdf'
  ),
  height = 8,
  width = 8,
  useDingbats = F
)
g4
dev.off()

## --------------------------------------------------------------------------
date()
sessionInfo()