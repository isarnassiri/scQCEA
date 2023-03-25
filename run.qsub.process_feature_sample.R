#!/usr/bin/env Rscript

library(knitr)
library(dplyr)
library(ggplot2)
library(optparse)
library(Matrix)
library(R.utils)
library(Seurat)

## Plot constants
plot.units <- "cm"
plot.size <- 10
plot.dpi <- 150

## Saving output files
saveHashOutput = function(where, folder, object, metadata, gex.data.dir, adt.data.dir, hto.data.dir, remove_hashing)
{
  # Data destination
  out_dir <- file.path(paste(where, folder, sep="/"))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  # Saving 10X tables (matrix, barcodes, genes) + metadata
  Matrix::writeMM(obj=object, file=paste(out_dir, "matrix.mtx", sep="/"))
  gzip(paste(out_dir, "matrix.mtx", sep="/"))
  write.table(colnames(object), sep="\t", file=gzfile(paste(out_dir, "barcodes.tsv.gz", sep="/")),
              quote=F, row.names=F, col.names=F)
  #write.table(rownames(object), sep="\t", file=gzfile(paste(out_dir, "features.tsv.gz", sep="/")),
  #            quote=F, row.names=F, col.names=F)
  write.table(metadata, sep="\t", file=gzfile(paste(out_dir, "metadata.tsv.gz", sep="/")),
              quote=F, row.names=T, col.names=T)
  
  gex.source.file <- paste(gex.data.dir, 'features.tsv.gz', sep='/')
  adt.source.file <- paste(adt.data.dir, 'features.tsv.gz', sep='/')
  hto.source.file <- paste(hto.data.dir, 'features.tsv.gz', sep='/')

  has.adt     <- exists('adt.source.file') && !is.null(adt.source.file) && file.exists(adt.source.file)
  has.hashing <- exists('hto.source.file') && !is.null(hto.source.file) && file.exists(hto.source.file)
  
  features.dest.file <- paste(out_dir, 'features.tsv', sep='/')
  if (isTRUE(remove_hashing)) {
    remove_hashing_string <- " | grep -vP 'Hash'"
  } else {
    remove_hashing_string <- ""
  }
  exec <- paste0("gunzip -c ", gex.source.file, remove_hashing_string, " > ", features.dest.file)
  if (isTRUE(has.adt)) {
    exec <- paste0(exec, " && gunzip -c ", adt.source.file, remove_hashing_string, " >> ", features.dest.file)
  }
  if (isTRUE(has.hashing) && identical(FALSE, remove_hashing) && adt.source.file != hto.source.file) {
    exec <- paste0(exec, " && gunzip -c ", hto.source.file, remove_hashing_string, " >> ", features.dest.file)
  }
  exec <- paste0(exec, " && gzip -f ", features.dest.file)
  system(exec)
}


# Parse arguments
option_list <- list(
  make_option(c("--project_id"), default=NULL,
              help = "Project ID to be used to name the main Seurat object."),
  make_option(c("--gex_library_id"), default=NULL,
              help = "The Library ID of the Gene Expression sequenced sample (re-sequenced samples will have the same Library ID)."),
  make_option(c("--gex_sequencing_id"), default=NULL,
              help = "Sequencing ID of the Gene Expression sequenced data. This ID will be used as a prefix when splitting by the hash tags."),
  make_option(c("--gex_data_dir"), default=NULL,
              help = "Path to the Gene Expression data directory."),
  make_option(c("--hto_library_id"), default=NULL,
              help = "The Library ID of the Hashing sequenced sample (re-sequenced samples will have the same Library ID)."),
  make_option(c("--hto_sequencing_id"), default=NULL,
              help = "Sequencing ID of the Hashing sequenced data."),
  make_option(c("--hto_data_dir"), default=NULL,
              help = "Path to the Hashing data directory."),
  make_option(c("--adt_library_id"), default=NULL,
              help = "The Library ID of the Protein Expression sequenced sample (re-sequenced samples will have the same Library ID)."),
  make_option(c("--adt_sequencing_id"), default=NULL,
              help = "Sequencing ID of the Protein Expression sequenced data."),
  make_option(c("--adt_data_dir"), default=NULL,
              help = "Path to the Antibody Expression data directory."),
  make_option(c("--output_dir"), default=NULL,
              help = "Path to the output directory."),
  make_option(c("--demux_method"), default=NULL,
              help = "Demultiplexing method to use for Hashing data. Available methods are: htodemux or multiseq."),
  make_option(c("--demux_ms_autothresh"), default=TRUE,
              help = "When MULTIseq method is chosen for demultiplexing, enable/disable autoThresh option."),
  make_option(c("--save_full_matrix"), default=TRUE,
              help = "Save the matrix after loading all samples, without any filtering."),
  make_option(c("--save_splitted_matrices"), default=FALSE,
              help = "Save one matrix per hashtag."),
  make_option(c("--save_filtered_matrix"), default=TRUE,
              help = "Save a unique matrix containing only the singlet barcodes after demultiplexing.")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Validating arguments
if (is.null(opt$project_id)) {
  stop(paste0("Argument not provided. --project_id"))
}
if (!is.null(opt$gex_library_id) && length(grep("[^a-zA-Z0-9_-]", opt$gex_library_id)) > 0) {
  stop(paste0("Invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --gex_library_id ", opt$gex_library_id))
}
if (is.null(opt$gex_sequencing_id) || length(grep("[^a-zA-Z0-9_-]", opt$gex_sequencing_id)) > 0) {
  stop(paste0("Argument not provided or invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --gex_sequencing_id"))
}
if (is.null(opt$gex_data_dir) || !dir.exists(opt$gex_data_dir)) {
  stop(paste0("Argument not provided or directory provided does not exist. --gex_data_dir ", opt$gex_data_dir))
}

if (!is.null(opt$hto_library_id) && length(grep("[^a-zA-Z0-9_-]", opt$hto_library_id)) > 0) {
  stop(paste0("Invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --hto_library_id ", opt$hto_library_id))
}
if (!is.null(opt$hto_sequencing_id) && length(grep("[^a-zA-Z0-9_-]", opt$hto_sequencing_id)) > 0) {
  stop(paste0("Invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --hto_sequencing_id"))
}
if (!is.null(opt$hto_data_dir) && !dir.exists(opt$hto_data_dir)) {
  stop(paste0("Directory provided does not exist. --hto_data_dir ", opt$hto_data_dir))
}

if (!is.null(opt$adt_library_id) && length(grep("[^a-zA-Z0-9_-]", opt$adt_library_id)) > 0) {
  stop(paste0("Invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --adt_library_id ", opt$adt_library_id))
}
if (!is.null(opt$adt_sequencing_id) && length(grep("[^a-zA-Z0-9_-]", opt$adt_sequencing_id)) > 0) {
  stop(paste0("Invalid characters detected. This option can only contain letters, numbers, dash (-) and underscore (_). --adt_sequencing_id"))
}
if (!is.null(opt$adt_data_dir) && !dir.exists(opt$adt_data_dir)) {
  stop(paste0("Directory provided does not exist. --adt_data_dir ", opt$adt_data_dir))
}

if (is.null(opt$output_dir)) {
  stop(paste0("Argument not provided. --output_dir ", opt$output_dir))
}
if (!is.null(opt$demux_method) && !(opt$demux_method %in% c("htodemux", "multiseq"))) {
  stop(paste0("Invalid demultiplexing option. Please, choose between 'htodemux' or 'multiseq'. --demux_method"))
}

project.id <- opt$project_id
gex.library.id <- opt$gex_library_id
gex.sequencing.id <- opt$gex_sequencing_id
gex.data.dir <- opt$gex_data_dir
hto.library.id <- opt$hto_library_id
hto.sequencing.id <- opt$hto_sequencing_id
hto.data.dir <- opt$hto_data_dir
adt.library.id <- opt$adt_library_id
adt.sequencing.id <- opt$adt_sequencing_id
adt.data.dir <- opt$adt_data_dir
output.dir <- opt$output_dir
demux.method <- opt$demux_method
demux.ms.autothresh <- opt$demux_ms_autothresh
save.full.matrix <- opt$save_full_matrix
save.splitted.matrices <- opt$save_splitted_matrices
save.filtered.matrix <- opt$save_filtered_matrix
            

# Testing data values
#
# Hash and Cite
# -------------
# project.id          <- 'P220447'
# gex.library.id      <- 'DAK10037A1'
# gex.sequencing.id   <- 'DAK10037A1'
# gex.data.dir  <- '/well/singlecell/P220447/10X-gex-grouped/DAK10037A1/outs/filtered_feature_bc_matrix'
# hto.library.id      <- 'DAK10037A3'
# hto.sequencing.id   <- 'DAK10037A3'
# hto.data.dir    <- '/well/singlecell/P220447/10X-feat-grouped/DAK10037A3/outs/raw_feature_bc_matrix'
# adt.library.id     <- 'DAK10037A3'
# adt.sequencing.id   <- 'DAK10037A3'
# adt.data.dir    <- '/well/singlecell/P220447/10X-feat-grouped/DAK10037A3/outs/raw_feature_bc_matrix'
# output.dir         <- '/well/singlecell/P220447/10X-aggregate'
# demux.method        <- 'htodemux'

# project.id <- "P190726"
# gex.library.id <- "FOR7865A1"
# gex.sequencing.id <- "747496_GX69"
# hto.library.id <- "FOR7865A7"
# hto.sequencing.id <- "747496_NX05"
# adt.library.id <- "FOR7865A7"
# adt.sequencing.id <- "747496_NX05"
# output.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-features")
# gex.data.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-count/", gex.sequencing.id,"/outs/filtered_feature_bc_matrix")
# hto.data.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-count/", hto.sequencing.id,"/outs/raw_feature_bc_matrix")
# adt.data.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-count/", adt.sequencing.id,"/outs/raw_feature_bc_matrix")
# demux.method <- "htodemux"
# save.full.matrix <- TRUE
# save.splitted.matrices <- TRUE
# save.filtered.matrix <- TRUE
#
# Only Cite
# -------------
# project.id <- "P190730"
# gex.library.id <- "BAS7763A17"
# gex.sequencing.id <- "750205_GX68"
# hto.library.id <- ""
# hto.sequencing.id <- ""
# adt.library.id <- "BAS7763A24"
# adt.sequencing.id <- "750205_NX61"
# output.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-features")
# gex.data.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-count/", gex.sequencing.id,"/outs/filtered_feature_bc_matrix")
# hto.data.dir <- ""
# adt.data.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-count/", adt.sequencing.id,"/outs/raw_feature_bc_matrix")
# demux.method <- ""
# save.full.matrix <- TRUE
# save.splitted.matrices <- TRUE
# save.filtered.matrix <- TRUE
#
# Only Hash
# -------------
# project.id <- "P190800"
# gex.library.id <- "TOF7909A2"
# gex.sequencing.id <- "750206_GX24"
# hto.library.id <- "TOF7909A1"
# hto.sequencing.id <- "750206_GX01"
# adt.library.id <- ""
# adt.sequencing.id <- ""
# output.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-features")
# gex.data.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-count/", gex.sequencing.id,"/outs/filtered_feature_bc_matrix")
# hto.data.dir <- paste0("/data/borrar/hashing.reports/", project.id,"/10X-count/", hto.sequencing.id,"/outs/raw_feature_bc_matrix")
# adt.data.dir <- ""
# demux.method <- "htodemux"
# save.full.matrix <- TRUE
# save.splitted.matrices <- TRUE
# save.filtered.matrix <- TRUE

# Create output directory
output.dir <- paste(output.dir, gex.sequencing.id, sep="/")
dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
# Create plots directory
plots.dir <- file.path(paste(output.dir, "plots", sep="/"))
dir.create(plots.dir, showWarnings=FALSE)
# Create objects directory
objects.dir <- file.path(paste(output.dir, "objects", sep="/"))
dir.create(objects.dir, showWarnings=FALSE)

# Define whether project has hashing and/or antibody data
has.hashing <- exists('hto.data.dir') && !is.null(hto.data.dir) && dir.exists(hto.data.dir)
has.adt     <- exists('adt.data.dir') && !is.null(adt.data.dir) && dir.exists(adt.data.dir)

# Read GEX 10X matrix
gex.data.counts <- Read10X(data.dir=gex.data.dir)

# Check if it's a multi-feature sample
if (class(gex.data.counts) == "list" && "Gene Expression" %in% names(gex.data.counts)) {
  gex.data.counts <- gex.data.counts$`Gene Expression`
}

# TODO: filter cells before proceeding!


# When there is a hashing sample, intersect the cell barcodes before
# creating the RNA Seurat object.
if (isTRUE(has.hashing)) {
  # Read 10X matrix for Hashing
  hto.data.counts <- Read10X(data.dir=hto.data.dir)

  # Check if it's a multi-feature sample
  if (class(hto.data.counts) == "list" && "Antibody Hashing" %in% names(hto.data.counts)) {
    hto.data.counts <- hto.data.counts$`Antibody Hashing`
  }
  
  # Remove the trailing string after the feature name (like _TotalSeqC, etc.)
  rownames(hto.data.counts) <- gsub("[_-].*", "", rownames(hto.data.counts))
  
  # Select cell barcodes detected by both RNA and HTO
  joint.bcs <- intersect(colnames(gex.data.counts), colnames(hto.data.counts))
  
  # Subset RNA and HTO counts by joint cell barcodes
  gex.data.counts <- gex.data.counts[, joint.bcs]
  hto.data.counts <- hto.data.counts[, joint.bcs]
}


# Define the RNA assay, and store raw counts for it
data <- CreateSeuratObject(counts=gex.data.counts, project=project.id)


# Add the gex.sequencing.id as a metadata field
meta.gex.sequencing.id <- rep(gex.sequencing.id, length(data@meta.data$orig.ident))
names(meta.gex.sequencing.id) <- rownames(data@meta.data)
data$gex.sequencing.id <- meta.gex.sequencing.id
rm("meta.gex.sequencing.id")

# Add the gex.library.id as a metadata field
meta.gex.library.id <- rep(gex.library.id, length(data@meta.data$orig.ident))
names(meta.gex.library.id) <- rownames(data@meta.data)
data$gex.library.id <- meta.gex.library.id
rm("meta.gex.library.id")

# Get mitochondrial gene content and add it as metadata
mito.genes <- grep(pattern="^MT-", x=rownames(data), ignore.case=TRUE, value=TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object=data, slot="counts")[mito.genes, ]) / Matrix::colSums(GetAssayData(object=data, slot="counts"))
data$percent.mito <- percent.mito
rm("mito.genes", "percent.mito")

# Get ribosomal gene content and add it as metadata
ribo.genes <- grep(pattern="^RP", x=rownames(data), ignore.case=TRUE, value=TRUE)
percent.ribo <- Matrix::colSums(GetAssayData(object=data, slot="counts")[ribo.genes, ]) / Matrix::colSums(GetAssayData(object=data, slot="counts"))
data$percent.ribo <- percent.ribo
rm("ribo.genes", "percent.ribo")

# When there is a features sample, match its cell barcodes against
# the RNA ones, before creating the ADT Seurat object.
if (isTRUE(has.adt)) {
  
  # Add the adt.sequencing.id as a metadata field
  meta.adt.sequencing.id <- rep(adt.sequencing.id, length(data@meta.data$orig.ident))
  names(meta.adt.sequencing.id) <- rownames(data@meta.data)
  data$adt.sequencing.id <- meta.adt.sequencing.id
  rm("meta.adt.sequencing.id")
  
  # Add the adt.library.id as a metadata field
  meta.adt.library.id <- rep(adt.library.id, length(data@meta.data$orig.ident))
  names(meta.adt.library.id) <- rownames(data@meta.data)
  data$adt.library.id <- meta.adt.library.id
  rm("meta.adt.library.id")

  # Read 10X matrix for ADT
  adt.data.counts <- Read10X(data.dir=adt.data.dir)

  # Check if it's a multi-feature sample
  if (class(adt.data.counts) == "list" && "Antibody Capture" %in% names(adt.data.counts)) {
    adt.data.counts <- adt.data.counts$`Antibody Capture`
  }
  
  # Remove the trailing string after the feature name (like _TotalSeqC, etc.)
  rownames(adt.data.counts) <- gsub("_.*", "", rownames(adt.data.counts))
  
  # Select cell barcodes detected by both RNA and ADT
  joint.bcs <- intersect(colnames(gex.data.counts), colnames(adt.data.counts))
  
  # Subset ADT counts by joint cell barcodes
  adt.data.counts <- adt.data.counts[, joint.bcs]
  
  # Define the ADT assay, and store raw counts for it
  data[["ADT"]] <- CreateAssayObject(counts=adt.data.counts)
}


# Process the hashing data
if (isTRUE(has.hashing)) {
  
  # Add the hto.sequencing.id as a metadata field
  meta.hto.sequencing.id <- rep(hto.sequencing.id, length(data@meta.data$orig.ident))
  names(meta.hto.sequencing.id) <- rownames(data@meta.data)
  data$hto.sequencing.id <- meta.hto.sequencing.id
  rm("meta.hto.sequencing.id")
  
  # Add the hto.library.id as a metadata field
  meta.hto.library.id <- rep(hto.library.id, length(data@meta.data$orig.ident))
  names(meta.hto.library.id) <- rownames(data@meta.data)
  data$hto.library.id <- meta.hto.library.id
  rm("meta.hto.library.id")
  
  # Remove any hashtag with less than 1000 reads to avoid HTODemux crashing
  # to avoid "Cells with zero counts exist as a cluster." error [Isar Solution]
  # hto.data.counts <- hto.data.counts[rowSums(hto.data.counts)>as.numeric(summary(rowSums(hto.data.counts))[1]), ]
  
  # to avoid "Cells with zero counts exist as a cluster." error [Fabeola Solutions]
  # hto.data.counts <- hto.data.counts[rowSums(hto.data.counts)>1000, ]
  
  # Remove any hashtag with less than 5% of the mean total reads to avoid HTODemux crashing
  #rowSums(hto.data.counts) / mean(rowSums(hto.data.counts)) > 0.05 -> keep
  #hto.data.counts <- hto.data.counts[keep,]

  # Define the HTO assay, and store raw counts for it
  data[["HTO"]] <- CreateAssayObject(counts = hto.data.counts)
  data <- NormalizeData(data, assay="HTO", normalization.method="CLR", display.progress=FALSE)
  data <- ScaleData(data, features = VariableFeatures(data))
  
  if (demux.method == 'multiseq') {

    data <- MULTIseqDemux(data, assay="HTO", autoThresh=TRUE, maxiter=10)
    
    # Create the missing metadata variables to make it compatible with HTODemux
    ## hash.ID
    data@meta.data$hash.ID <- data@meta.data$MULTI_ID
    ## HTO_classification
    data@meta.data$HTO_classification <- data@meta.data$MULTI_classification
    ## HTO_classification.global
    data@meta.data$HTO_classification.global <- as.character(data@meta.data$HTO_classification)
    data@meta.data$HTO_classification.global[grepl("_", data@meta.data$HTO_classification.global)] <- "Doublet"
    data@meta.data$HTO_classification.global[!(data@meta.data$HTO_classification.global %in% c("Doublet","Negative"))] <- "Singlet"
    data@meta.data$HTO_classification.global <- factor(data@meta.data$HTO_classification.global)
    data@meta.data$HTO_maxID <- gsub('_.*','',as.character(data@meta.data$HTO_classification))    
    
    # Remove the previous objects
    data@meta.data$MULTI_ID <- NULL
    data@meta.data$MULTI_classification <- NULL
    
    print(kable(table(data$HTO_classification.global), col.names = c("Cells Classification", "Total")))
    print(kable(table(data$hash.ID), col.names = c("Classification", "Cells Count")))
    
  } else if (demux.method == 'htodemux') {
    data <- HTODemux(data, assay="HTO", positive.quantile=0.99, kfunc="clara", verbose=TRUE)
    
    # Find and scale variable features
    pbmc.hashtag <- HTODemux(data, assay = "HTO", positive.quantile = 0.99, verbose=TRUE)

    library("knitr")
    print(kable(table(data$HTO_classification.global), col.names = c("Cells Classification", "Total")))
    print(kable(table(data$hash.ID), col.names = c("Classification", "Cells Count")))
    
  }
  #print(kable(data@meta.data %>% group_by(HTO_classification.global) %>% summarise(Total=n()), caption="Cells Classification"))
  #print(kable(data@meta.data %>% group_by(hash.ID) %>% summarise(Total=n()), caption="Cell Counts"))

  # Heatmap
  # if (nrow(data@meta.data) > 5000) {
  #   plot = HTOHeatmap(data, assay="HTO", ncells=5000)
  # } else {
  #   plot = HTOHeatmap(data, assay="HTO", ncells=nrow(data@meta.data)-1)
  # }
  # ggsave(paste(plots.dir, "HTO_HTOHeatmap.png", sep="/"), plot=plot, height=2*plot.size, width=2*plot.size, dpi=plot.dpi, units=plot.units)
  # 
  # RidgePlot on hash.ID, it includes the categories Doublet/Negative
  # It dropped plotting on HTO_maxID, because it forces an identity to
  # doublets/negatives.
  # hto.features <- rownames(x=data[["HTO"]])
  # num.rows <- round(length(hto.features) / 2, 0)
  # plot = RidgePlot(object=data, assay="HTO", features=hto.features, ncol=2)
  # ggsave(paste(plots.dir, "HTO_RidgePlot.png", sep="/"), plot=plot, height=num.rows*plot.size, width=2*plot.size, dpi=plot.dpi, units=plot.units)
  # 
  # 
  # # Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
  # for (i in 1:(length(hto.features)-1)) {
  #   for (j in (i+1):length(hto.features)) {
  #     plot = FeatureScatter(object=data, feature1=paste0("hto_", hto.features[i]), feature2=paste0("hto_", hto.features[j]))
  #     ggsave(paste(plots.dir, paste0("HTO_FeatureScatter_", hto.features[i], "_", hto.features[j], ".png"), sep="/"), plot=plot, height=2*plot.size, width=2*plot.size, dpi=plot.dpi, units=plot.units)
  #   }
  }

# Save the merged raw data matrix and data object
save.full.matrix <- TRUE
if (isTRUE(save.full.matrix)) {
  # Save sparse matrix
  raw.dir <- "raw_matrix"
  cells.meta <- data@meta.data
  cells.data <- gex.data.counts
  if (has.adt) {
    cells.data <- rbind(cells.data, adt.data.counts)
  }
  if (isTRUE(has.hashing)) {
    cells.data <- rbind(cells.data, hto.data.counts)
  }
  saveHashOutput(objects.dir, raw.dir, cells.data, cells.meta, gex.data.dir, adt.data.dir, hto.data.dir, FALSE)
  # Save RDS object
  saveRDS(cells.data, file=paste(objects.dir, raw.dir, "seurat.rds", sep="/"))
  # Remove temporary objects
  rm("cells.meta", "cells.data")
}

  # Print some metrics per HTO antibody
  print(as.matrix(data@assays$HTO@counts) %>% rowSums(.))
  print(as.matrix(data@assays$HTO@counts) %>% rowMeans(.))
  print(apply(as.matrix(data@assays$HTO@counts), 1, sd))
  print(apply(as.matrix(data@assays$HTO@counts), 1, summary))
  
  # Remove data object, because it will be re-created based on the splitted matrices
  rm("data")

# Remove data.counts objects, because they are not needed any more
suppressWarnings(rm("gex.data.counts", "hto.data.counts", "adt.data.counts"))

