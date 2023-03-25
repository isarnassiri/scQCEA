
# module add R-bundle-Bioconductor/3.12-foss-2020b-R-4.0.3

suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))

#---------------------------------- Parse arguments
option_list <- list(
  make_option(c("--project_id"), default=NULL,
              help = "Project ID to be used to name the main Seurat object."),
  make_option(c("--backend.data.dir"), default=NULL,
              help = "Path to the backend directory.")
)
opt <- parse_args(OptionParser(option_list=option_list))


#---------------------------------- example input
project_id <- "P220453"
backend.data.dir <- "/well/singlecell/references/reference_gene_sets/"
SamlesMetadata = fread(paste0('/well/singlecell/P220453/samples.metadata'), stringsAsFactors = F, header = T)
#----------------------------------

project_id <- opt$project_id
backend.data.dir <- opt$backend.data.dir
#-----------------

SamlesMetadata = fread(paste0('/well/singlecell/', project_id, '/samples.metadata'), stringsAsFactors = F, header = T)
print( paste0("Used Genome: ", unique(SamlesMetadata$Genome)) )
GenomeType <- unique(SamlesMetadata$Genome)

if(GenomeType == 'GRCm38' | GenomeType == 'GRCm38-premrna'){organism = 'mmusculus'}else if(GenomeType == 'GRCh38' | GenomeType == 'GRCh38-premrna'){organism = 'hsapiens'}else{organism = 'others'}

setwd(paste0('/well/singlecell/', project_id))

list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

#list files recursive up to a certain level in R
dirs <- list.dirs.depth.n(".", n = 2)

#Include pattern in list.dirs
dirs <- grep("10X-gex", dirs, value = TRUE)

if(length(grep("^./10X-gex$", dirs, value = TRUE))>0){dirs <- dirs[-which(grepl("^./10X-gex$", dirs))]}
if(length(grep("Inputs", dirs, value = TRUE))>0){dirs <- dirs[-which(grepl("Inputs", dirs))]}
if(length(grep("^./10X-gex-grouped$", dirs, value = TRUE))>0){dirs <- dirs[-which(grepl("^./10X-gex-grouped$", dirs))]}

OP = if(organism == "hsapiens"){paste0(backend.data.dir, 'human')}else if(organism == "mmusculus"){paste0(backend.data.dir, 'mouse')}else{'others'}

if(length(grep('/10X-gex-grouped', dirs))>0 & OP != 'others')
{
  input.dir = dirs[grep('/10X-gex-grouped/', dirs)]
  input.dir = paste0('/well/singlecell/', project_id, '/', input.dir)
  input.dir = gsub('/./', '/', input.dir)
  
  guide_file = data.frame(
    project_id = replicate(length(input.dir),project_id),
    sample_id = gsub('.*\\/','',input.dir),
    inpute_dir = input.dir,
    backend_data_dir = replicate(length(input.dir),OP),
    organism = replicate(length(input.dir),organism))
  
  write.table(guide_file, 'gex_grouped_aggregation', quote = F, col.names = F, row.names = F, sep = '\t')
}

if(length(grep('/10X-gex', dirs))>0 & OP != 'others')
{
  input.dir = dirs[grep('/10X-gex/', dirs)]
  input.dir = paste0('/well/singlecell/', project_id, '/', input.dir)
  input.dir = gsub('/./', '/', input.dir)
  
  guide_file = data.frame(
    project_id = replicate(length(input.dir),project_id),
    sample_id = gsub('.*\\/','',input.dir),
    inpute_dir = input.dir,
    backend_data_dir = replicate(length(input.dir),OP),
    organism = replicate(length(input.dir),organism))
  
  write.table(guide_file, 'gex_ungrouped_aggregation', quote = F, col.names = F, row.names = F, sep = '\t')
}

if(length(grep('/10X-gex', dirs))==0 | OP == 'others')
{
  write.table("others", 'gex_others_aggregation', quote = F, col.names = F, row.names = F, sep = '\t')
}
  
