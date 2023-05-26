
#############################################################################################
###################################### Generate Interactive QC Report #######################
#############################################################################################

#'@import rstudioapi
#'@import devtools
#'@import zip
#'@import DiagrammeR
#'@import stringr
#'@import bsselectR
#'@import kableExtra
#'@import DT
#'@import downloadthis
#'@import ggplot2
#'@import data.table
#'@import dplyr
#'@import readr
#'@import rmarkdown
#'@import R.utils
#'@export
#'@name GenerateInteractiveQCReport
#'@title Generate Interactive QC Report
#'@description An interactive QC report automatically will be generated in one HTML file, including four sections: experimental workflow, data processing workflow, sample information and QC metrics, data analysis and quality control.
#'@author {Isar Nassiri}
#'@param InputDir
#'Path of a folder including input files, called Inputs/
#'@return You can find the results in R object under title of 'RESULTsGenomicFeatures' and 'RESULTsChromatinState'.
#'@examples
#'library("scQCEA")
#'InputDir=system.file("extdata", package = "scQCEA")
#'GenerateInteractiveQCReport()
#'@export
GenerateInteractiveQCReport <- NULL
GenerateInteractiveQCReport <- function(InputDir)
{
  
  setwd(InputDir)
  SamplesMetadata = fread('Inputs/samples.metadata', stringsAsFactors = FALSE, header = TRUE);
  
  #- For human or mouse input files we run enrichment analysis
  if(length(grep('GRCh38|GRCm38', SamplesMetadata$Genome)) > 0)
  {
    invisible(file.copy(paste0(system.file("extdata", package = "scQCEA"),'/RMarkDownR/','RMarkDown.Rmd'), InputDir ))
    invisible(file.copy(paste0(system.file("extdata", package = "scQCEA"),'/RMarkDownR/RMarkDown/'), InputDir, recursive=TRUE))
  }else{
    invisible(file.copy(paste0(system.file("extdata", package = "scQCEA"),'/RMarkDownR/','RMarkDown-nonHumanMouse.Rmd'), InputDir ))
    invisible(file.rename(paste0(InputDir, '/RMarkDown-nonHumanMouse.Rmd'), paste0(InputDir, '/RMarkDown.Rmd')))
    invisible(file.copy(paste0(system.file("extdata", package = "scQCEA"),'/RMarkDownR/RMarkDown/'), InputDir, recursive=TRUE))
  }

  setwd(InputDir)
  render("RMarkDown.Rmd", quiet = TRUE);
  
  invisible(file.rename('RMarkDown.html', 'CLICK_ME.html'));
  invisible(file.copy(from = 'CLICK_ME.html', to = "Outputs/CLICK_ME.html"));
  invisible(file.remove(paste0(getwd(), '/CLICK_ME.html')));
  
  #--- Delete extra files
  setwd(paste0(InputDir, "/Outputs")); 

  list.dirs.depth.n <- function(p, n) {
    res <- list.dirs(p, recursive = FALSE)
    if (n > 1) {
      add <- list.dirs.depth.n(res, n-1)
      c(res, add)
    } else {
      res
    }
  }
  
  listdir = list.dirs.depth.n(".", n = 5)
  toMatch <- c('analysis', 'raw_feature_bc_matrix')
  matches <- unique(grep(paste(toMatch,collapse="|"), listdir, value=TRUE))
  invisible(unlink(matches, recursive = TRUE, force = TRUE));
  
  listfiles = list.files(pattern = 'read_count.csv', recursive = TRUE)
  invisible(unlink(listfiles, recursive = TRUE, force = TRUE));
  
  listfiles = list.files(pattern = '*.pdf', recursive = TRUE)
  invisible(unlink(listfiles, recursive = TRUE, force = TRUE));
  
  #--- generate zip file
  setwd(InputDir); 
  PInf = fread('Inputs/PInf.txt', stringsAsFactors = FALSE, header = FALSE);
  
  setwd(paste0(InputDir, "/Outputs"));
  zip(zipfile = paste0('OGC_Interactive_QC_Report_', gsub('.*=','',PInf$V1[1]),'.zip'), files = c('CLICK_ME.html', 'Inputs'), recurse = TRUE, include_directories = TRUE);
  
  invisible(file.remove(paste0(InputDir,'/','RMarkDown.Rmd')))
  invisible(unlink(paste0(InputDir,'/RMarkDown/'), recursive=TRUE))
  
  # ----------------------------------
  cat(paste0("\033[0;", 47, "m", "You can find the Interactive QC Report in: ", "\033[0m","\n", paste0(InputDir, '/Outputs/') ))
  
  # date()
  # sessionInfo()
}

