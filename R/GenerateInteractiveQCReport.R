
#############################################################################################
###################################### Generate Interactive QC Report #######################
#############################################################################################
#'@export
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
#'@return You can find the results in R object under title of 'RESULTsGenomicFeatures' and 'RESULTsChromatinState'.
#'@examples
#'library("scQCEA")
#'GenerateInteractiveQCReport()
#'@export
GenerateInteractiveQCReport <- NULL
GenerateInteractiveQCReport <- function()
{
  setwd(system.file("extdata", package = "scQCEA")); 
  render("SourceCode.Rmd", quiet = T);
  
  invisible(file.rename('SourceCode.html', 'CLICK_ME.html'));
  invisible(file.copy(from = 'CLICK_ME.html', to = "Outputs/CLICK_ME.html"));
  invisible(file.remove(paste0(getwd(), '/CLICK_ME.html')));
  
  #--- Delete extra files
  setwd(system.file("extdata/Outputs", package = "scQCEA")); 

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
  invisible(unlink(matches, recursive = T, force = T));
  
  listfiles = list.files(pattern = 'read_count.csv', recursive = T)
  invisible(unlink(listfiles, recursive = T, force = T));
  
  listfiles = list.files(pattern = '*.pdf', recursive = T)
  invisible(unlink(listfiles, recursive = T, force = T));
  
  #--- generate zip file
  setwd(system.file("extdata/", package = "scQCEA")); 
  PInf = fread('Inputs/PInf.txt', stringsAsFactors = F, header = F);
  
  setwd(system.file("extdata/Outputs", package = "scQCEA")); 
  zip(zipfile = paste0('OGC_Interactive_QC_Report_', gsub('.*=','',PInf$V1[1]),'.zip'), files = c('CLICK_ME.html', 'Inputs'), recurse = TRUE, include_directories = TRUE);
  
  # ----------------------------------
  cat(paste0("\033[0;", 47, "m", "You can find the Interactive QC Report in: ", "\033[0m","\n", system.file("extdata/Outputs/", package = "scQCEA")))
  
  # date()
  # sessionInfo()
}

