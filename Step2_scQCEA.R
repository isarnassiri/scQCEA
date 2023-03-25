# install the scQCEA package
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# library(devtools)
# devtools::install_github("walkerke/bsselectR")
# install_github("isarnassiri/scQCEA")

library("scQCEA")
InputDir='/Folder/Path/' # this folder includes Inputs/ folder - the output of Step1_Inputs.sh script.
GenerateInteractiveQCReport(InputDir)

Notes:
1. Run this script on your local machine.
2. You can use sftp to download Inputs/ folder and upload Outputs/*.zip file to the project folder on the server.
