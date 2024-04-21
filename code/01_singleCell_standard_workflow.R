here::i_am("code/01_singleCell_standard_workflow.R")

#script to perform standard workflow steps to analyze single cell RNA-Seq data
#data: 10k Human PBMCs, Multiome v1.0, Chromium Controller

#load libraries
library(Seurat)
library(tidyverse)

#load the NSCLC dataset