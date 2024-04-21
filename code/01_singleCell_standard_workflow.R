here::i_am("code/01_singleCell_standard_workflow.R")

#script to perform standard workflow steps to analyze single cell RNA-Seq data
#data: 10k Human PBMCs, Multiome v1.0, Chromium Controller

#load libraries
library(Seurat)
library(tidyverse)

#load the NSCLC dataset
nsclc_sparse_m <- Read10X_h5(filename = here::here("raw_data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"))
#note in console: "Genome matrix has multiple modalities, returning a list of matrices for this genome"
str(nsclc_sparse_m) # see the list of modalities in the data:
  #lists 3: gene expression, antibody capture, multiplexing capture
cts <- nsclc_sparse_m$"Gene Expression"
  #set to just gene expression -> our modality of interest


#initialize the Seurat object with the raw (non-normalized data)
nsclc_seurat_obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
  #cts = counts variable we just made, project = naming, min.cells = keep features that have are expressed in at least 3 cells,
  #min.features = keep cells that have at least 200 genes/features
str(nsclc_seurat_obj)
# 29552 features across 42081 samples

# 1. QC
  #want to filter out the low quality cells and only keep high quality cells
  #to determine quality, look at number of features/genes in a cell, and number of total molecules
    #poor quality cells will have low number of genes and total molecules
  #we also look at % of mitochondrial genes
    #in low quality cells, we see a higher percentage of mitochondrial genes

#look at metadata
view(nsclc_seurat_obj@meta.data)

## % MT reads
nsclc_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(nsclc_seurat_obj, pattern = "^MT-")
  #calculating the % MT reads for all genes which start with "MT"
view(nsclc_seurat_obj@meta.data)

#visualize these features as a violin plot
VlnPlot(nsclc_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #results:
    ## alot of cells have a variety of genes
    ## a lot of cells have a high number of molecules detected
    ## a lot of cells have a high % of MT genes -> need to be filtered out

FeatureScatter(nsclc_seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
#good quality cells should have BOTH high number of genes and total molecules,
  ## x=number of transcripts, y=number of genes
  ## good quality cells will follow the straight line
  ## generally good quality, but there is a point of plateau indicating these are poor quality cells\
  ## lower right corner = very low quality (only captured same genes many times)
  ## upper left corner = very low quality (only captured a minimal amount of transcripts)

# 2. filtering
nsclc_seurat_obj <- subset(nsclc_seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                             percent.mt < 5)
  #we have now filtered down to 24708 cells

# 3. normalize data
  #need to normalize so we can compare gene expression across cells
nsclc_seurat_obj <- NormalizeData(nsclc_seurat_obj)
  #performs log normalization

  #how to go back to Seurat object
  str(nsclc_seurat_obj)
    #look at commands section; saves parameters of commands as well

# 4. identify highly variable features
nsclc_seurat_obj <- FindVariableFeatures(nsclc_seurat_obj, selection.method = "vst", nfeatures = 2000)

#identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc_seurat_obj), 10)
view(top10)






