# PROJECT 2 ###

# Yagmur Akarsu 
# Lizzie Schmitz

# 1.1 Setting up the environment #

set.seed(42) # set seed

install.packages("devtools")
BiocManager::install(c("rhdf5", "motifmatchr", "chromVAR", "Rsamtools")) # install dependenies dirst 
devtools::install_github("GreenleafLab/ArchR", ref = "master")
install.packages("Cairo")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("gridExtra")

library(gridExtra)
library(ggplot2)
library(Cairo)
library(ArchR)
library(patchwork)

addArchRThreads(threads = 1) # number of threads should be 1/2 or 3/4 of the total available cores, 16 is default, 1 for window OS

addArchRGenome("hg38")

#1.2 

# Defining our directory to get our data.İ
setwd("c:/Users/ArchR")
list.files("./scbi_p2/")

# Paths to your fragment files
frag_files <- c(
  "./scbi_p2/hft_ctx_w21_dc1r3_r1_atac_fragments.tsv.gz",
  "./scbi_p2/hft_ctx_w21_dc2r2_r1_atac_fragments.tsv.gz",
  "./scbi_p2/hft_ctx_w21_dc2r2_r2_atac_fragments.tsv.gz"
)

# create arrow files for the three samples

ArrowFiles <- createArrowFiles(
  inputFiles = frag_files,               # Use fragment file paths
  sampleNames = c("hft_ctx_w21_dc1r3_r1_atac_fragments", "hft_ctx_w21_dc2r2_r1_atac_fragments", "hft_ctx_w21_dc2r2_r2_atac_fragments"),       # Use the names of the files
  minTSS = 4,                            # Replace filterTSS
  minFrags = 500,                        # Replace filterFrags
  addGeneScoreMat = TRUE,                
  addTileMat = TRUE,                     # Add tiling matrix
  TileMatParams = list(tileSize = 1000), # Set tile size
  verbose = TRUE
)

ArrowFiles <- c("hft_ctx_w21_dc1r3_r1_atac_fragments.arrow", "hft_ctx_w21_dc2r2_r1_atac_fragments.arrow", "hft_ctx_w21_dc2r2_r2_atac_fragments.arrow")
ArrowFiles # check to see if data loaded properly 

# Create ArchR project
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "ArchR_Project")

### IGNORE THIS IF U HAVE A MAC #### 
# Retrieve all cell metadata
allCellData <- getCellColData(proj, select = c("Sample", "nFrags", "TSSEnrichment"))
# Generate log10(nFrags) vs. TSSEnrichment plots
log10nFragsPlots <- lapply(ArrowFiles, function(file) {
  # Extract sample name
  sampleName <- gsub("\\.arrow$", "", basename(file))
  # Subset data for the current sample
  cellData <- allCellData[allCellData$Sample == sampleName, ]
  
  # Add log10nFrags column
  cellData$log10nFrags <- log10(cellData$nFrags)
  
  # Create scatter plot with density gradient
  ggplot(cellData, aes(x = log10nFrags, y = TSSEnrichment)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(
      x = "log10(Unique Nuclear Fragments)",
      y = "TSSEnrichment",
      title = paste("log10(nFrags) vs. TSSEnrichment\n", sampleName)
    ) +
    theme(plot.title = element_text(size = 10))
})

# Combine log10nFrags plots for display
combinedLog10nFragsPlot <- wrap_plots(log10nFragsPlots, ncol = 1)
print(combinedLog10nFragsPlot)

# Generate fragment size distribution plots
fragSizePlots <- lapply(ArrowFiles, function(file) {
  # Plot the fragment sizes for each Arrow file
  plotFragmentSizes(proj) + 
    ggtitle(paste("Fragment Size Distribution\n", basename(file)))
})

# Combine fragment size distribution plots for display
combinedFragSizePlot <- wrap_plots(fragSizePlots, ncol = 1)
print(combinedFragSizePlot)

# Display combined plots
print(combinedLog10nFragsPlot)
print(combinedFragSizePlot)

# Save plots
ggsave("Log10nFrags_vs_TSSEnrichment_AllSamples.pdf", plot = combinedLog10nFragsPlot, height = 12, width = 6)
ggsave("FragmentSize_Distribution_AllSamples.pdf", plot = combinedFragSizePlot, height = 12, width = 6)
### IGNORE THIS IF U HAVE A MAC #### 

# 1.3) Identify Doublets

# Add doublet scores to the project
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,            # Number of nearest neighbors to consider for pseudo-doublets
  knnMethod = "UMAP", # Embedding method for nearest neighbor search
  LSIMethod = 1      # Method for inferring doublets (set to 1 as default)
)

# Visualize the doublet scores
plotEmbedding(proj, colorBy = "cellColData", name = "DoubletScore")

# Set threshold for doublets
proj <- proj[!proj$DoubletScore > 0.25,]

# Check the number of cells after removal
ncol(proj)

# Number of doublets detected before and after removal
sum(proj$DoubletScore > 0.25)  # Number of detected doublets
ncol(proj)  # Number of remaining cells
