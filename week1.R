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
install.packages("viridis")


library(viridis)
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

ArrowFiles # check to see if data loaded properly 

# 1.3) Identify Doublets

# Add doublet scores to the project
doubScores <- addDoubletScores( # adds to each arrow file
  input = ArrowFiles,
  k = 10,            # Number of nearest neighbors to consider for pseudo-doublets
  knnMethod = "UMAP", # Embedding method for nearest neighbor search
  LSIMethod = 1      # Method for inferring doublets (set to 1 as default)
)

# 1.4) Create ArchR project
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "ArchR_Project")

# 1.4) Inspect the project and cell metadata
proj

colnames(proj@cellColData) # metadata names, frags = nFrags and TSS = ReadsInTSS

# Getting peak set from data

proj_peaks <- getPeakSet(proj) #GRanges object

head(proj_peaks) # glimpse at peaks

proj <- addPeakSet(
  ArchRProj = proj,
  peakSet = proj_peaks,
  genomeAnnotation = getGenomeAnnotation(proj),
  force = TRUE # allows us to overwrite old peaks 
)

# Filter out doublets based on the threshold
proj_filtered <- proj[proj$DoubletScore < 0.25, ] # threshold = 0.25, the default, or whatever is fit

# Check the filtered project
proj_filtered

colnames(proj_filtered@cellColData) # metadata names, frags = nFrags and TSS = ReadsInTSS

# save the filtered project to a new file
saveArchRProject(proj_filtered, "c:/Users/ArchR/proj_doubletfiltered")

#### This is data after QC filtering and doublet removal ## 
# 1.4)  Number of cells in the project
num_cells <- sum(table(proj_filtered@cellColData$Sample))
cat("Number of cells in the project: \n")
print(num_cells)

# 1.4) Median TSS value and median number of fragments
median_tss <- median(proj_filtered@cellColData$ReadsInTSS)
median_frags <- median(proj_filtered@cellColData$nFrags)
cat("Median TSS value:", median_tss, "\n")
cat("Median number of fragments:", median_frags, "\n")

# 1.4) Dimensions of the peak set (number of peaks)
num_peaks <- nrow(proj_filtered@peakSet)
cat("Number of peaks in the dataset:", num_peaks, "\n")

# 1.5) Number of cells per sample
cells_per_sample <- table(proj_filtered@cellColData$Sample)
cat("Number of cells per sample:\n")
print(cells_per_sample)

# 1.5) Plot the fragment length distribution of all samples in one plot
plotFragmentLengths(proj_filtered, 
                    colorBy = "Sample", 
                    pal = viridis(100), 
                    plotTitle = "Fragment Length Distribution"
                    )

# 1.5)´ Plot the distribution of TSS enrichment scores for each sample
plotTSSEnrichment(
  ArchRProj = proj_filtered,
  groupBy = "Sample",
  chromSizes = getChromSizes(proj),
  TSS = getTSS(proj),
  pal = viridis(100),
  returnDF = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("plotTSSEnrichment")
)

# log10(nFrags) v TSS enrichment score plots generated upon arrow file creation


####################################
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
####################################

