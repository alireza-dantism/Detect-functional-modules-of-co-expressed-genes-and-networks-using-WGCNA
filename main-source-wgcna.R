# Detect functional modules of co-expressed genes and networks using WGCNA (R package)
# Created on Tue May 12 2022

library(WGCNA)      # weighted gene correlation networks for analysis
library(tidyverse)  # for using ggplot2
library(magrittr)   # provides the %>% operator

data <- readr::read_delim("GSE128078.txt", delim = "\t") # Load and clean data => convert txt data to frame data 
data <- data[1:43497,1:5] # only use 4 samples

# change title of tracking_id to Id and Sample_G1-* to sample*
names(data)[1] = "Id"
names(data)[2] = "sample1"
names(data)[3] = "sample2"
names(data)[4] = "sample3"
names(data)[5] = "sample4"



col_sel = names(data)[-1] # Get all but first column name

# Group all of date base on 4 samples
mdata <- data %>%
  tidyr::pivot_longer(
    ., col = all_of(col_sel)
  ) %>% mutate(group = gsub("-.*","", name) %>% gsub("[.].*","", .))


# Plot groups (Sample Groups vs RNA Seq Counts) to identify outliers
(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = samples, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)          # Rotate sample text
    ) +
    labs(x = "Sample Groups", y = "RNA Seq Counts") +
    facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")
)

# ----------- Normalize data ----------------------

library(DESeq2) # Use this library to normalize the counts before sending to WGCNA

# Prepare DESeq input, which is expecting a matrix of integers.
de_input = as.matrix(data[,-1]) 
row.names(de_input) = data$Id # base on Samples Id


# Group samples and their types
meta_df <- data.frame( Sample = names(data[-1])) %>% mutate( Type = gsub("-.*","", Sample) %>% gsub("[.].*","", .) )

dds <- DESeqDataSetFromMatrix(round(de_input), meta_df, design = ~1) # converting counts to integer mode
dds <- DESeq(dds) # estimating size factor, dispersions and ...

vsd <- varianceStabilizingTransformation(dds) # calculate variance


library(genefilter) # Use thi library for filtering genes from high-throughput experiments


wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)


q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn , ]


# Create a plot which normalized
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "samples",
    y = "normalized expression"
  )





# ------- Transpose the data and prepare the dataset for WGCNA --------

input_mat = t(expr_normalized) # transpose data

allowWGCNAThreads() # allow multi-threading (optional)

powers = c(c(1:10), seq(from = 12, to=20, by=2)) # Choose a set of soft-thresholding powers

sft = pickSoftThreshold(input_mat, powerVector = powers, verbose = 1) # Call the network topology analysis function


# Creating Scale independece and Mean connectivity plots
par(mfrow = c(1,2))

cex1 = 0.9

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)

text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "blue" )

abline(h = 0.90, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)

text(sft$fitIndices[, 1],sft$fitIndices[, 5],labels = powers, cex = cex1, col = "red")



# Pick a soft threshold power near the curve of the plot

picked_power = 18 # we could pick 16, 18 or 20.
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat, # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power, 
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 1266,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == 
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 2)


cor <- temp_cor     # Return cor function to original namespace


mergedColors = labels2colors(netwk$colors) # Convert labels to colors for plotting

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )




# -------------- figure out which modules are associated with each trait/treatment group -----------

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add samples names
MEs0$samples = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-samples) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=samples, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


