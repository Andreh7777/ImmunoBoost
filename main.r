# Create a new directory and set the working directory
dir.create("~/meta_viz")
setwd("~/meta_viz") 

# Load the required libraries
library(tidyverse)
library(vegan)
library(tidyr)
library(microbiome)
library(knitr)
library(ggdist)
library(ggplot2)

# Load data files (taxonomy, metadata, OTU)
taxonomy = read_tsv("C:\\path\\dataframe.txt")
bug_df = read_delim("C:\\path\\OTU.txt", delim = " ")
metadata = read_tsv("C:\\path\\metadata.tsv")


# ---- Relative Abundance Calculation ----
# Convert the OTU data to a dataframe
bug_df <- as.data.frame(bug_df)
colnames(bug_df)[1] <- "OTU"  # Rename the first column to "OTU"
row.names(bug_df) <- bug_df[['OTU']]
bug_df <- subset(bug_df, select = -OTU)  # Remove the OTU column as row names are now set
bug.matrix = as.matrix(bug_df)  # Convert to matrix
OTU = otu_table(bug.matrix, taxa_are_rows = TRUE)  # Create OTU table for phyloseq object

# Prepare metadata
sample_names = c("ERR2162205_metaphlan","ERR2162202_metaphlan","ERR2162207_metaphlan","ERR2162208_metaphlan","ERR2162209_metaphlan","ERR2162211_metaphlan","ERR2162215_metaphlan","ERR2162214_metaphlan","ERR2162216_metaphlan","ERR2162217_metaphlan","ERR2162219_metaphlan","ERR2162218_metaphlan","ERR2162220_metaphlan","ERR2162222_metaphlan","ERR2162223_metaphlan","ERR2162224_metaphlan","ERR2162201_metaphlan","ERR2162200_metaphlan","ERR2162203_metaphlan","ERR2162204_metaphlan","ERR2162206_metaphlan","ERR2162210_metaphlan","ERR2162212_metaphlan","ERR2162213_metaphlan","ERR2162221_metaphlan")  # List of sample names
metadata = as.data.frame(metadata)
rownames(metadata) = sample_names
meta = as.matrix(metadata)
META = sample_data(metadata)  # Create sample data for phyloseq object

# Prepare taxonomy data
taxonomy <- as.data.frame(taxonomy)
taxonomy <- subset(taxonomy, select = -Sottospecie)  # Remove "Sottospecie" column if not needed
taxonomy <- unique(taxonomy)  # Ensure unique rows
row.names(taxonomy) <- taxonomy$OTU
taxonomy.matrix = as.matrix(taxonomy)
taxonomy <- subset(taxonomy, select = -OTU)
TAX = tax_table(taxonomy.matrix)  # Create taxonomy table for phyloseq object

# Create phyloseq object
physeq <- phyloseq(OTU, TAX, META)

# Normalize the OTU table by sample counts (relative abundance calculation)
otu_table(physeq) <- transform_sample_counts(otu_table(physeq), function(x) x / sum(x))

# Aggregate taxa at the Phylum level
phy.comp <- microbiome::aggregate_taxa(physeq, "Phylum")

# Create a color palette for the plot
palette = c("#641E16", "#C0392B", "#D98880", "#512E5F", "#8E44AD", "#D2B4DE",
            "#154360", "#2980B9", "#7FB3D5", "#D6EAF8", "#0E6251", "#1ABC9C", "#A3E4D7",
            "#138D75", "#27AE60", "#7D6608", "#F1C40F", "#F7DC6F", "#F9E79F", "#784212",
            "#E67E22", "#D35400", "#F5CBA7", "#7B7D7D", "#D0D3D4", "#FBFCFC", "#4D5656", "#1B2631", "#5D6D7E", "#17202A")

# Plot relative abundance at the Phylum level
g <- plot_bar(phy.comp, fill = "Phylum", facet_grid = "R_NR") +
  geom_bar(aes(fill = Phylum), stat = "identity") +
  scale_fill_manual(values = palette) +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~R_NR, scales = "free") +
  theme(panel.background = element_blank())

# Save the plot to a TIFF file
tiff(paste("Phylum_relative_abundance_status.tiff"), height = 20, width = 30, units = "cm", res = 300)
print(g)
dev.off()


# ---- Alpha Diversity Calculation ----
# Calculate alpha diversity using various indices (Shannon, Inverse Simpson, Chao1, etc.)
alpha_div <- microbiome::alpha(physeq, index = "all")

# Add alpha diversity metrics to the metadata
physeq.meta <- meta(physeq)
physeq.meta$Shannon <- alpha_div$diversity_shannon 
physeq.meta$InverseSimpson <- alpha_div$diversity_inverse_simpson
physeq.meta$Observed <- alpha_div$observed
physeq.meta$Chao1 <- alpha_div$chao1

# Display the first few rows of the metadata table with alpha diversity metrics
kable(head(physeq.meta))

# Create a boxplot for Chao1, shannon, simpson and observed index by R_NR grouping (change y = "index" )
p1 <- ggplot(physeq.meta, aes(x = R_NR, y = Chao1, fill = R_NR)) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5) +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) +
  labs(x = "R_NR", y = "Chao1") +
  theme_minimal()

# Save the alpha diversity plots (Chao1, shannon, simpson, observed) to a TIFF file. For example:
tiff("Alpha_diversity_Chao1.tiff", height = 20, width = 30, units = "cm", res = 300)
print(p1)
dev.off()


# ---- Beta Diversity Calculation ----

# Read in the bug data from a TSV file into a data frame
bug_df   = read_tsv("C:/path/bugs.tsv")

# Read in the metadata from a TSV file into a data frame
metadata = read_tsv("C:/path/metadata.tsv")

# Prepare the bug data matrix:
# - Convert the "OTU" column to row names
# - Convert the dataframe to a matrix
# - Transpose the matrix
bug_mat = bug_df |> 
  column_to_rownames("OTU") |> 
  as.matrix() |>
  t()

# Compute the distance matrix using the vegan package's vegdist function
dist_mat = vegdist(bug_mat)

# Perform Principal Coordinates Analysis (PCA) using cmdscale
# - k: Number of dimensions for the result (one less than the number of rows)
# - eig: Return eigenvalues
cmd_res = cmdscale(dist_mat, 
                   k = (nrow(bug_mat) - 1),
                   eig = TRUE)

# Display the structure of the PCA results
str(cmd_res)

# Create a tibble (data frame) with the first two principal coordinates
pcoa_df = tibble(PC1 = cmd_res$points[,1], 
                 PC2 = cmd_res$points[,2])

# Create a scatter plot of the principal coordinates
p = ggplot(pcoa_df, aes(x = PC1, y = PC2)) + 
  geom_point()
p

# Combine the PCA results with metadata side-by-side
pcoa_meta = bind_cols(pcoa_df, metadata)

# Display the combined data frame
pcoa_meta

# Create a scatter plot of the principal coordinates, colored by the "R_NR" variable
p_diag = ggplot(pcoa_meta,
                aes(x = PC1, y = PC2, color = R_NR)) + 
  geom_point() 

# Show the plot
p_diag

# Save the beta diversity plot to a TIFF file.
tiff("pcoa_plot.tiff", width = 800, height = 600, res = 300, compression = "lzw")
print(p_diag)
dev.off()