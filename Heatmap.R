########################### 
# Heatmap Figure Supplementary Figure 1
# Syrena Whitner
# December 2025 
########################### 

# Load libraries
library(dplyr)
library(stringr)
library(pheatmap)

common1 <- read_csv("/Users/syrenawhitner/Desktop/overlays/parsed/parsed_domtbls_common/GVOG_summary_common.csv")

#--------------------------------------------
# Clean 'origin' column
# Remove everything from "_combined" onward
#--------------------------------------------
common1$origin <- str_replace(common1$origin, "_combined.*", "")
df<-common1
df_subset<-common1

#--------------------------------------------
# Subset to selected origins
#--------------------------------------------
origin_set <- c(
  "1002442_sa182533", "1002610_sa182534", "1003183_sa150016", "1003189_sa149994",
  "1003356_sa149980", "UH_1000136_sa123987", "UH_1000339_sa123914", "UH_1000590_sa144727",
  "UH_1001061_sa129894", "UH_1001187_sa123975", "UH_1001220_sa145563", "UH_1001271_sa144725",
  "UH_1001376_sa129866", "UH_1001502_sa132352", "UH_1001523_sa129663", "UH_1001535_sa132328",
  "UH_1001556_sa132306", "UH_1002031_sa145543"
)

df_subset <- df %>%
  filter(origin %in% origin_set)

#--------------------------------------------
# Prepare matrix for heatmap
#--------------------------------------------
# Set 'origin' as rownames and drop that column
mat <- df_subset %>%
  column_to_rownames("origin") %>%
  as.matrix()

# Ensure numeric
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- df_subset$origin

#--------------------------------------------
# Plot heatmap (white → teal)
#--------------------------------------------
heat_common <- pheatmap(
  mat,
  color = colorRampPalette(c("white", "#0092a0"))(100),
  cluster_rows = TRUE,          
  cluster_cols = FALSE,          
  display_numbers = TRUE,        
  number_format = "%.0f",
  number_color = "black",
  fontsize_number = 8,
  border_color = "lightgrey",
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45,
  main = "GVOG Counts Across Genomes"
)

heat_common

# Save to a PDF, e.g. 8 inches wide × 6 inches tall
pdf("/Users/syrenawhitner/Desktop/viruses/Sept2025/int_figs/GVOG_heatmap_common.pdf", width = 10, height = 6)

# Print the heatmap object to the PDF
print(heat_common)

# Close the PDF device
dev.off()
