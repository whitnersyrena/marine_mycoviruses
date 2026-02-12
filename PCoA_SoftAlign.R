########################### 
# PCoA Figure 3 
# converts soft align results into distance matrix and generates PCoA 
# NOTE: this is only for one of the GVOGs, each was run indidually (same script, just swapped out the input file)
# all individual GVOGs were then combined into one single plot in inkscape
# Syrena Whitner
# December 2025 
########################### 

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Input files - note this is just for one GVOG, each was done indiviaully and combined later 
# similarity table (SimScore-based)
similarity_file <- "/Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0023/GVOGm0023_output.tsv"

# metadata CSV with at least two columns: id, source
metadata_file   <- "/Users/syrenawhitner/Downloads/Pairwise_IDs-GVOGm0023(1).csv"

# 2. Read input data 
df <- read_tsv(similarity_file, show_col_types = FALSE)
meta <- read_csv(metadata_file, show_col_types = FALSE)

meta$id <- meta$ID

# Prepare similarity pairs
df <- df %>%
  select(QueryID, HitID, SimScore) %>%
  filter(!is.na(SimScore), QueryID != HitID)

ids <- sort(unique(c(df$QueryID, df$HitID)))

S <- df %>%
  group_by(QueryID, HitID) %>%
  summarize(SimScore = max(SimScore), .groups = "drop") %>%
  pivot_wider(names_from = HitID, values_from = SimScore, values_fn = max) %>%
  rename(id = QueryID) %>%
  right_join(tibble(id = ids), by = "id") %>%
  arrange(id)

S_mat <- as.matrix(S[,-1])
rownames(S_mat) <- S$id
missing_cols <- setdiff(ids, colnames(S_mat))
if (length(missing_cols) > 0) {
  add <- matrix(NA_real_, nrow = nrow(S_mat), ncol = length(missing_cols),
                dimnames = list(rownames(S_mat), missing_cols))
  S_mat <- cbind(S_mat, add)
}
S_mat <- S_mat[ids, ids, drop = FALSE]

S_sym <- pmax(S_mat, t(S_mat), na.rm = TRUE)
diag(S_sym) <- 1
S_sym[is.na(S_sym)] <- 0

# Distance matrix + PCoA 
D <- 1 - S_sym
Ddist <- as.dist(D)
mds <- cmdscale(Ddist, eig = TRUE, k = 2, add = TRUE)

coords <- data.frame(
  id   = rownames(mds$points),
  Dim1 = mds$points[,1],
  Dim2 = mds$points[,2],
  stringsAsFactors = FALSE
)

# Merge metadata 
coords <- coords %>%
  left_join(meta, by = "id")

# Plot with ggplot

ggplot(coords, aes(x = Dim1, y = Dim2, color = source)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 13) +
  labs(title = "PCoA (GVOGm0013 / SFII)",
       x = "Axis 1", y = "Axis 2", color = "Source") +
  theme(panel.grid = element_line(color = "gray90"),
        legend.position = "right")

my_colors <- c(
  "fungi"  = "#90c24e",
  "NCLDV"  = "#5BA3FF",
  "syrena" = "#a83e9a",
  "Mycodnoviridae" = "#5BA3FF"
)

pcoa <- ggplot(coords, aes(x = Dim1, y = Dim2, color = source)) +
  geom_point(size = 3, alpha = 0.9, shape = 16) +  # shape 16 = solid circle, no outline
  theme_minimal(base_size = 13) +
  labs(
    title = "PCoA - GVOGm0023 / RNAPL",
    x = "Axis 1",
    y = "Axis 2",
    color = "Source"
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right"
  ) +
  scale_color_manual(values = my_colors, na.value = "grey80")


# Save to a PDF, e.g. 8 inches wide × 6 inches tall
pdf("/Users/syrenawhitner/Desktop/viruses/Sept2025/int_figs/pcoa_0023.pdf", width = 7, height = 5.2)

# Print the heatmap object to the PDF
print(pcoa)

# Close the PDF device
dev.off()

