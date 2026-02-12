#### PACO analysis
library(ape)
library(dplyr)
library(tidyr)
library(bipartite)
library(vegan)
library(GUniFrac)
library(RPANDA)  
library(paco)
library(phangorn)

###############################################
# PACo analysis for host–parasite co-phylogeny
###############################################

#-------------------------------------------------
# 1. Read host (H) and parasite (P) trees with branch lengths
#-------------------------------------------------
H <- read.tree("/Users/syrenawhitner/Desktop/trees/4gvog_fungi.contree")
P <- read.tree("/Users/syrenawhitner/Desktop/trees/my_subtree_for_PACo.nwk")

# Confirm both have branch lengths (stop if not)
stopifnot(!is.null(H$edge.length), !is.null(P$edge.length))

#-------------------------------------------------
# 2. Root both trees consistently (here: midpoint root)
#-------------------------------------------------
H <- midpoint.root(H)
P <- midpoint.root(P)

# Clean up host tip labels if they have extra suffixes
H$tip.label <- gsub("_scaffolds\\.proteins$", "", H$tip.label)
P$tip.label <- gsub("_virus_proteins$", "", P$tip.label)

# Optional check
H$tip.label
P$tip.label

#-------------------------------------------------
# 3. Keep only shared taxa between both trees
#    (drop any mismatched tips so both trees contain the same labels)
#-------------------------------------------------
ids <- intersect(H$tip.label, P$tip.label)
H <- drop.tip(H, setdiff(H$tip.label, ids))
P <- drop.tip(P, setdiff(P$tip.label, ids))

#-------------------------------------------------
# 4. Tidy any non-positive edge lengths (very rare)
#    Some tree tools require all branch lengths > 0
#-------------------------------------------------
eps <- min(c(H$edge.length, P$edge.length)[c(H$edge.length, P$edge.length) > 0]) * 1e-6
H$edge.length[H$edge.length <= 0] <- eps
P$edge.length[P$edge.length <= 0] <- eps

#-------------------------------------------------
# 5. Compute pairwise phylogenetic distances (cophenetic matrices)
#-------------------------------------------------
DH <- cophenetic(H)
DP <- cophenetic(P)

#-------------------------------------------------
# 6. Build the host–parasite association matrix (L)
#    Here: 1-to-1 mapping (each host has exactly one parasite)
#-------------------------------------------------
L <- diag(1, nrow = length(ids), ncol = length(ids))
rownames(L) <- colnames(L) <- ids

#-------------------------------------------------
# 7. Prepare PACo input object (distance matrices + associations)
#    Apply Euclidean correction ("cailliez" or "lingoes") for ordination
#-------------------------------------------------
D <- prepare_paco_data(DH, DP, L)
D <- add_pcoord(D, correction = "cailliez")

#-------------------------------------------------
# 8. Run PACo with permutation test (global fit)
#-------------------------------------------------
set.seed(123)
res <- PACo(D, nperm = 9999, method = "r0")

# View global test statistic and permutation p-value
res$g   # contains m² (ss) and p

#-------------------------------------------------
# 9. Extract Procrustes results and compute residuals per link
#    Residual = Euclidean distance between host and parasite positions
#-------------------------------------------------
p <- res$proc
link_res <- sqrt(rowSums((p$Yrot - p$X)^2))
names(link_res) <- rownames(p$X)

#-------------------------------------------------
# 10. Quick diagnostics on residuals
#-------------------------------------------------
is.numeric(link_res)   # should be TRUE
length(link_res)       # number of associations
head(link_res)         # preview first few

# Top 5 most congruent (lowest residuals)
sort(link_res)[1:5]

# Top 5 least congruent (highest residuals)
sort(link_res, decreasing = TRUE)[1:5]

# Plot residual distributions
hist(link_res, breaks = 20,
     main = "PACo residuals",
     xlab = "Procrustes residuals (lower = more congruent)")

barplot(sort(link_res),
        las = 2, cex.names = 0.6,
        main = "PACo residuals per host–parasite pair",
        ylab = "Residual (lower = more congruent)")

#-------------------------------------------------
# 11. confirm naming alignment between residuals and PACo object
#-------------------------------------------------
head(rownames(p$X))
setequal(rownames(p$X), names(link_res))  # should be TRUE
all(rownames(p$X) == names(link_res))     # ideally TRUE

#-------------------------------------------------
# 12. Leave-one-out (LOO) sensitivity analysis
#     Tests whether removing any association changes the global p-value
#-------------------------------------------------
run_paco_subset <- function(keep_ids, H, P, L_full,
                            correction = "cailliez",
                            nperm = 9999, method = "r0", seed = 1) {
  Hk <- drop.tip(H, setdiff(H$tip.label, keep_ids))
  Pk <- drop.tip(P, setdiff(P$tip.label, keep_ids))
  Lk <- L_full[keep_ids, keep_ids, drop = FALSE]
  
  DHk <- cophenetic(Hk)
  DPk <- cophenetic(Pk)
  
  Dk <- prepare_paco_data(DHk, DPk, Lk)
  Dk <- add_pcoord(Dk, correction = correction)
  set.seed(seed)
  PACo(Dk, nperm = nperm, method = method)
}

# Build full L matrix for all current associations
ids_all <- rownames(res$proc$X)
L_full  <- diag(1, length(ids_all))
rownames(L_full) <- colnames(L_full) <- ids_all

# Compute LOO p-values
loo <- sapply(ids_all, function(bad) {
  keep <- setdiff(ids_all, bad)
  r <- run_paco_subset(keep, H, P, L_full)
  r$g$p
})
sort(loo)  # shows which removal improves global congruence most

#-------------------------------------------------
# 13. Export residuals table for external use or plotting
#-------------------------------------------------
residual_table <- data.frame(
  Association = names(link_res),
  Residual = round(link_res, 3)
) |>
  arrange(desc(Residual))

print(residual_table, row.names = FALSE)
# write.csv(residual_table, "PACo_link_residuals.csv", row.names = FALSE, quote = FALSE)

###############################################

library(viridis)

# Use the same host (H) and parasite (P) trees from your PACo run
# Make sure they contain only the shared tip labels
ids_match <- rownames(res$proc$X)   # PACo’s internal order (same as link_res names)

# Prune to those IDs
Hk <- drop.tip(H, setdiff(H$tip.label, ids_match))
Pk <- drop.tip(P, setdiff(P$tip.label, ids_match))

# Reorder both trees so that tips match PACo's internal order
Hk <- keep.tip(Hk, ids_match)
Pk <- keep.tip(Pk, ids_match)

# reorder tip labels explicitly to match PACo order
Pk$tip.label <- factor(Pk$tip.label, levels = ids_match)
Pk <- keep.tip(Pk, ids_match)
Pk$tip.label <- ids_match  # enforce consistent naming order

# confirm
all(Pk$tip.label == ids_match)
all(Hk$tip.label == ids_match)

co <- cophylo(Hk, Pk, rotate = TRUE)


left_order  <- co$trees[[1]]$tip.label   # host tip order (left)
right_order <- co$trees[[2]]$tip.label   # parasite tip order (right)

# Make sure residuals are named by IDs
stopifnot(!is.null(names(link_res)))

# Reorder residuals to the final order used in the plot
link_res_plot <- link_res[left_order]

# Build color vector from this aligned version
res_scaled <- (link_res_plot - min(link_res_plot)) /
  (max(link_res_plot) - min(link_res_plot))
pal <- colorRampPalette(c("purple", "#42a5f5", "#ffe082", "#e53935"))
edge_cols <- pal(100)[floor(res_scaled * 99) + 1]

plot(co,
     link.lwd = 2,
     link.col = edge_cols)

legend("topleft",
       fill = pal(5),
       legend = sprintf("%.1f",
                        seq(min(link_res_plot), max(link_res_plot), length.out = 5)),
       title = "Residual",
       cex = 0.8, bty = "n")

co <- cophylo(Hk, Pk, rotate = TRUE)

# Plot with *solid* black lines (lty = 1 ensures solid)
plot(
  co,
  link.lwd = 1.5,          # thicker lines for clarity
  link.col = "black",    # black color
  link.lty = 1,          # <--- solid lines
  fsize = 0.6
)

title("Host–Parasite Co-phylogeny (solid black links)")

residual_table <- data.frame(
  Association = names(link_res),
  Residual = round(link_res, 3)
) |>
  dplyr::arrange(desc(Residual))

breaks_custom <- c(2.3, 3.2, 4.2, 5.2, 7.5, Inf)
labels_custom <- c("Very congruent (2.4–3.2)",
                   "Congruent (3.2–4.2)",
                   "Intermediate (4.2–5.2)",
                   "Incongruent (5.2–7.5)",
                   "Strongly incongruent (>7.5)")

residual_table <- data.frame(
  Association = names(link_res),
  Residual = as.numeric(link_res)
) %>%
  mutate(
    Bin_custom = cut(Residual, breaks = breaks_custom, labels = labels_custom, right = TRUE, include.lowest = TRUE)
  )

resid_z <- scale(link_res)
resid_z


# extract z-scores as numeric vector
resid_vec <- as.numeric(resid_z)
names(resid_vec) <- names(link_res)  # ensure consistent naming

# reorder to host order in the left tree
resid_vec <- resid_vec[co$trees[[1]]$tip.label]

library(RColorBrewer)

pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))  # red = high z, blue = low z
z_min <- min(resid_vec, na.rm = TRUE)
z_max <- max(resid_vec, na.rm = TRUE)

# create colors
edge_cols <- pal(100)[as.numeric(cut(resid_vec, breaks = 100))]

plot(
  co,
  link.lwd = 2.5,
  link.col = edge_cols,
  link.lty = 1,
  fsize = 0.6
)

# add color legend
legend("topleft",
       fill = pal(5),
       legend = sprintf("%.1f", seq(z_min, z_max, length.out = 5)),
       title = "Z-score",
       cex = 0.8, bty = "n")

ncldv_cols <- read.csv("/Users/syrenawhitner/Downloads/colors - Sheet1.csv")

unique_pairs <- ncldv_cols %>%
  distinct(Order, color)
