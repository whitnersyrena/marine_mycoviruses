### making the contig coverage plots 
# looking at the contig intron/gvog hit locations 
library(tidyverse)

##########################################################################################
# For a single genome  
##########################################################################################

#  Load data 
aug <- read_csv("/Users/syrenawhitner/Desktop/overlays/parsed/parsed_gff/UH_1001523_sa129663_combined_augustus_features.csv")
gvog <- read_csv("/Users/syrenawhitner/Desktop/overlays/parsed/parsed_domtbls_full/UH_1001523_sa129663_combined_gvog.complete_gvog_domains_genomic.csv")

gvog <- gvog %>%
  filter(as.numeric(evalue) <= 1e-3)

# Define GVOGs of interest 
target_gvogs <- c(
  "GVOGm0003","GVOGm0013","GVOGm0022","GVOGm0023",
  "GVOGm0054","GVOGm0172","GVOGm0461","GVOGm0760","GVOGm0890"
)

# Classify GVOG hits into target vs other 
gvog <- gvog %>%
  mutate(
    group = ifelse(feature %in% target_gvogs, "target", "other"),
    group = factor(group, levels = c("target", "other")) # ensure correct factor levels
  )

# Choose contigs to plot (optional) 
# Either use all contigs or only those containing ≥1 target GVOG
target_contigs <- gvog %>%
  filter(group == "target") %>%
  pull(contig) %>%
  unique()

aug_sub  <- aug  %>% filter(contig %in% target_contigs)
gvog_sub <- gvog %>% filter(contig %in% target_contigs)

# Color palette 
feature_colors <- c(
  "CDS"        = "#0072B2",   # blue
  "intron"     = "orange",   # dark gray
  "intergenic" = "#d9d9d9",   # light gray
  "target"     = "#E7298A",   # bright pink
  "other"      = "#B3A2C9"    # muted purple
)

# normalize strand position
gvog_sub <- gvog_sub %>%
  mutate(
    dom_start_norm = pmin(dom_genomic_start, dom_genomic_end),
    dom_end_norm   = pmax(dom_genomic_start, dom_genomic_end)
  )

# plot!

ggplot() +
  # contig backbone
  geom_segment(
    data = aug_sub %>%
      group_by(contig) %>%
      summarize(maxlen = max(end, na.rm = TRUE)),
    aes(x = 0, xend = maxlen, y = 0.01, yend = 0.01),
    color = "black", linewidth = 1, alpha = 0.3
  ) +
  # AUGUSTUS features (mapped color aesthetic)
  geom_segment(
    data = filter(aug_sub, feature %in% c("CDS", "intron", "intergenic")),
    aes(x = start, xend = end, y = 0.1, yend = 0.1, color = feature),
    linewidth = 6, alpha = 0.9
  ) +
  
  # GVOG hits (mapped color aesthetic)
  geom_segment(
    data = gvog_sub,
    aes(x = dom_genomic_start, xend = dom_genomic_end,
        y = 0.25, yend = 0.25, color = group),
    linewidth = 6
  ) +
  
  # Label only target GVOGs
  geom_text(
    data = filter(gvog_sub, group == "target"),
    aes(x = (dom_genomic_start + dom_genomic_end)/2, y = 0.35, label = feature),
    size = 2.5, color = feature_colors["target"], angle = 45, hjust = 0
  ) +
  
  # Shared color scale for both AUGUSTUS + GVOG layers
  scale_color_manual(
    values = feature_colors,
    name   = "Annotations",
    breaks = c("CDS", "intron", "intergenic", "target", "other"),
    labels = c("ORFs", "Introns", "Intergenic Regions",
               "NCLDV Hallmarks", "Other GVOGs")
  ) +
  
  facet_wrap(~ contig, scales = "free_x", ncol = 1) +
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing.y = unit(0.5, "lines")
  ) +
  labs(
    x = "Contig Length (bp)",
    title = "1002442 - AUGUSTUS (bottom) vs GVOG hits (top)"
  )



##########################################################################################
# looping through all the directories to build all the contig plots 
##########################################################################################

library(tidyverse)

# directories
aug_dir  <- "/Users/syrenawhitner/Desktop/overlays/parsed/parsed_gff"
gvog_dir <- "/Users/syrenawhitner/Desktop/overlays/parsed/parsed_domtbls_full"
out_dir  <- "/Users/syrenawhitner/Desktop/overlays/contig_plots_2"

dir.create(out_dir, showWarnings = FALSE)

# GVOG list
target_gvogs <- c(
  "GVOGm0003","GVOGm0013","GVOGm0022","GVOGm0023",
  "GVOGm0054","GVOGm0172","GVOGm0461","GVOGm0760","GVOGm0890"
)

# colors
feature_colors <- c(
  "CDS" = "#0072B2", 
  "intron" = "orange", 
  "intergenic" = "#d9d9d9",
  "target" = "#E7298A", 
  "other" = "#B3A2C9"
)

# just to combine everything, can normall skip # 
# list all gvog files
gvog_files <- list.files(gvog_dir, pattern = "_gvog.complete_gvog_domains_genomic.csv$", full.names = TRUE)

# read and combine all gvog files
all_gvogs <- gvog_files %>%
  map_dfr(function(f) {
    df <- suppressMessages(read_csv(f))
    
    # add filename as parent id
    df <- df %>%
      filter(as.numeric(evalue) <= 1e-5) %>%
      mutate(
        parent_file = basename(f) %>% gsub("_gvog.complete_gvog_domains_genomic.csv$", "", .),
        group = ifelse(feature %in% target_gvogs, "target", "other"),
        group = factor(group, levels = c("target", "other")),
        dom_start_norm = pmin(dom_genomic_start, dom_genomic_end),
        dom_end_norm   = pmax(dom_genomic_start, dom_genomic_end)
      )
    
    return(df)
  })

# bring in metadata 
gvog_annotat <-read_tsv("/Users/syrenawhitner/Downloads/gvog.complete.annot.tsv")

all_gvogs <- all_gvogs %>%
  left_join(
    gvog_annotat %>%
      select(GVOG, NCVOG_descs),
    by = c("feature" = "GVOG")
  ) %>%
  rename(desc = NCVOG_descs)

descriptions <- all_gvogs %>%
  select(feature, desc) %>%   # keep only GVOG ID and description
  distinct() %>%              # remove duplicates
  arrange(feature) 

write_csv(descriptions, "/Users/syrenawhitner/Desktop/overlays/descriptions2.csv")

target_contigs <- all_gvogs %>%
  filter(group == "target") %>%      # group == "target" marks hallmark GVOGs
  pull(contig) %>%
  unique()

# subset all gvogs that occur on any of those contigs
gvogs_same_contig <- all_gvogs %>%
  filter(contig %in% target_contigs)

gvogs_same_contig <- gvogs_same_contig %>%
  left_join(
    gvog_annotat %>%
      select(GVOG, NCVOG_descs),
    by = c("feature" = "GVOG")
  ) %>%
  rename(desc = NCVOG_descs)

descriptions2 <- gvogs_same_contig %>%
  select(feature, desc) %>%   # keep only GVOG ID and description
  distinct() %>%              # remove duplicates
  arrange(feature) 

write_csv(descriptions2, "/Users/syrenawhitner/Desktop/overlays/descriptions2.csv")

gvog_counts <- gvogs_same_contig %>%
  count(feature, name = "count") %>%
  arrange(desc(count)) %>% 
  slice_head(n = 15) 

# Quick look
head(gvog_counts)

# Plot, just looking at counts of differnt GVOGs
ggplot(gvog_counts, aes(x = reorder(feature, -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  labs(
    x = "GVOG feature",
    y = "Count (occurrences)",
    title = "Number of each GVOG detected across all genomes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5)
  )

###########

# list augustus files
aug_files <- list.files(aug_dir, pattern = "_augustus_features.csv$", full.names = TRUE)


for (aug_file in aug_files) {
  # derive matching gvog file
  base <- gsub("_augustus_features.csv$", "", basename(aug_file))
  gvog_file <- file.path(gvog_dir, paste0(base, "_gvog.complete_gvog_domains_genomic.csv"))
  if (!file.exists(gvog_file)) next
  
  # load
  aug  <- suppressMessages(read_csv(aug_file))
  gvog <- suppressMessages(read_csv(gvog_file))
  
  # filter by evalue
  gvog <- gvog %>% filter(as.numeric(evalue) <= 1e-5)
  
  # classify gvogs
  gvog <- gvog %>%
    mutate(
      group = ifelse(feature %in% target_gvogs, "target", "other"),
      group = factor(group, levels = c("target", "other")),
      dom_start_norm = pmin(dom_genomic_start, dom_genomic_end),
      dom_end_norm   = pmax(dom_genomic_start, dom_genomic_end)
    )
  
  # subset contigs that contain at least one target GVOG
  target_contigs <- gvog %>% filter(group == "target") %>% pull(contig) %>% unique()
  if (length(target_contigs) == 0) next
  
  aug_sub  <- aug  %>% filter(contig %in% target_contigs)
  gvog_sub <- gvog %>% filter(contig %in% target_contigs)
  
  # skip if empty
  if (nrow(aug_sub) == 0 | nrow(gvog_sub) == 0) next
  
  
  # plot
  p <- ggplot() +
    geom_segment(
      data = aug_sub %>%
        group_by(contig) %>%
        summarize(maxlen = max(end, na.rm = TRUE)),
      aes(x = 0, xend = maxlen, y = 0.01, yend = 0.01),
      color = "black", linewidth = 1, alpha = 0.3
    ) +
    geom_segment(
      data = filter(aug_sub, feature %in% c("CDS", "intron", "intergenic")),
      aes(x = start, xend = end, y = 0.1, yend = 0.1, color = feature),
      linewidth = 6, alpha = 0.9
    ) +
    geom_segment(
      data = filter(gvog_sub, group == "other"),
      aes(x = dom_start_norm, xend = dom_end_norm, y = 0.25, yend = 0.25, color = group),
      linewidth = 5, alpha = 0.8
    ) +
    geom_segment(
      data = filter(gvog_sub, group == "target"),
      aes(x = dom_start_norm, xend = dom_end_norm, y = 0.4, yend = 0.4, color = group),
      linewidth = 5
    ) +
    geom_text(
      data = filter(gvog_sub, group == "target"),
      aes(x = (dom_start_norm + dom_end_norm)/2, y = 0.47, label = feature),
      size = 2.5, color = feature_colors["target"], angle = 45, hjust = 0
    ) +
    scale_color_manual(
      values = feature_colors,
      name   = "Annotations",
      breaks = c("CDS","intron","intergenic","target","other"),
      labels = c("ORFs","Introns","Intergenic Regions","NCLDV Hallmarks","Other GVOGs")
    ) +
    facet_wrap(~ contig, scales = "free_x", ncol = 1) +
    coord_cartesian(ylim = c(0, 0.55)) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.spacing.y = unit(0.5, "lines")
    ) +
    labs(
      x = "Genomic position (bp)",
      title = paste0(base, " — AUGUSTUS (bottom), other GVOGs (middle), hallmark GVOGs (top)")
    )
  
  # dynamic height scaling
  n_contigs <- length(target_contigs)
  fig_height <- min(max(3, n_contigs * 2.5), 14)  # min 3", max 14"
  
  out_file <- file.path(out_dir, paste0(base, "_contig_plot.png"))
  ggsave(out_file, p, width = 9, height = fig_height, dpi = 300)
}