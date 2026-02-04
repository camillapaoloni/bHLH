library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
# === File paths: edit if needed ===
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)
normalize_classes <- function(df) {
  if ("Ledent2002+Simionato2007" %in% names(df)) {
    df <- dplyr::rename(df, Ledent2002.Simionato2007 = `Ledent2002+Simionato2007`)
  }
  if (!"Ledent2002.Simionato2007" %in% names(df)) {
    stop("Class column not found in LS_classes")
  }
  df$Ledent2002.Simionato2007 <- ifelse(is.na(df$Ledent2002.Simionato2007) | df$Ledent2002.Simionato2007 == "",
                                          "NC", df$Ledent2002.Simionato2007)
  df
}


path_df <- p("data", "intermediate", "bHLH_StartEnd_withISO.csv")
path_df_classes <- p("data", "raw", "LS_classes.csv")
path_df_dym <- p("data", "raw", "Lambert_bHLH.csv")

# === Read ===
df <- read.csv(path_df, stringsAsFactors = FALSE)
df_classes <- read.csv(path_df_classes, stringsAsFactors = FALSE)
df_dym <- read.csv(path_df_dym, stringsAsFactors = FALSE) %>% select(HGNC.symbol, Binding.mode) %>% distinct()

# Normalize missing class to NC (original behavior)
df_classes$Ledent2002.Simionato2007 <- ifelse(df_classes$Ledent2002.Simionato2007 == "", "NC", df_classes$Ledent2002.Simionato2007)

# === Build bin_sizes (same logic as the working script) ===
bin_sizes <- split(df, paste(df$HGNC.symbol, df$ensembl_transcript_id)) %>%
  lapply(function(tdf){
    data.frame(
      HGNC.symbol = rep(tdf$HGNC.symbol, 3),
      ensembl_transcript_id = rep(tdf$ensembl_transcript_id, 3),
      bin = c('pre',
              plyr::mapvalues(x = tdf$HGNC.symbol,
                              from = df_classes$HGNC.symbol,
                              to = df_classes$Ledent2002.Simionato2007, warn_missing = FALSE),
              'post'),
      size = c(tdf$interpro_start,
               tdf$interpro_end - tdf$interpro_start,
               tdf$protein_length - tdf$interpro_end),
      stringsAsFactors = FALSE
    )
  }) %>% do.call('rbind', .)

# Add class column explicitly
bin_sizes$class <- plyr::mapvalues(x = bin_sizes$HGNC.symbol,
                                   from = df_classes$HGNC.symbol,
                                   to = df_classes$Ledent2002.Simionato2007, warn_missing = FALSE)


# 0) Ensure 'bin' factor has logical order (pre .. domain .. post)
# (This avoids 'post' being drawn before other bins in the stack)
bin_sizes$bin <- factor(bin_sizes$bin, levels = c('pre','A','B','C','D','E','NC','post'))

# 1) Keep the original Ensembl ID for joins/ordering
bs <- bin_sizes %>%
  mutate(orig_ensembl = ensembl_transcript_id,
         class = plyr::mapvalues(HGNC.symbol,
                                 from = df_classes$HGNC.symbol,
                                 to   = df_classes$Ledent2002.Simionato2007,
                                 warn_missing = FALSE))

# 2) Compute absolute lengths (B) before changing labels
transcript_lengths <- bs %>%
  group_by(orig_ensembl) %>%
  summarise(len = sum(size, na.rm = TRUE), .groups = "drop")

# 3) Define desired class order (edit if you want a different order)
desired_class_order <- intersect(c("A","B","C","D","E","NC"), unique(bs$class))

# 4) Build transcript order: group by class -> gene -> transcript
transcript_order <- bs %>%
  distinct(class, HGNC.symbol, orig_ensembl) %>%
  mutate(class = factor(class, levels = desired_class_order)) %>%
  arrange(class, HGNC.symbol, orig_ensembl) %>%
  pull(orig_ensembl)

# 5) Build combined labels "HGNC | ENST..." (shown on y-axis)
label_map <- bs %>%
  distinct(orig_ensembl, HGNC.symbol) %>%
  arrange(match(orig_ensembl, transcript_order)) %>%
  mutate(ensembl_label = paste0(HGNC.symbol, " | ", orig_ensembl))

# 6) Build a plotting df: add 'ensembl_label' and lengths
bs2 <- bs %>%
  mutate(ensembl_label = factor(orig_ensembl, levels = transcript_order, labels = label_map$ensembl_label)) %>%
  left_join(transcript_lengths, by = c("orig_ensembl"))

# 7) Define palette for 'bin' levels: pre/post = gray, A..NC = class colors
class_colors <- c("A" = "#7AC36A", "B" = "#F7CB45", "C" = "#9063CD",
                  "D" = "#D83F3D", "E" = "#4C86C6", "NC"= "#cb4289")
# Bin palette (includes pre/post)
fill_bins_all <- c(pre = "lightgrey",
                   post = "lightgrey",
                   class_colors[intersect(names(class_colors), unique(bs2$class))])
# Ensure names exist
fill_bins_all <- fill_bins_all[!is.na(names(fill_bins_all))]

# 8) Compute positions for separator lines (discrete y coordinates)
group_counts <- bs2 %>% distinct(class, orig_ensembl, ensembl_label) %>%
  mutate(class = factor(class, levels = desired_class_order)) %>%
  arrange(class) %>%
  count(class, name = "n") %>%
  mutate(cum = cumsum(n))
sep_positions <- group_counts$cum + 0.5
if (length(sep_positions) > 0) sep_positions <- sep_positions[-length(sep_positions)]

# 9) Build the unified plot
p_single <- ggplot(bs2, aes(y = ensembl_label, x = size)) +
  # 9a) Stacked bar normalized per transcript (0..1) — fill = bin (A, B, pre, post...)
  geom_bar(aes(fill = bin), stat = "identity", position = position_fill(reverse = TRUE),
           color = "black", linewidth = 0.2, show.legend = FALSE) +
  # 9b) Bin color scale: A..NC get class colors, pre/post gray
  scale_fill_manual(values = fill_bins_all,
                    # Show only class letters in the legend (exclude pre/post)
                    breaks = intersect(names(class_colors), unique(as.character(bs2$bin))),
                    labels = intersect(names(class_colors), unique(as.character(bs2$bin)))) +
  # 9c) Point indicating absolute transcript length (B)
  geom_point(data = distinct(bs2, ensembl_label, len),
             aes(x = 1.05, y = ensembl_label, color = len, size = len),
             inherit.aes = FALSE, show.legend = TRUE) +
  scale_color_viridis_c(option = "turbo", name = "Transcript\nlength (aa)") +
  scale_size_continuous(range = c(1.5, 4), guide = "none") +
  # 9d) Axes, limits, styling
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), labels = c("0","25%","50%","75%","100%"),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 1.15)) +
  # 9e) Visual separators between classes
  { if(length(sep_positions)>0) geom_hline(yintercept = sep_positions, color = "black", linewidth = 0.3) } +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = "Relative position of bHLH domain", fill = "Class (domain color)")

# 10) Force discrete legend for 'fill' (show only classes A..NC)
p_single <- p_single + guides(fill = guide_legend(order = 1))

# 11) Save with calibrated height
n_transcripts <- length(unique(bs2$ensembl_label))
per_iso_inch <- 0.12
total_height <- per_iso_inch * n_transcripts + 1
total_width <- 12
out_file <- "bHLH_singleplot_isoforms.svg"
output_dir <- p("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(output_dir, "bHLH_singleplot_isoforms.svg")
ggsave(out_file, p_single, width = total_width, height = total_height, units = "in", dpi = 300)
message("Saved: ", out_file, " (w=", total_width, "in, h=", total_height, "in)")
