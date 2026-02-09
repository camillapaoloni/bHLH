# isoform_plot_merged_final.R - robust version without plyr conflicts
# Load packages (avoid plyr)
if("package:plyr" %in% search()) detach("package:plyr", unload = TRUE, character.only = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(viridisLite) # viridis palette backend (scale_color_viridis_c uses viridisLite)
# (do not load 'plyr' to avoid conflicts)
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))
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



# === File paths ===
path_df <- intermediate_csv_path("bHLH_StartEnd_withISO.csv")
path_df_classes <- p("data", "raw", "LS_classes.csv")
path_df_dym <- p("data", "raw", "Lambert_bHLH.csv")

# === Read data ===
df <- read.csv(path_df, stringsAsFactors = FALSE)
df_classes <- read.csv(path_df_classes, stringsAsFactors = FALSE)
df_dym <- read.csv(path_df_dym, stringsAsFactors = FALSE) %>%
  select(HGNC.symbol, Binding.mode) %>% distinct()

# Normalize missing classes to "NC"
df_classes$Ledent2002.Simionato2007 <- ifelse(df_classes$Ledent2002.Simionato2007 == "" | is.na(df_classes$Ledent2002.Simionato2007),
                                              "NC",
                                              df_classes$Ledent2002.Simionato2007)

# Add class directly to df for easier manipulation
df2 <- df %>%
  left_join(df_classes %>% select(HGNC.symbol, Ledent2002.Simionato2007),
            by = "HGNC.symbol") %>%
  rename(class = Ledent2002.Simionato2007)

# === Build bin_sizes (vectorized, 3 rows per transcript) ===
# size_mid = interpro_end - interpro_start (domain), pre = interpro_start, post = protein_length - interpro_end
bin_sizes <- df2 %>%
  mutate(
    size_pre  = ifelse(is.na(interpro_start), 0, interpro_start),
    size_mid  = ifelse(is.na(interpro_end) | is.na(interpro_start), 0, interpro_end - interpro_start),
    size_post = ifelse(is.na(protein_length) | is.na(interpro_end), 0, protein_length - interpro_end),
    orig_ensembl = ensembl_transcript_id
  ) %>%
  select(HGNC.symbol, ensembl_transcript_id, class, orig_ensembl, size_pre, size_mid, size_post) %>%
  pivot_longer(cols = starts_with("size_"),
               names_to = "bin_short",
               values_to = "size") %>%
  mutate(
    bin = case_when(
      bin_short == "size_pre"  ~ "pre",
      bin_short == "size_mid"  ~ class,  # here the bin name is the class (A,B,...,NC)
      bin_short == "size_post" ~ "post",
      TRUE ~ as.character(bin_short)
    )
  ) %>%
  select(HGNC.symbol, ensembl_transcript_id, orig_ensembl, bin, size, class)

# If you want to remove rows with size <= 0 (optional)
bin_sizes <- bin_sizes %>% mutate(size = ifelse(is.na(size), 0, size))

# Ensure bin levels are ordered as desired:
# class letters present
class_levels_present <- sort(unique(na.omit(bin_sizes$class)))
desired_class_order <- intersect(c("A","B","C","D","E","NC"), class_levels_present)

# bin levels: pre, (class letters...), post
bin_levels <- c("pre", desired_class_order, "post")
# Force levels (if a class is missing, fct_expand handles it)
bin_sizes$bin <- factor(bin_sizes$bin, levels = bin_levels)

# === Compute transcript lengths (absolute) ===
transcript_lengths <- bin_sizes %>%
  filter(!is.na(orig_ensembl) & orig_ensembl != "") %>%
  group_by(orig_ensembl) %>%
  summarise(len = sum(size, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame(stringsAsFactors = FALSE)

# === Define transcript ordering and labels ===
transcript_order <- bin_sizes %>%
  distinct(class, HGNC.symbol, orig_ensembl) %>%
  mutate(class = factor(class, levels = desired_class_order)) %>%
  arrange(class, HGNC.symbol, orig_ensembl) %>%
  pull(orig_ensembl)

label_map <- bin_sizes %>%
  distinct(orig_ensembl, HGNC.symbol) %>%
  arrange(match(orig_ensembl, transcript_order)) %>%
  mutate(ensembl_label = paste0(HGNC.symbol, " | ", orig_ensembl))

# Merge label and lengths
bs2 <- bin_sizes %>%
  mutate(ensembl_label = factor(orig_ensembl, levels = transcript_order, labels = label_map$ensembl_label)) %>%
  left_join(transcript_lengths, by = "orig_ensembl")

# === Define colors ===
class_colors <- c("A" = "#7AC36A", "B" = "#F7CB45", "C" = "#9063CD",
                  "D" = "#D83F3D", "E" = "#4C86C6", "NC"= "#cb4289")

# Build named vector for all bins (pre/post + classes)
fill_bins_all <- c(pre = "lightgrey", post = "lightgrey")
# Add colors for present classes
for(cl in desired_class_order){
  if(cl %in% names(class_colors)) fill_bins_all[cl] <- class_colors[cl]
}
# Ensure name order matches bs2$bin levels
fill_bins_all <- fill_bins_all[as.character(bin_levels)]
# Remove any NA entries (classes not present)
fill_bins_all <- fill_bins_all[!is.na(names(fill_bins_all))]

# === Separator lines between classes ===
group_counts <- bs2 %>%
  distinct(class, orig_ensembl, ensembl_label) %>%
  mutate(class = factor(class, levels = desired_class_order)) %>%
  arrange(class) %>%
  count(class) %>%
  mutate(cum = cumsum(n))

sep_positions <- group_counts$cum + 0.5
if(length(sep_positions) > 0) sep_positions <- sep_positions[-length(sep_positions)]

# === Build unified plot ===
p_single <- ggplot(bs2, aes(y = ensembl_label, x = size, fill = bin)) +
  geom_col(position = position_fill(reverse = TRUE), color = "black", linewidth = 0.2, show.legend = TRUE) +
  scale_fill_manual(values = fill_bins_all, na.value = "grey50") +
  # Points for transcript length (placed slightly to the right of fill bars)
  geom_point(data = distinct(bs2, ensembl_label, len),
             aes(x = 1.05, y = ensembl_label, color = len, size = len),
             inherit.aes = FALSE, show.legend = TRUE) +
  geom_text(data = distinct(bs2, ensembl_label, len),
            aes(x = 1.07, y = ensembl_label, label = len),
            inherit.aes = FALSE, hjust = 0, size = 2.8) +
  scale_color_viridis_c(option = "turbo", name = "Transcript\nlength (aa)") +
  scale_size_continuous(range = c(1.5, 4), guide = "none") +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = c("0","25%","50%","75%","100%"),
                     expand = expansion(mult = c(0,0.02))) +
  coord_cartesian(xlim = c(0,1.15)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = "Relative position of bHLH domain", fill = "Bin / Class")

# Add separators if present
if(length(sep_positions) > 0) p_single <- p_single + geom_hline(yintercept = sep_positions, color = "black", linewidth = 0.3)

# Force discrete legend (order)
p_single <- p_single + guides(fill = guide_legend(order = 1))

# === Save plot ===
n_transcripts <- length(unique(bs2$ensembl_label))
per_iso_inch <- 0.12
total_height <- per_iso_inch * n_transcripts + 1
total_width <- 12
output_dir <- p("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(output_dir, "bHLH_singleplot_isoforms.svg")

ggsave(out_file, p_single, width = total_width, height = total_height, units = "in", dpi = 300)
message("Saved: ", out_file, " (w=", total_width, "in, h=", total_height, "in)")
