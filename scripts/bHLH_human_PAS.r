library(tidyverse)
library(patchwork)
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))


# ---------- paths ----------
# Prefer data/intermediate/table_input_PAS.csv, fallback to data/intermediate/table_input_withPAS.csv if missing.
input_table <- intermediate_csv_path("table_input_PAS.csv")
if (!file.exists(input_table) && file.exists(p("data", "intermediate", "table_input_withPAS.csv"))) {
  input_table <- p("data", "intermediate", "table_input_withPAS.csv")
}
output_dir <- p("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_plot_svg <- file.path(output_dir, "bHLH_human_position_PAS.svg")
output_plot_pdf <- file.path(output_dir, "bHLH_human_position_PAS.pdf")

# ---------- load data ----------
df <- read_csv(input_table, show_col_types = FALSE)

# check minimal expected columns
required_cols <- c("HGNC.symbol", "Sequence", "interpro", "interpro_start", "interpro_end", "length_aa", "number.of.isoforms", "class")
missing <- setdiff(required_cols, colnames(df))
if(length(missing)>0) stop("Missing required columns in input CSV: ", paste(missing, collapse = ", "))

# ---------- normalize types ----------
df <- df %>%
  mutate(
    interpro_start = as.numeric(interpro_start),
    interpro_end   = as.numeric(interpro_end),
    length_aa      = as.numeric(length_aa),
    number.of.isoforms = as.numeric(number.of.isoforms),
    class = as.character(class)
  )

# ensure class values are canonical (A B C D E NC)
df$class[is.na(df$class) | df$class == ""] <- "NC"

# ---------- Prepare plotting data (single source df) ----------
# gene order: keep df order but reversed for vertical layout
gene_levels <- df %>% pull(HGNC.symbol) %>% unique() %>% rev()
df <- df %>% mutate(HGNC.symbol = factor(HGNC.symbol, levels = gene_levels))

# genes_df: one row per gene (for p2/p3)
genes_df <- df %>%
  group_by(HGNC.symbol, class) %>%
  summarise(
    prot_length = if(all(is.na(length_aa))) NA_real_ else max(length_aa, na.rm = TRUE),
    number_of_isoforms = if(all(is.na(number.of.isoforms))) NA_real_ else max(number.of.isoforms, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(HGNC.symbol = factor(HGNC.symbol, levels = gene_levels))

# compute numeric y positions to align layers that need numeric coords
genes_df <- genes_df %>% mutate(y_num = as.numeric(HGNC.symbol))

# domain_df (PAS = IPR000014) absolute coords (aa)
domain_abs <- df %>%
  filter(interpro == "IPR000014") %>%
  filter(!is.na(interpro_start) & !is.na(interpro_end)) %>%
  group_by(HGNC.symbol) %>%
  summarise(domain_start = min(as.numeric(interpro_start), na.rm = TRUE),
            domain_end   = max(as.numeric(interpro_end), na.rm = TRUE),
            .groups = "drop") %>%
  mutate(HGNC.symbol = factor(HGNC.symbol, levels = gene_levels)) %>%
  left_join(genes_df %>% select(HGNC.symbol, prot_length, y_num), by = "HGNC.symbol") %>%
  filter(!is.na(y_num)) %>%
  mutate(domain_start = pmax(domain_start, 1),
         domain_end   = pmin(domain_end, ifelse(is.na(prot_length), domain_end, prot_length)),
         valid = domain_end >= domain_start) %>%
  filter(valid)

# domain_df relative coords (for p1 alignment 0-1)
domain_rel <- domain_abs %>%
  mutate(xmin_rel = domain_start / prot_length,
         xmax_rel = domain_end / prot_length)

# highlight_df for bHLH (IPR011598) relative positions 0-1
highlight_bhlh <- df %>%
  filter(interpro == "IPR011598") %>%
  filter(!is.na(Sequence) | !is.na(length_aa)) %>%
  mutate(seq_len = ifelse(!is.na(length_aa), length_aa, nchar(Sequence))) %>%
  filter(!is.na(interpro_start) & !is.na(interpro_end) & seq_len > 0) %>%
  rowwise() %>%
  mutate(xmin = interpro_start / seq_len, xmax = interpro_end / seq_len) %>%
  ungroup() %>%
  select(HGNC.symbol, xmin, xmax, class) %>%
  mutate(class = ifelse(is.na(class) | class == "", "NC", class)) %>%
  filter(!is.na(HGNC.symbol))

# highlight_df_pas for PAS relative positions (for p1 overlay)
highlight_pas_rel <- domain_rel %>%
  select(HGNC.symbol, xmin_rel, xmax_rel, y_num) %>%
  left_join(df %>% select(HGNC.symbol, class) %>% distinct(), by = "HGNC.symbol") %>%
  mutate(class = ifelse(is.na(class) | class == "", "NC", class))

# ---------- plotting aesthetics ----------
class_colors_vec <- c(
  "A" = "#7AC36A", "B" = "#F7CB45", "C" = "#9063CD",
  "D" = "#D83F3D", "E" = "#4C86C6", "NC"= "grey"
)
domain_ipr000014_color <- "#FF8C00"

# combined fill for class legend only (domain drawn as separate geom with constant fill)
combined_fill <- c(class_colors_vec, "PAS domain" = domain_ipr000014_color)


# ---------- p1: relative positions (0-1) with both bHLH and PAS (relative) ----------
p1 <- ggplot() +
  # bHLH rectangles filled by class
  geom_rect(data = highlight_bhlh,
            aes(xmin = xmin, xmax = xmax,
                ymin = as.numeric(HGNC.symbol) - 0.4, ymax = as.numeric(HGNC.symbol) + 0.4,
                fill = class),
            color = NA) +
  # PAS relative rectangles (same y positions), fixed color (overlaid)
  geom_rect(data = highlight_pas_rel,
            aes(xmin = xmin_rel, xmax = xmax_rel,
                ymin = y_num - 0.4, ymax = y_num + 0.4, fill = "PAS domain"),
            color = NA, alpha = 0.9) +

  geom_vline(xintercept = 0.5, color = 'gray30', linetype = 'dashed') +
  scale_y_continuous(breaks = seq_along(gene_levels), labels = gene_levels, expand = c(0,0)) +
  scale_fill_manual(values = combined_fill, name = "Class") + 
  labs(x = "Relative position within sequence (0-1)", y = "Gene symbol") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 7, face = "bold"),
        axis.ticks.y = element_blank())

# ---------- p2: horizontal bars showing prot_length (aa) + absolute PAS overlay ----------
# We'll draw bars with y = HGNC.symbol (factor) so bars are horizontal and line up with p1 labels
p2 <- ggplot() +
  geom_col(data = genes_df, aes(x = prot_length, y = HGNC.symbol, fill = class), width = 0.6) +
  # overlay PAS in absolute aa coordinates; need numeric vertical positions, compute center +/- offset
  
  geom_rect(data = domain_abs,
            inherit.aes = FALSE,
            aes(xmin = domain_start, xmax = domain_end,
              ymin = y_num - 0.3, ymax = y_num + 0.3, fill = "PAS domain"),
            color = NA, alpha = 0.9) +

  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = combined_fill, guide = "none") +
  labs(x = "Length of the longest isoform (aa)", y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 8))

# ---------- p3: number of isoforms per gene (horizontal bars) ----------
p3 <- ggplot(genes_df, aes(x = number_of_isoforms, y = HGNC.symbol, fill = class)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = combined_fill, guide = "none") +
  labs(x = "Number of isoforms per gene", y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 8))

# ---------- Combine plots: use guides collected and shared legend for class ------
# We want the legend only for class (not for PAS) so we collect from p1 (which has fill mapped to class)
combined <- (p1 + p2 + p3) + plot_layout(ncol = 3, widths = c(2.5, 1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 8))

print(combined)

# ---------- Save ----------
ggsave(output_plot_svg, plot = combined, width = 40, height = 24, units = "cm")
ggsave(output_plot_pdf, plot = combined, width = 40, height = 24, units = "cm")
message("Plot saved.")
