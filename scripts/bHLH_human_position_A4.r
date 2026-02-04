# Tests to generate the file in horizontal A4
library(ggplot2)
library(tidyverse)
library(patchwork)
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



# Output directory
output_dir <- p("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
df <- read.csv(p("data", "intermediate", "table_input.csv"))
df_classes <- read.csv(p("data", "raw", "LS_classes.csv"))
df_classes <- normalize_classes(df_classes)

# Merge and order
df_all <- df %>%
  left_join(df_classes, by = "HGNC.symbol")

df$HGNC.symbol <- factor(df$HGNC.symbol, levels = rev(unique(df$HGNC.symbol)))
df_all$HGNC.symbol <- factor(df_all$HGNC.symbol, levels = levels(df$HGNC.symbol))

# Define class colors
class_colors <- scale_fill_manual(values = c(
  "A" = "#7AC36A",
  "B" = "#F7CB45",
  "C" = "#9063CD",
  "D" = "#D83F3D",
  "E" = "#4C86C6",
  "NC" = "grey"
))

# Build highlight_df
highlight_df <- df %>%
  mutate(
    length = nchar(Sequence),
    relative_position = map(length, ~ seq(0, 1, length.out = .x))
  ) %>%
  unnest(relative_position) %>%
  mutate(
    index = round(relative_position * (length - 1)) + 1,
    highlight = index >= interpro_start & index <= interpro_end
  ) %>%
  filter(highlight) %>%
  group_by(HGNC.symbol) %>%
  summarise(
    xmin = min(relative_position),
    xmax = max(relative_position),
    .groups = "drop"
  ) %>%
  left_join(df_classes, by = "HGNC.symbol") %>%
  mutate(class = ifelse(Ledent2002.Simionato2007 == "", "NC", Ledent2002.Simionato2007)) %>%
  filter(!is.na(class))

# Plot 1: Highlight rectangles (domain)
p1 <- ggplot(highlight_df, aes(
  xmin = xmin, xmax = xmax,
  ymin = as.numeric(factor(HGNC.symbol, levels = levels(df$HGNC.symbol))) - 0.4,
  ymax = as.numeric(factor(HGNC.symbol, levels = levels(df$HGNC.symbol))) + 0.4,
  fill = class
)) +
  geom_rect() +
  geom_vline(xintercept = 0.5, color = 'gray30', linetype = 'dashed') +
  scale_y_continuous(
    breaks = seq_along(levels(df$HGNC.symbol)),
    labels = levels(df$HGNC.symbol),
    expand = c(0, 0)
  ) +
  labs(x = "Relative Position of the bHLH domain", y = "Gene Symbol", fill = "Atchley(1997) Ledent(2002) classification") +
  class_colors +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 5, face = "bold"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 15, face = "bold"),
    legend.position = "none"
  )

# Plot 2: Sequence length bar
p2 <- df_all %>%
  distinct(HGNC.symbol, length, Ledent2002.Simionato2007) %>%
  mutate(class = ifelse(Ledent2002.Simionato2007 == "", "NC", Ledent2002.Simionato2007)) %>%
  ggplot(aes(y = HGNC.symbol, x = length, fill = class)) +
  geom_bar(stat = "identity") +
  labs(x = "Length of the Longest Isoform", y = NULL) +
  class_colors +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0),
    axis.text.x = element_text(size = 12, ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 30, face = "bold"),
    legend.position = "none"
  )

# Plot 3: Isoform count bar
p3 <- df_all %>%
  distinct(HGNC.symbol, number.of.isoforms, Ledent2002.Simionato2007) %>%
  mutate(class = ifelse(Ledent2002.Simionato2007 == "", "NC", Ledent2002.Simionato2007)) %>%
  ggplot(aes(y = HGNC.symbol, x = number.of.isoforms, fill = class)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of Isoforms per Gene", y = NULL, fill = "Class") +
  class_colors +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.3),

    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 50, face = "bold"),
    legend.position = "right",  # or "right"
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

# Combine plots
final_plot <- p1 + p2 + p3 + plot_layout(ncol = 3, widths = c(2.5, 1, 1))

print(final_plot)

# For horizontal A4 (29.7 x 21 cm => ~11.7 x 8.3 in)
ggsave(file.path(output_dir, "bHLH_human_position_A4.svg"), plot = final_plot, width = 11.7, height = 8.3, units = "in", dpi = 300)
