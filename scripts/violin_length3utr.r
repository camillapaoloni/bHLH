library(ggplot2)
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



input_path <- p("data", "intermediate", "rbp", "Violin_plot_3UTR.csv")
output_dir <- p("outputs", "rbp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdata <- read.csv(input_path)
if ("Ledent2002+Simionato2007" %in% names(pdata)) {
  names(pdata)[names(pdata) == "Ledent2002+Simionato2007"] <- "Ledent2002.Simionato2007"
}

# Boxplot of the median
p_median <- ggplot(pdata, aes(x = Ledent2002.Simionato2007, y = UTR_length_median, fill = Ledent2002.Simionato2007)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.4) +
  labs(title = "Median 3'UTR length") +
  theme_classic() +
  theme(legend.position = "none")

# Boxplot of the mean
p_mean <- ggplot(pdata, aes(x = Ledent2002.Simionato2007, y = UTR_length_mean, fill = Ledent2002.Simionato2007)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.4) +
  labs(title = "Mean 3'UTR length") +
  theme_classic() +
  theme(legend.position = "none")

# Boxplot of the maximum
p_max <- ggplot(pdata, aes(x = Ledent2002.Simionato2007, y = UTR_length_max, fill = Ledent2002.Simionato2007)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.4) +
  labs(title = "Max 3'UTR length") +
  theme_classic() +
  theme(legend.position = "none")
my_colors <- c(
  "A" = "#7AC36A",   # green
  "B" = "#F7CB45",   # yellow
  "C" = "#9063CD",   # purple
  "D" = "#D83F3D",   # red
  "E" = "#4C86C6",   # blue
  "NC" = "grey"      # gray
)

p_median <- p_median + scale_fill_manual(values = my_colors)
p_mean   <- p_mean   + scale_fill_manual(values = my_colors)
p_max    <- p_max    + scale_fill_manual(values = my_colors)
# Combine plots
final_plot <- p_median + p_mean + p_max + plot_layout(ncol = 3)
print(final_plot)
# Save as SVG
ggsave(file.path(output_dir, "UTR_length_boxplots.svg"), plot = final_plot, width = 30, height = 10)
