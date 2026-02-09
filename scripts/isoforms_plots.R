# === Load Libraries ===
library(ggplot2)
library(tidyverse)
library(patchwork)
library(tidyr)
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



# Output directory
output_dir <- p("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data files
df <- read.csv(intermediate_csv_path("bHLH_StartEnd_withISO.csv")) # nolint
df_classes <- read.csv(p("data", "raw", "LS_classes.csv")) # nolint
df_classes <- normalize_classes(df_classes)
df_dym <- read.csv(p("data", "raw", "Lambert_bHLH.csv")) # nolint
df_dym <- df_dym %>% select(HGNC.symbol, Binding.mode) %>% distinct()

df_classes$Ledent2002.Simionato2007 <- ifelse(df_classes$Ledent2002.Simionato2007 == "", "NC", df_classes$Ledent2002.Simionato2007)

split(df, paste(df$HGNC.symbol, df$ensembl_transcript_id)) %>%
lapply(function(tdf){
  data.frame(
    HGNC.symbol = tdf$HGNC.symbol,
    ensembl_transcript_id = tdf$ensembl_transcript_id,
    bin = c('pre', plyr::mapvalues(x = tdf$HGNC.symbol, 
                                   from = df_classes$HGNC.symbol,
                                   to = df_classes$Ledent2002.Simionato2007, warn_missing = F), 'post'),
    size = c(tdf$interpro_start,
             tdf$interpro_end - tdf$interpro_start,
             tdf$protein_length - tdf$interpro_end)
  )
  }
) %>% do.call('rbind', .) -> bin_sizes


# Define class color mapping
class_colors <- c(
  "A" = "#7AC36A",   # Green
  "B" = "#F7CB45",   # Yellow
  "C" = "#9063CD",   # Purple
  "D" = "#D83F3D",   # Red
  "E" = "#4C86C6",   # Blue
  "NC" = "#cb4289"      # Unclassified in pink
)

bin_sizes$bin <- factor(bin_sizes$bin, levels = c('pre', names(class_colors), 'post'))
bin_sizes$class = plyr::mapvalues(x = bin_sizes$HGNC.symbol, 
                                  from = df_classes$HGNC.symbol,
                                  to = df_classes$Ledent2002.Simionato2007, warn_missing = F)

table(bin_sizes$bin)


split(bin_sizes, bin_sizes$class) %>%
lapply(function(pdata){

  pdata.size <- pdata %>% group_by(HGNC.symbol, ensembl_transcript_id) %>% summarise(len = sum(size))
ggplot(pdata, aes(y = ensembl_transcript_id, x = size, fill = bin)) +
  
  geom_bar(color = 'black', linewidth = 0.25, stat = 'identity', position = position_fill(reverse = TRUE), show.legend = F) +
  
  geom_point(data = pdata.size, 
             inherit.aes = T, 
             mapping = aes(fill = NULL, color = len, size = len, x = 1.05)) +
  
  geom_text(data = pdata.size, show.legend = F,
       size = 2,
            inherit.aes = T,
            mapping = aes(fill = NULL, label = len, x = 1.1)) + 

  scale_size_continuous(range = c(1,3)) + 
  scale_color_viridis_c(option='turbo') +
  scale_fill_manual(values = c('pre' = 'lightgrey', 'post' = 'lightgrey', class_colors), guide = guide_none()) + 
  facet_grid(rows = vars(HGNC.symbol), space = 'free', scales = 'free',
             switch = 'y') +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1, 1.1),
                     labels = c(0,'25%','50%','75%','100%', 'Transcript\nlength (aa)'),
                     #labels = scales::label_percent(),
                     expand = expansion(mult = c(0,0.05))) + 
  scale_y_discrete(expand = expansion(mult = c(0,0))) + 
  theme(strip.placement = 'outside', 
  panel.spacing.y = unit(0.3, "lines"),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, color = 'black'),
    strip.text.y.left = element_text(angle = 0, hjust = 1, vjust=0.5)) + 
  labs(fill = 'Gene class', x = 'Relative position of bHLH domain', y = 'Transcript ID', 
       color = 'Transcript\nlength (aa)', size = 'Transcript\nlength (aa)')
}
) -> plots_transcripts_human


ggsave(file.path(output_dir, "bHLH_human_transcript_classA_A3.svg"), plot =plots_transcripts_human$A, width = 12, height = 12, units = "in")
plots_transcripts_human$A
# Save plot for Class B
ggsave(file.path(output_dir, "bHLH_human_transcript_classB_A3.svg"), plot = plots_transcripts_human$B, width = 12, height = 12, units = "in")
plots_transcripts_human$B
# Save plot for Class C
ggsave(file.path(output_dir, "bHLH_human_transcript_classC_A3.svg"), plot = plots_transcripts_human$C, width = 12, height = 12, units = "in")
plots_transcripts_human$C
# Save plot for Class D
ggsave(file.path(output_dir, "bHLH_human_transcript_classD_A3.svg"), plot = plots_transcripts_human$D, width = 12, height = 12, units = "in")
plots_transcripts_human$D
# Save plot for Class E
ggsave(file.path(output_dir, "bHLH_human_transcript_classE_A3.svg"), plot = plots_transcripts_human$E, width = 12, height = 12, units = "in")
plots_transcripts_human$E
# Save plot for Class NC (Not Classified)
ggsave(file.path(output_dir, "bHLH_human_transcript_classNC_A3.svg"), plot = plots_transcripts_human$NC, width = 12, height = 12, units = "in")
plots_transcripts_human$NC
