library(GGally)
library(dplyr)
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


species_order <- c(
  "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae", "Neurospora_crassa",
  "Caenorhabditis_elegans", "Anopheles_gambiae", "Drosophila_melanogaster",
  "Tribolium_castaneum", "Helobdella_robusta",
  "Danio_rerio", "Oryzias_latipes", "Lepisosteus_oculatus",
  "Xenopus_tropicalis", "Gallus_gallus", "Anolis_carolinensis",
  "Monodelphis_domestica", "Bos_taurus", "Canis_lupus_familiaris",
  "Mus_musculus", "Rattus_norvegicus",
  "Macaca_mulatta", "Gorilla_gorilla", "Pan_troglodytes", "Homo_sapiens"
) %>% tolower()

classes= read.csv(p("data", "raw", "LS_classes.csv"))
classes <- normalize_classes(classes)
data= read.csv(p("data", "intermediate", "orthologs", "annotated_bHLH_merged_data_with_gene_names.csv"))
data <- data %>% left_join(classes %>% select(HGNC.symbol, Ledent2002.Simionato2007), by = "HGNC.symbol")
data <- data %>% subset(homology_type == 'ortholog_one2one') %>%
mutate(midpoint = (Rel_start_T + Rel_end_T) / 2)

## Remove genes with duplicated bHLH domains in the same species
dup.data <- data %>% group_by(HGNC.symbol, target_species) %>% summarise(n=n()) %>% arrange(desc(n))
dup.genes <- dup.data %>% subset(n > 1) %>% select(HGNC.symbol) %>% unlist() %>% unique()
data <- data %>% subset(!HGNC.symbol %in% dup.genes)
# data %>% group_by(HGNC.symbol, target_species) %>% summarise(n=n()) %>% arrange(desc(n))


data %>% select(HGNC.symbol, Rel_start_Q, Rel_end_Q) %>% unique %>%
 mutate(homo_sapiens = (Rel_start_Q + Rel_end_Q) / 2, 
        Rel_start_Q = NULL, Rel_end_Q = NULL) -> data.human

data.wide <- reshape2::dcast('HGNC.symbol ~ target_species', value.var = 'midpoint', data = data)
data.merged <- merge(data.human, data.wide, by = 'HGNC.symbol', all = TRUE)

sel_species <- species_order[species_order%in%colnames(data.merged)]
data.merged <- data.merged[,c('HGNC.symbol', sel_species)]


# p <- ggparcoord(data.merged,
#     columns = 2:6
#    # columns = 2:ncol(data.wide)
#     ) 
# ggsave(p, filename = "parallel_coordinates_plot.png", width = 10, height = 6)
# print(p)


library(plotly)

data.merged[is.na(data.merged)] <- -.5
colnames(data.merged)[-1] %>% lapply(function(x) {
    list(label = x, values = data.merged[[x]], range = c(-.5,1))
}) -> dimensions


classes$Ledent2002.Simionato2007 <- ifelse(classes$Ledent2002.Simionato2007 == "", "NC", classes$Ledent2002.Simionato2007)
data.merged$class = plyr::mapvalues(x = data.merged$HGNC.symbol, 
                                  from = classes$HGNC.symbol,
                                  to = classes$Ledent2002.Simionato2007, warn_missing = F)


table(data.merged$class)

class_colors <- list(
  c(0.0, "#7AC36A"),   # Green
  c(0.2, "#F7CB45"),   # Yellow
  c(0.4, "#9063CD"),   # Purple
  c(0.6, "#D83F3D"),   # Red
  c(0.8, "#4C86C6"),   # Blue
  c(1, "#cb4289")      # Unclassified in pink
)

data.merged$class_num <- factor(data.merged$class, levels = names(class_colors)) %>% as.numeric
data.merged$class_fr <- (data.merged$class_num - min(data.merged$class_num)) / (max(data.merged$class_num) - min(data.merged$class_num))
# 
# data.merged[,c('class_num', 'class')] %>% 
#  unique() %>% arrange(class_num) %>%
#  apply(1, function(x) {
#     c(x[1], class_colors[x[2]])
# }) %>% lapply(function(x){x})-> class_colors_list
# 

# colorscale <- lapply(names(class_colors), function(x){c(x, class_colors[[x]])})



fig <- data.merged %>% plot_ly(type = 'parcoords',
          line = list(color = ~class_fr,
                      colorscale = class_colors),
          dimensions = dimensions)
fig





### try 2

library(GGally)
library(dplyr)
species_order <- c(
  "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae", "Neurospora_crassa",
  "Caenorhabditis_elegans", "Anopheles_gambiae", "Drosophila_melanogaster",
  "Tribolium_castaneum", "Helobdella_robusta",
  "Danio_rerio", "Oryzias_latipes", "Lepisosteus_oculatus",
  "Xenopus_tropicalis", "Gallus_gallus", "Anolis_carolinensis",
  "Monodelphis_domestica", "Bos_taurus", "Canis_lupus_familiaris",
  "Mus_musculus", "Rattus_norvegicus",
  "Macaca_mulatta", "Gorilla_gorilla", "Pan_troglodytes", "Homo_sapiens"
) %>% tolower()

classes= read.csv(p("data", "raw", "LS_classes.csv"))
classes <- normalize_classes(classes)
data= read.csv(p("data", "intermediate", "orthologs", "annotated_bHLH_merged_data_with_gene_names.csv"))
data <- data %>% left_join(classes %>% select(HGNC.symbol, Ledent2002.Simionato2007), by = "HGNC.symbol")
data <- data %>% mutate(midpoint = (Rel_start_T + Rel_end_T) / 2)


non.o2o_genes <- data %>% group_by(HGNC.symbol, target_species) %>% summarise(n=n()) %>% subset(n > 1) %>% select(HGNC.symbol) %>% unlist() %>% unique()
data.o2o <- data %>% subset(!HGNC.symbol %in% non.o2o_genes) %>% select(-source_seq,-target_seq)
data.o2o <- data.o2o %>% subset(!is.na(midpoint))
subset(data.o2o, is.na(midpoint)) %>% select(Rel_start_T, Rel_end_T, midpoint, Start_T, Stop_T)

head(data.o2o)


classes$Ledent2002.Simionato2007 <- ifelse(classes$Ledent2002.Simionato2007 == "", "NC", classes$Ledent2002.Simionato2007)
data.o2o$class = plyr::mapvalues(x = data.o2o$HGNC.symbol, 
                                  from = classes$HGNC.symbol,
                                  to = classes$Ledent2002.Simionato2007, warn_missing = F)


subset(data.o2o, target_species == 'drosophila_melanogaster') %>% subset(class == 'C')

ggplot(data.o2o, aes(x = target_species, y = midpoint, group = HGNC.symbol, color = class)) + 
  geom_line(alpha = 0.6) + 
  geom_point(alpha = 0.6, size = 1) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) +
  labs(title = "One-to-One Orthologs Midpoint Distribution",
       x = "Target Species",
       y = "Midpoint")