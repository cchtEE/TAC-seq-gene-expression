library(tidyverse)
library(plotly)
library(recipes)
library(embed)
library(pheatmap)


# read data ---------------------------------------------------------------

targets <- read_tsv("TAC-seq-gene-expression/data/targets/READY76_target_set.tsv")

controls <- "TAC-seq-gene-expression/data/controls/READY65_control_set.tsv" %>%
  set_names(nm = basename(.)) %>%
  map_dfr(read_tsv, .id = "file")

patients <- list.files("TAC-seq-gene-expression/data/counts/",
                       pattern = "TAC-seq_counts.tsv", full.names = TRUE) %>%
  set_names(nm = basename(.)) %>%
  map_dfr(read_tsv, .id = "file") %>%
  filter(!str_detect(sample, "Undetermined"),  # remove undetermined samples
         locus != "unmatched") %>%  # remove unmatched loci
  left_join(targets, by = c("locus" = "target"))


# plot biomarkers ---------------------------------------------------------

patients %>%
  filter(type == "biomarker") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus), show.legend = FALSE) +
  scale_y_log10() +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "Biomarkers", x = NULL, y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
ggplotly()


# plot spike-ins ----------------------------------------------------------

patients %>%
  filter(type == "spike_in") %>%
  ggplot(aes(sample, molecule_count, color = locus)) +
  geom_point() +
  geom_line(aes(group = locus)) +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "ERCC spike-ins", x = NULL, y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
ggplotly()


# plot housekeeper --------------------------------------------------------

geo_mean <- function(x) {
  # Compute the sample geometric mean.
  exp(mean(log(x)))
}

patients %>%
  filter(type == "housekeeper") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_point(aes(color = locus)) +
  stat_summary(aes(statistic = "geometric mean"), fun = geo_mean,
               geom = "crossbar") +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "Housekeeping genes", x = NULL, y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
ggplotly()


# normalize molecule counts -----------------------------------------------

biomarkers <- patients %>%
  group_by(file, sample) %>%
  mutate(norm_factor = geo_mean(molecule_count[type == "housekeeper"])) %>%
  ungroup() %>%
  filter(type == "biomarker") %>%
  mutate(norm_molecule_count = molecule_count / norm_factor) %>%
  filter(is.finite(norm_molecule_count)) %>%
  pivot_wider(id_cols = c(file, sample), names_from = locus,
              values_from = norm_molecule_count) %>%
  bind_rows(controls)


# plot PCA ----------------------------------------------------------------

controls %>%
  recipe() %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 2) %>%
  prep(strings_as_factors = FALSE) %>%
  bake(new_data = biomarkers) %>%
  ggplot(aes(PC1, PC2, color = label, sample = sample, file = file)) +
  geom_point() +
  labs(title = "PCA of biomarkers", color = NULL)
ggplotly()


# plot UMAP ---------------------------------------------------------------

controls %>%
  recipe() %>%
  step_normalize(all_numeric()) %>%
  step_umap(all_numeric(), seed = c(1, 1)) %>%
  # step_string2factor(label) %>%
  # step_umap(all_numeric(), outcome = vars(label), seed = c(1, 1)) %>%
  prep(strings_as_factors = FALSE) %>%
  bake(new_data = biomarkers) %>%
  ggplot(aes(umap_1, umap_2, color = label, sample = sample, file = file)) +
  geom_point() +
  labs(title = "UMAP of biomarkers", color = NULL)
ggplotly()


# plot heatmap ------------------------------------------------------------

biomarkers %>%
  column_to_rownames("sample") %>%
  select(where(is.numeric)) %>%
  t() %>%
  na_if(0) %>%
  log() %>%
  pheatmap(treeheight_row = 0, main = "Biomarkers")
