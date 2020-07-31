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

patients <- list.files("TAC-seq-gene-expression/data/counts/", pattern = "TAC-seq_counts2.tsv", full.names = TRUE) %>%
  set_names(nm = basename(.)) %>%
  map_dfr(read_tsv, .id = "file") %>%
  filter(
    !str_detect(sample, "Undetermined"),  # remove undetermined samples
    locus != "unmatched"  # remove unmatched loci
  ) %>%
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
  stat_summary(aes(statistic = "geometric mean"), fun = geo_mean, geom = "crossbar") +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "Housekeeping genes", x = NULL, y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
ggplotly()


# normalize molecule counts -----------------------------------------------

# lst <- df %>%
#   select(sample, locus, molecule_count, type) %>%
#   spread(sample, molecule_count) %>%
#   split(.$type) %>%
#   map(select, -type) %>%
#   map(column_to_rownames, "locus") %>%
#   map(as.matrix)
#
# bm_raw <- lst$biomarker
# hk_raw <- lst$housekeeper
#
# norm_factor <- exp(apply(log(hk_raw), 2, mean))  # housekeepers geometric mean
# bm <- sweep(bm_raw, 2, norm_factor, "/")

patient_biomarkers <- patients %>%
  group_by(file, sample) %>%
  mutate(norm_factor = geo_mean(molecule_count[type == "housekeeper"])) %>%
  ungroup() %>%
  filter(type == "biomarker") %>%
  mutate(norm_molecule_count = molecule_count / norm_factor)

patient_biomarkers_clean <- patient_biomarkers %>%
  filter(is.finite(norm_molecule_count))

setdiff(patient_biomarkers$sample, patient_biomarkers_clean$sample)

patient_biomarkers_wide <- patient_biomarkers_clean %>%
  pivot_wider(id_cols = c(file, sample), names_from = locus, values_from = norm_molecule_count)

biomarkers <- bind_rows(controls, patient_biomarkers_wide)


# remove housekeeper outliers ---------------------------------------------

# hk_zeros <- which(norm_factor == 0)
# if (length(hk_zeros) > 0) {
#   outliers <- names(hk_zeros)
#   cat("removed samples:", outliers, sep = "\n")  # housekeepers geometric mean is 0
#   bm <- bm[, -hk_zeros]
# }


# plot PCA ----------------------------------------------------------------

# df <- bm %>%
#   t() %>%
#   as_tibble(rownames = "sample") %>%
#   bind_rows(controls)
#
# mat <- df %>%
#   select(-label) %>%
#   column_to_rownames("sample") %>%
#   as.matrix()
#
# pca <- mat %>%
#   prcomp(scale. = T)
#
# pca$x %>%
#   as_tibble(rownames = "sample") %>%
#   left_join(df) %>%
#   ggplot(aes(PC1, PC2, color = label, group = sample)) +
#   geom_point() +
#   labs(title = "Biomarkers", subtitle = "normalized molecule counts")
# ggplotly()

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

# mat[mat == 0] <- NA  # replace 0 with NA
# pheatmap(log(t(mat)), scale = "row", main = "Biomarkers", treeheight_row = 0)
biomarkers %>%
  column_to_rownames("sample") %>%
  select(where(is.numeric)) %>%
  t() %>%
  na_if(0) %>%
  log() %>%
  pheatmap(treeheight_row = 0, main = "Biomarkers")
