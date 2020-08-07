library(tidyverse)
library(fs)
library(plotly)
library(recipes)
library(embed)
library(pheatmap)


# read data ---------------------------------------------------------------

targets <- read_tsv("TAC-seq-gene-expression/data/targets/READY_61targets.tsv")

counts <- dir_ls("TAC-seq-gene-expression/data/counts/", glob = "*TAC-seq_counts.tsv") %>%
  map_dfr(read_tsv, .id = "file") %>%
  mutate(file = path_file(file)) %>%
  filter(!str_detect(sample, "Undetermined"),  # remove undetermined samples
         locus != "unmatched") %>%  # remove unmatched loci
  right_join(targets, by = c("locus" = "target")) %>%
  complete(nesting(file, sample), nesting(locus, type)) %>%
  group_by(sample) %>%
  filter(n() == nrow(targets),  # remove duplicated samples
         !any(is.na(molecule_count))) %>%  # remove samples with missing targets
  mutate(hk_geo_mean = exp(mean(log(molecule_count[type == "housekeeper"]))),  # geometric mean of housekeeping genes
         norm_molecule_count = molecule_count / hk_geo_mean) %>%
  ungroup()

controls <- read_tsv("TAC-seq-gene-expression/data/controls/READY65_control_set.tsv") %>%
  select(sample, label, !!targets$target[targets$type == "biomarker"])


# plot counts -------------------------------------------------------------

counts %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus)) +
  scale_y_log10() +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "Counts", x = NULL, y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
ggplotly()


# plot housekeeper --------------------------------------------------------

counts %>%
  filter(type == "housekeeper") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_point(aes(color = locus)) +
  geom_errorbar(aes(y = hk_geo_mean, ymin = hk_geo_mean, ymax = hk_geo_mean),
                size = 1) +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "Housekeeping genes", x = NULL, y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
ggplotly()


# biomarker counts --------------------------------------------------------

norm_biomarkers <- counts %>%
  filter(type == "biomarker") %>%
  pivot_wider(id_cols = c(file, sample), names_from = locus,
              values_from = norm_molecule_count) %>%
  filter(across(where(is.numeric), is.finite))  # remove samples with geometric mean of zero


# train and test data -----------------------------------------------------

train_data <- controls
test_data <- bind_rows(norm_biomarkers, controls)


# plot heatmap ------------------------------------------------------------

test_data %>%
  column_to_rownames("sample") %>%
  select(where(is.numeric)) %>%
  t() %>%
  na_if(0) %>%
  log() %>%
  pheatmap(treeheight_row = 0, main = "Biomarkers")


# plot PCA ----------------------------------------------------------------

train_data %>%
  recipe() %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 2) %>%
  prep(strings_as_factors = FALSE) %>%
  bake(new_data = test_data) %>%
  ggplot(aes(PC1, PC2, color = label, sample = sample)) +
  geom_point() +
  labs(title = "PCA of biomarkers", color = NULL)
ggplotly()


# plot UMAP ---------------------------------------------------------------

train_data %>%
  recipe() %>%
  step_normalize(all_numeric()) %>%
  # step_umap(all_numeric(), seed = c(1, 1)) %>%
  step_string2factor(label) %>%
  step_umap(all_numeric(), outcome = vars(label), seed = c(1, 1)) %>%
  prep(strings_as_factors = FALSE) %>%
  bake(new_data = test_data) %>%
  ggplot(aes(umap_1, umap_2, color = label, sample = sample)) +
  geom_point() +
  labs(title = "UMAP of biomarkers", color = NULL)
ggplotly()
