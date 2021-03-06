library(tidyverse)
library(fs)
library(plotly)
library(recipes)
library(heatmaply)
library(RColorBrewer)
library(embed)


# read data ---------------------------------------------------------------

targets <- read_tsv("TAC-seq-gene-expression/data/targets/READY_72targets.tsv")

counts <- dir_ls("TAC-seq-gene-expression/data/counts/", glob = "*TAC-seq_counts_umi1.tsv") %>%
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

controls <- read_tsv("TAC-seq-gene-expression/data/controls/READY_HRT_controls_400k_reads.tsv") %>%
  select(sample, group, !!targets$target[targets$type == "biomarker"]) %>%
  mutate(group = factor(group, c("pre-receptive",
                                 "early-receptive",
                                 "receptive",
                                 "receptive HRT",
                                 "late-receptive",
                                 "post-receptive",
                                 "menstrual blood",
                                 "polyp")))


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
  na_if(0) %>%
  mutate(across(where(is.numeric), log)) %>%
  select(-file) %>%
  column_to_rownames("sample") %>%
  heatmaply(colors = rev(brewer.pal(n = 7, name = "RdYlBu")), scale = "column",
            show_dendrogram = c(TRUE, FALSE), hide_colorbar = TRUE)


# plot PCA ----------------------------------------------------------------

train_data %>%
  recipe() %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 2) %>%
  prep(strings_as_factors = FALSE) %>%
  bake(new_data = test_data) %>%
  ggplot(aes(PC1, PC2, color = group, shape = group, sample = sample)) +
  geom_point() +
  scale_shape_manual(values = 1:7) +
  labs(title = "PCA of biomarkers")
ggplotly()


# plot UMAP ---------------------------------------------------------------

train_data %>%
  recipe() %>%
  step_normalize(all_numeric()) %>%
  # step_umap(all_numeric(), seed = c(1, 1)) %>%
  step_string2factor(group) %>%
  step_umap(all_numeric(), outcome = vars(group), seed = c(1, 1)) %>%
  prep(strings_as_factors = FALSE) %>%
  bake(new_data = test_data) %>%
  ggplot(aes(umap_1, umap_2, color = group, shape = group, sample = sample)) +
  geom_point() +
  scale_shape_manual(values = 1:7) +
  labs(title = "UMAP of biomarkers")
ggplotly()
