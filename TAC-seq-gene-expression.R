library(tidyverse)
library(plotly)
library(pheatmap)


# read data ---------------------------------------------------------------

targets <- read_tsv("TAC-seq-gene-expression/data/targets/READY_targets.tsv")
controls <- read_tsv("TAC-seq-gene-expression/data/controls/READY_controls.tsv")

df <- read_tsv("TAC-seq-gene-expression/data/counts/TAC-seq_counts.tsv") %>%
  filter(
    !str_detect(sample, "Undetermined"),  # remove undetermined samples
    locus != "unmatched"  # remove unmatched loci
  ) %>%
  left_join(targets, by = c("locus" = "id"))

# plot biomarkers ---------------------------------------------------------

df %>%
  filter(type == "biomarker") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus)) +
  scale_y_log10() +
  labs(title = "Biomarkers", subtitle = "raw molecule counts", y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))

# plot spike-ins ----------------------------------------------------------

df %>%
  filter(type == "spike_in") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus)) +
  labs(title = "ERCC spike-ins", subtitle = "raw molecule counts", y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))

# plot housekeeper --------------------------------------------------------

df %>%
  filter(type == "housekeeper") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus)) +
  labs(title = "Housekeeping genes", subtitle = "raw molecule counts", y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))

# normalize molecule counts -----------------------------------------------

lst <- df %>%
  select(sample, locus, molecule_count, type) %>%
  spread(sample, molecule_count) %>%
  split(.$type) %>%
  map(select, -type) %>%
  map(column_to_rownames, "locus") %>%
  map(as.matrix)

bm_raw <- lst$biomarker
hk_raw <- lst$housekeeper

norm_factor <- exp(apply(log(hk_raw), 2, mean))  # geometric mean of housekeepers
bm_norm <- sweep(bm_raw, 2, norm_factor, "/")

# remove housekeeper outliers ---------------------------------------------

hk_zeros <- which(norm_factor == 0)
if (length(hk_zeros) > 0) {
  bm_norm <- bm_norm[, -hk_zeros]
  cat("removed sample(s) where geometric mean of housekeeping genes is zero:", names(hk_zeros), sep = "\n")
}

# plot PCA ----------------------------------------------------------------

test_ctrl <- bm_norm %>%
  t() %>%
  as_tibble(rownames = "sample") %>%
  bind_rows(controls)

pca <- test_ctrl %>%
  select(-label) %>%
  column_to_rownames("sample") %>%
  as.matrix() %>%
  prcomp(scale. = T)

pca$x %>%
  as_tibble(rownames = "sample") %>%
  left_join(test_ctrl) %>%
  ggplot(aes(PC1, PC2, color = label, group = sample)) +
  geom_point() +
  labs(title = "Biomarkers", subtitle = "normalized molecule counts")
ggplotly()

# plot heatmap ------------------------------------------------------------

bm_norm[bm_norm == 0] <- NA  # replace 0 with NA
pheatmap(log(bm_norm), scale = "row", main = "Biomarkers", treeheight_row = 0)
