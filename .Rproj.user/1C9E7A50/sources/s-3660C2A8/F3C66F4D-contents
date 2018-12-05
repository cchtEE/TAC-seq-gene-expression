library(tidyverse)
library(plotly)
library(pheatmap)


# read data ---------------------------------------------------------------

targets <- read_tsv("TAC-seq-gene-expression/data/targets/READY_targets.tsv")
controls <- read_tsv("TAC-seq-gene-expression/data/controls/READY_controls.tsv")

df <- list.files("TAC-seq-gene-expression/data/counts/", pattern = "TAC-seq_counts.tsv", full.names = T) %>%
# df <- list.files("../../../Desktop/", "READY", full.names = T) %>%
  set_names(nm = basename(.)) %>%
  map_dfr(read_tsv, .id = "file") %>%

  filter(
    !str_detect(sample, "Undetermined"),  # remove undetermined samples
    locus != "unmatched"  # remove unmatched loci
  ) %>%
  left_join(targets, by = c("locus" = "target"))

# plot biomarkers ---------------------------------------------------------

df %>%
  filter(type == "biomarker") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus)) +
  scale_y_log10() +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "Biomarkers", subtitle = "raw molecule counts", y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))

# plot spike-ins ----------------------------------------------------------

df %>%
  filter(type == "spike_in") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus)) +
  facet_wrap(vars(file), scales = "free") +
  labs(title = "ERCC spike-ins", subtitle = "raw molecule counts", y = "molecule count") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90))

# plot housekeeper --------------------------------------------------------

df %>%
  filter(type == "housekeeper") %>%
  ggplot(aes(sample, molecule_count)) +
  geom_boxplot() +
  geom_point(aes(color = locus)) +
  facet_wrap(vars(file), scales = "free") +
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

norm_factor <- exp(apply(log(hk_raw), 2, mean))  # housekeepers geometric mean
bm <- sweep(bm_raw, 2, norm_factor, "/")

# remove housekeeper outliers ---------------------------------------------

hk_zeros <- which(norm_factor == 0)
if (length(hk_zeros) > 0) {
  outliers <- names(hk_zeros)
  cat("removed samples:", outliers, sep = "\n")  # housekeepers geometric mean is 0
  bm <- bm[, -hk_zeros]
}

# plot PCA ----------------------------------------------------------------

df <- bm %>%
  t() %>%
  as_tibble(rownames = "sample") %>%
  bind_rows(controls)

mat <- df %>%
  select(-label) %>%
  column_to_rownames("sample") %>%
  as.matrix()

pca <- mat %>%
  prcomp(scale. = T)

pca$x %>%
  as_tibble(rownames = "sample") %>%
  left_join(df) %>%
  ggplot(aes(PC1, PC2, color = label, group = sample)) +
  geom_point() +
  labs(title = "Biomarkers", subtitle = "normalized molecule counts")
ggplotly()

# plot heatmap ------------------------------------------------------------

mat[mat == 0] <- NA  # replace 0 with NA
pheatmap(log(t(mat)), scale = "row", main = "Biomarkers", treeheight_row = 0)
