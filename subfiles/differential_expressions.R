library(tidyverse)
library(broom)
library(DESeq2, lib.loc="/data/users/mccallke0364/.conda/envs/final/lib/R/library/")

counts_df <- read_tsv("/data/users/mccallke0364/Final_Project/Sm_Mira_IvT/subfiles/star_counts_2.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))
write.csv(counts_df, "plots/counts_df.csv")

counts_summary <- counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting_sub_2/dedup_sub/star.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)
write.csv(counts_summary, "plots/counts_summary.csv")

sample_summary <- counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting_sub_2/dedup_sub/star.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)
write.csv(sample_summary, "plots/sample_summary.csv")
genes_to_remove = sample_summary$Geneid

counts_filt <- counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

counts_m <- counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(counts_m) <- counts_filt$Geneid
write.csv(counts_m, "plots/counts_m.csv")

dists <- dist(t(counts_m))

dists_df <- as.matrix(dists) |>
    as_tibble(rownames = 'sample')
write.csv(dists_df, "plots/dists_df.csv")

dist_plot <- dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
ggsave("plots/dist_plot.png")

pca_fit <- t(log10(counts_m + 1)) |> 
  prcomp(scale = TRUE)
pca_fit

pca_fit |>
  augment(t(counts_m)) |>
  dplyr::rename(sample = .rownames) |>
  mutate(sample = str_remove(sample, '[0-9]')) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = sample)) + 
  geom_point(size = 4)
ggsave("plots/pca_fit.png")

metadata <- data.frame(sample_id = colnames(counts_m)) |>
    mutate(sample = str_sub(sample_id, 1, 3),
           rep = str_sub(sample_id, 5))
rownames(metadata) <- metadata$sample_id
metadata <- select(metadata, -sample_id)
metadata

all(rownames(metadata) == colnames(counts_m))


dds2 <- DESeqDataSetFromMatrix(countData = counts_m,
                               colData = metadata,
                               design = ~ sample)
dds2 <- DESeq(dds2)

result <- results(dds2)

volcano_data <- as_tibble(result, rownames = "gene_id")

volcano_plot <- volcano_data |> 
    ggplot( aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2) +
    theme_minimal()
volcano_plot
ggsave("plots/volcano_plot.png")

volcano_plot_2 <- volcano_data |> 
    mutate(colour.grp = ifelse(log2FoldChange > 1.6, yes = "above", no = "nothing")) |>
    mutate(colour.grp = ifelse(log2FoldChange < -1.6, yes = "below", no = colour.grp))|>
    ggplot( aes(x = log2FoldChange, y = -log10(pvalue), colour=colour.grp)) +
    scale_colour_manual(values = c("green", "red", "grey"),
                        labels = c("Upregulated","Downregulated", "Not significant"))+
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2) +
    ggtitle("Gene Expression")+
    theme_gray()
volcano_plot_2
ggsave("plots/volcano_plot_2.png")

vsd <- varianceStabilizingTransformation(dds2)
write.csv(assay(vsd), "counting_sub_2/vsd.csv")
pca_desq <- plotPCA(vsd, intgroup = "sample")
ggsave("plots/pca_desq.png")
write_csv(volcano_data, "plots/volcano_data.csv")

counts_df <- read_tsv("/data/users/mccallke0364/Final_Project/Sm_Mira_IvT/subfiles/star_counts_2.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))

counts_save <- counts_df |>
    select(Geneid, contains("star.bam")) |>
    rename_with(~str_remove(., "counting_sub_2/star.bam:"), everything()) |>
    rename(Geneid = "gene_id")

write_csv(counts_save, "plots/counts_clean.csv")



















