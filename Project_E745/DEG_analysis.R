library( "DESeq2" )
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

E745 <- read.table("/Users/chenjinglian/Desktop/UU25 master/Genome_analysis/Project_E745/05_DEG/counts_E745.txt",
           sep = "\t",header = T,row.names = 2)
E745 <- E745[,-1]
hist(log10(gene_total_counts + 1),
     breaks = 50,
     main = "Distribution of log10 read counts per gene",
     xlab = "log10(total counts + 1)",
     ylab = "Number of genes")
E745 <- E745[!grepl("^__", rownames(E745)), ]
E745 <- E745[rowSums(E745) > 0, ]
# get the sample names from the count_data matrix
SampleName <- c(colnames(E745))
# specify the conditions for each sample
# In my sample names, DMSO_1 and DMSO_2 are control replicated and TCPMOH_1 and TCPMOH_2 are treated replicates
condition <- c("Serum", "Serum", "Serum", "BH", "BH", "BH")
# generate the metadata data frame
meta_data <- data.frame(SampleName, condition)
# make the sample name the row id
meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var="SampleName")

# create deseq data set object
dds <- DESeqDataSetFromMatrix(countData = E745,
                              colData = meta_data,
                              design = ~ condition)
# filter any counts less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# set the reference to be the BH(control, full-nutrient)
dds$condition <- relevel(dds$condition, ref = 'BH')
# get normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
# Run Differential Expression
dds <- DESeq(dds)
res <- results(dds)
summary(res)
# Set p-value < 0.05
res05 <- res[!is.na(res$padj) & res$padj < 0.05, ]
summary(res05)
# Visualize result
data <- data.frame(res05)
# PCA
rld <- rlog(dds)
plotPCA(rld)
# add an additional column that identifies a gene as unregulated, downregulated, or unchanged
# note the choice of pvalue and log2FoldChange cutoff. 
data <- data %>%
  mutate(
    Expression = case_when(log2FoldChange > 2 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < -2 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

top <- 10
# we are getting the top 10 up and down regulated genes by filtering the column Up-regulated and Down-regulated and sorting by the adjusted p-value. 
top_genes <- bind_rows(
  data %>%
    filter(Expression == 'Up-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top),
  data %>%
    filter(Expression == 'Down-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top)
)
# create a datframe just holding the top 10 genes
Top_Hits = head(arrange(data,pvalue),10)
Top_Hits

# Volcano Plot
data$label = if_else(rownames(data) %in% rownames(Top_Hits), rownames(data), "")
# with labels for top 10 sig overall
p1 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf)
ggsave("/Users/chenjinglian/Desktop/UU25 master/Genome_analysis/Project_E745/05_DEG/p1.pdf", 
       plot = p1, width = 9, height = 6)

# Annotation
anno <- read.table("/Users/chenjinglian/Desktop/UU25 master/Genome_analysis/Project_E745/05_DEG/E745.emapper_annotation_simplified.tsv",
                   sep = "\t",header = T)
anno$gene_label <- ifelse(anno$Preferred_name == "-" | anno$Preferred_name == "" | is.na(anno$Preferred_name),
                          anno$query,
                          anno$Preferred_name)
anno <- anno[,-c(1,2)]
data$gene_label <- rownames(data)
# Merge
deg_anno <- merge(data, anno, by = "gene_label", all.x = TRUE)

# Seperate up-regulated/down-regulated
up_deg <- subset(deg_anno, !is.na(padj) & padj < 0.05 & log2FoldChange > 1)
down_deg <- subset(deg_anno, !is.na(padj) & padj < 0.05 & log2FoldChange < -1)

# Seperate bg/up/down
bg_cog <- deg_anno %>%
  select(gene_label, COG_category) %>%
  filter(!is.na(COG_category), COG_category != "-", COG_category != "") %>%
  mutate(COG_category = strsplit(COG_category, "")) %>%
  unnest(COG_category)

up_cog <- up_deg %>%
  select(gene_label, COG_category) %>%
  filter(!is.na(COG_category), COG_category != "-", COG_category != "") %>%
  mutate(COG_category = strsplit(COG_category, "")) %>%
  unnest(COG_category)

down_cog <- down_deg %>%
  select(gene_label, COG_category) %>%
  filter(!is.na(COG_category), COG_category != "-", COG_category != "") %>%
  mutate(COG_category = strsplit(COG_category, "")) %>%
  unnest(COG_category)

# Fisher's exact test
cog_enrichment <- function(target_df, background_df, target_genes, background_genes) {
  categories <- sort(unique(background_df$COG_category))
  
  results <- lapply(categories, function(cat) {
    target_with_cat <- sum(unique(target_df$gene_label[target_df$COG_category == cat]) %in% target_genes)
    target_without_cat <- length(target_genes) - target_with_cat
    
    bg_with_cat <- sum(unique(background_df$gene_label[background_df$COG_category == cat]) %in% background_genes)
    bg_without_cat <- length(background_genes) - bg_with_cat
    
    mat <- matrix(c(target_with_cat,
                    target_without_cat,
                    bg_with_cat - target_with_cat,
                    bg_without_cat - target_without_cat),
                  nrow = 2, byrow = TRUE)
    
    ft <- fisher.test(mat, alternative = "greater")
    
    data.frame(
      COG_category = cat,
      target_with_cat = target_with_cat,
      background_with_cat = bg_with_cat,
      pvalue = ft$p.value
    )
  })
  
  res <- bind_rows(results)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$padj, res$pvalue), ]
  res
}
up_cog_res <- cog_enrichment(
  target_df = up_cog,
  background_df = bg_cog,
  target_genes = unique(up_deg$gene_label),
  background_genes = unique(deg_anno$gene_label)
)

down_cog_res <- cog_enrichment(
  target_df = down_cog,
  background_df = bg_cog,
  target_genes = unique(down_deg$gene_label),
  background_genes = unique(deg_anno$gene_label)
)
subset(up_cog_res, padj < 0.05)
subset(down_cog_res, padj < 0.05)
# Visualization
up_plot <- up_cog_res %>%
  filter(padj < 0.05) %>%
  mutate(direction = "Up",
         score = -log10(padj))
down_plot <- down_cog_res %>%
  filter(padj < 0.05) %>%
  mutate(direction = "Down",
         score = log10(padj)) 
cog_plot <- bind_rows(up_plot, down_plot)
ggplot(cog_plot, aes(x = reorder(COG_category, score), y = score, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "firebrick3", "Down" = "dodgerblue3")) +
  xlab("COG category") +
  ylab(expression(paste("Enrichment score (", -log[10], " adjusted p-value)"))) +
  ggtitle("COG enrichment in up- and down-regulated genes")

# GO enrichment
background_genes <- unique(deg_anno$gene_label)
bg_go <- deg_anno %>%
  select(gene_label, GOs) %>%
  filter(!is.na(GOs), GOs != "-", GOs != "") %>%
  separate_rows(GOs, sep = ",")
up_go <- up_deg %>%
  select(gene_label, GOs) %>%
  filter(!is.na(GOs), GOs != "-", GOs != "") %>%
  separate_rows(GOs, sep = ",")
down_go <- down_deg %>%
  select(gene_label, GOs) %>%
  filter(!is.na(GOs), GOs != "-", GOs != "") %>%
  separate_rows(GOs, sep = ",")
 
# 
go_enrichment <- function(target_df, background_df, target_genes, background_genes) {
  go_terms <- sort(unique(background_df$GOs))
  
  results <- lapply(go_terms, function(term) {
    target_with_term <- sum(unique(target_df$gene_label[target_df$GOs == term]) %in% target_genes)
    target_without_term <- length(target_genes) - target_with_term
    
    bg_with_term <- sum(unique(background_df$gene_label[background_df$GOs == term]) %in% background_genes)
    bg_without_term <- length(background_genes) - bg_with_term
    
    mat <- matrix(c(target_with_term,
                    target_without_term,
                    bg_with_term - target_with_term,
                    bg_without_term - target_without_term),
                  nrow = 2, byrow = TRUE)
    
    ft <- fisher.test(mat, alternative = "greater")
    
    data.frame(
      GO = term,
      target_with_term = target_with_term,
      background_with_term = bg_with_term,
      pvalue = ft$p.value
    )
  })
  
  res <- bind_rows(results)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$padj, res$pvalue), ]
  res
}
up_go_res <- go_enrichment(
  target_df = up_go,
  background_df = bg_go,
  target_genes = unique(up_deg$gene_label),
  background_genes = unique(deg_anno$gene_label)
)

down_go_res <- go_enrichment(
  target_df = down_go,
  background_df = bg_go,
  target_genes = unique(down_deg$gene_label),
  background_genes = unique(deg_anno$gene_label)
)
sig_up_go <- subset(up_go_res, pvalue < 0.05)
sig_down_go <- subset(down_go_res, padj< 0.1)[1:11,]

# Visulaize
library(GO.db)
library(AnnotationDbi)
library(dplyr)
go_info <- AnnotationDbi::select(
  GO.db,
  keys = unique(sig_down_go$GO),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
)
sig_down_go2 <- sig_down_go2 %>%
  filter(!is.na(TERM), !is.na(ONTOLOGY))
sig_down_go2 <- sig_down_go %>%
  left_join(go_info, by = c("GO" = "GOID"))
sig_down_go2$label <- paste0(sig_down_go2$TERM, " (", sig_down_go2$ONTOLOGY, ")")
sig_down_go2$label_short <- ifelse(
  nchar(sig_down_go2$label) > 50,
  paste0(substr(sig_down_go2$label, 1, 50), "..."),
  sig_down_go2$label
)
ggplot(sig_down_go2, aes(x = reorder(label_short, -log10(padj)), y = -log10(padj), fill = ONTOLOGY)) +
  geom_col() +
  coord_flip() +
  xlab("GO term") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Enriched GO terms in down-regulated genes")




