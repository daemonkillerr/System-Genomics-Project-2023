#setwd("Documents")
samples <- list.files("rsem")
expr <- sapply(samples, function(sample){
  file <- paste0("rsem/",sample,"/",sample,".genes.results")
  quant <- read.csv(file, sep="\t", header=T)
  tpm <- setNames(quant$TPM, quant$gene_id)
  return(tpm)
})
library(dplyr)
meta <-
  read.csv("SraRunTable.txt", header=T) %>%
  dplyr::select(c("Run",
           "Age",
           "Sample.Name",
           "Organ",
           "condition"))
expr <- expr[,meta$Run]

#sra_header <- read.csv("SraRunTable.txt", header=T)

library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "mmusculus_gene_ensembl")
meta_genes <- getBM(attributes = c("ensembl_gene_id",
                                   "ensembl_gene_id_version",
                                   "mgi_symbol",
                               "description",
                               "chromosome_name",
                               "start_position",
                               "end_position",
                               "strand"),
                    filters = "ensembl_gene_id_version",
                    values = rownames(expr),
                    mart = ensembl) %>%
  right_join(data.frame(ensembl_gene_id_version = rownames(expr)),
             by = "ensembl_gene_id_version") %>%
  distinct(ensembl_gene_id_version, .keep_all = TRUE)

dim(expr)
#expr2 <- expr[meta_genes$ensembl_gene_id_version]
#dim(expr2)
avg_expr <- rowMeans(expr)
layout(matrix(1, nrow=1))

#give histograms good titles and axis names
hist(avg_expr, main = "Histogram of average gene expression across all samples", xlab = "Average Gene expression", ylab = "Frequency")
hist(log10(avg_expr + 1), main = "Histogram of average gene expression across all samples", xlab = "Average Gene expression (log)", ylab = "Frequency")
library(ggplot2)
# make it clear on axis label if lgged (bot y and x if true)
ggplot(data.frame(avg_expr), aes(x=avg_expr)) + geom_histogram(bins = 50) + scale_x_continuous(breaks = c(0,1,10,100,1000,10000,20000), trans="log1p"
           , expand=c(0,0)) + scale_y_continuous(breaks = c(0,500), expand=c(0,0), trans="log1p") +
  theme_minimal() + labs(title = "Histogram of average gene expression across all samples", x = "Average Gene expression (log)", y = "Frequency (log)")
num_det <- rowSums(expr > 0)
hist(num_det, main = "Histogram of number of detected genes across all samples", xlab = "Number of detected genes", ylab = "Frequency")
expressed <- rowMeans(expr > 0) >= 0.5 | rowMeans(expr) >= 1
#first way
#expr <- expr[which(expressed),]
#meta_genes <- meta_genes[which(expressed),]
#second way
meta_genes$expressed <- expressed
# add an "id" column (just row number) to meta
meta$id <- 1:nrow(meta)
corr_pearson <- cor(log1p(expr[meta_genes$expressed,]))
corr_spearman <- cor(expr[meta_genes$expressed,], method = "spearman")
hcl_pearson <- hclust(as.dist(1 - corr_pearson))
hcl_spearman <- hclust(as.dist(1 - corr_spearman))
layout(matrix(1,nrow=1))
plot(hcl_pearson, paste(meta$Organ, meta$condition, meta$id, sep = " "))
plot(hcl_spearman, paste(meta$Organ, meta$condition, meta$id, sep = " "))
layout(matrix(1,nrow=1))
plot(hcl_spearman, paste(meta$Organ, meta$condition, meta$id, sep = " "))
plot(hcl_spearman, paste(meta$Organ, meta$condition, meta$id, sep = " "))
pca <- prcomp(log1p(t(expr[meta_genes$expressed,])), center = TRUE, scale. = TRUE)
eigs <- pca$sdev^2
# add proper axis labels and a title
plot(1:length(eigs), eigs, main = "Scree plot of PCA", xlab = "Principal component", ylab = "Eigenvalue")
ggplot(data.frame(pca$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = Organ, shape = condition), size = 5)
estimate_variability <- function(expr){
  means <- apply(expr, 1, mean)
  vars <- apply(expr, 1, var)
  cv2 <- vars / means^2
  minMeanForFit <- unname(median(means[which(cv2 > 0.3)]))
  useForFit <- means >= minMeanForFit
  fit <- glm.fit(x = cbind(a0 = 1, a1tilde = 1/means[useForFit]),
                 y = cv2[useForFit],
                 family = Gamma(link = "identity"))
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])
  xg <- exp(seq(min(log(means[means>0])), max(log(means)), length.out=1000))
  vfit <- a1/xg + a0
  df <- ncol(expr) - 1
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  res <- data.frame(mean = means,
                    var = vars,
                    cv2 = cv2,
                    useForFit = useForFit,
                    pval = pval,
                    padj = p.adjust(pval, method="BH"),
                    row.names = rownames(expr))
  return(res)
}
var_genes <- estimate_variability(expr[meta_genes$expressed,])
meta_genes$highvar <- meta_genes$ensembl_gene_id_version %in% rownames(var_genes)[which(var_genes$padj < 0.01)]
corr_spearman_highvar <- cor(expr[meta_genes$highvar,], method = "spearman")
hcl_spearman_highvar <- hclust(as.dist(1 - corr_spearman_highvar))
layout(matrix(1,nrow=1))
plot(hcl_spearman_highvar, paste(meta$Organ, meta$condition, meta$id, sep = " "))
plot(hcl_spearman_highvar, paste(meta$Organ, meta$condition, meta$id, sep = " "))

pca_highvar <- prcomp(log1p(t(expr[meta_genes$highvar,])), center = TRUE, scale. = TRUE)
eigs_highvar <- pca_highvar$sdev^2
plot(1:length(eigs_highvar), eigs_highvar, main = "Scree plot of PCA using only high variance genes", xlab = "Principal component", ylab = "Eigenvalue")

ggplot(data.frame(pca_highvar$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = Organ, shape = condition), size = 5)
library(sva)
expr_combat <- ComBat_seq(counts = expr,
                          batch = meta$Organ)

corr_spearman_combat <- cor(expr_combat[meta_genes$expressed,], method = "spearman")
hcl_spearman_combat <- hclust(as.dist(1 - corr_spearman_combat))
layout(matrix(1:2,nrow=1))
plot(hcl_spearman_combat, paste(meta$Organ, meta$condition, meta$id, sep = " "))
plot(hcl_spearman_combat, paste(meta$Organ, meta$condition, meta$id, sep = " "))

pca_combat <- prcomp(log1p(t(expr_combat[meta_genes$expressed,])), center = TRUE, scale. = TRUE)
eigs_combat <- pca_combat$sdev^2
plot(1:length(eigs_combat), eigs_combat, main = "Scree plot of PCA after accounting for batch effects", xlab = "Principal component", ylab = "Eigenvalue")
ggplot(data.frame(pca_combat$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = Organ, shape = condition), size = 5)

#Commented out version of the DE analysis that uses just 1 gene and compares the wrong models
#dat <- data.frame(y = log1p(as.numeric(expr["ENSMUSG00000000001.5",])), meta)
#m1 <- lm(y ~ Organ + condition, data = dat)
#m0 <- lm(y ~ Organ, data = dat)
#test <- anova(m1, m0)
#pval <- test$Pr[2]
library(pbapply)

# pvals <- pbapply(expr[meta_genes$expressed,], 1, function(e){
#   dat <- data.frame(y = log1p(e),
#                     meta)
#   m1 <- lm(y ~ Organ + condition, data = dat)
#   m0 <- lm(y ~ condition, data = dat)
#   test <- anova(m1, m0)
#   pval <- test$Pr[2]
#   return(unname(pval))
# })
# padj <- p.adjust(pvals, method = "bonferroni")
# fc <- pbapply(expr[meta_genes$expressed,], 1, function(e){
#   avg_layers <- tapply(log1p(e), meta$Organ, mean)
#   return(exp(max(avg_layers) - min(avg_layers)))
# })

meta$condition <- factor(meta$condition, levels = c("steady_state", "LCMVd8")) #TODO see if I reverse the ordering

DE_test <- function(expr,
                    cond,
                    ctrl = NULL,
                    covar = NULL,
                    padj_method = "holm"){
  pval_fc <- data.frame(t(pbapply(expr, 1, function(e){
    dat <- data.frame(y = log1p(e),
                      cond = cond)
    if (! is.null(covar))
      dat <- data.frame(dat, covar)
    
    m1 <- lm(y ~ ., data = dat)
    m0 <- lm(y ~ . - cond, data = dat)
    # print features used in m1 and m0
    
    test <- anova(m1, m0)
    pval <- test$Pr[2]
    
    avgs <- tapply(log1p(e), cond, mean)
    if (! is.null(ctrl) && sum(cond %in% ctrl) > 0){
      fc <- exp(max(avgs[names(avgs) != ctrl]) - avgs[ctrl])
    } else{
      fc <- exp(max(avgs) - min(avgs))
    }
    
    return(c(pval = unname(pval), fc = unname(fc)))
  })), row.names = rownames(expr))
  padj <- p.adjust(pval_fc$pval, method = padj_method)
  return(data.frame(pval_fc, "padj" = padj)[,c("pval","padj","fc")])
}

#Ab hier update on github

res_DE <- DE_test(expr = expr[meta_genes$expressed,],
                  cond = meta$condition,
                  covar = meta %>% dplyr::select(Organ),
                  padj_method = "holm") %>%
  tibble::rownames_to_column("gene")
res_DE <- res_DE %>%
  mutate(DE = padj < 0.1 & fc > 2) %>% 
  mutate(DEG = ifelse(DE, meta_genes$mgi_symbol, NA))


library(ggrepel)
ggplot(res_DE, aes(x = log(fc), y = -log10(padj), col = DE, label = DEG )) +
  geom_point() +
  geom_text_repel() +
  geom_vline(xintercept=c(log(2), 0), col="#303030", linetype="dotted") +
  geom_hline(yintercept=-log10(0.1), col="#303030", linetype="dotted") +
  scale_color_manual(values=c("#909090", "red")) +
  theme_minimal()

#Ab here needs to be checked and run and adjusted
#library(biomaRt)
#ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

tx2gene <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version"), 
                 filters = "ensembl_gene_id_version", values = rownames(expr), mart = ensembl) %>% 
  dplyr::select(ensembl_transcript_id_version, ensembl_gene_id_version)


#BiocManager::install("tximport")

library(tximport)
samples <- list.files("rsem")
files <- file.path("rsem", samples, paste0(samples,".isoforms.results"))
#file <- paste0("rsem/",sample,"/",sample,".genes.results")
txi <- tximport(files, type = "rsem", tx2gene = tx2gene)
#BiocManager::install("DESeq2")

library(DESeq2)

dds <- DESeqDataSetFromTximport(txi, colData = meta, 
                                 design = ~ Organ + condition)

dds_filtered <- dds[intersect(rownames(expr)[meta_genes$expressed], rownames(dds)),]
dds_filtered <- DESeq(dds_filtered, test="LRT", reduced= ~ Organ)
res_DESeq2 <- results(dds_filtered)

cor(res_DESeq2$padj, 
    res_DE %>% filter(gene %in% rownames(res_DESeq2)) %>% 
      pull(padj), method="spearman", use="complete.obs")

layout(matrix(1, nrow=1))
# Adjusting the margins (if necessary)
# par(mar=c(5.1, 4.1, 4.1, 2.1)) # default
par(mar=c(5.1, 4.1, 4.1, 2.1) + 0.1) # Increase if labels are cut off
plot(-log10(res_DESeq2$pvalue), 
     -log10(res_DE %>% filter(gene %in% rownames(res_DESeq2)) %>% pull(pval)),
     xlab = "-log10(pval DESeq2)", ylab = "-log10(pval DE)", pch=16)

smoothScatter(-log10(res_DESeq2$pvalue), -log10(res_DE %>% filter(gene %in% rownames(res_DESeq2)) %>% pull(pval)),
              xlab = "-log10(pval DESeq2)", ylab = "-log10(pval DE)", pch=16)

table(p.adjust(res_DESeq2$pvalue, method="bonferroni") < 0.1, res_DE %>% filter(gene %in% rownames(res_DESeq2)) %>% pull(padj) < 0.1) 

fold_change_threshold <- 1  # Adjust this threshold as needed

# get rid of NaNs
#logical_vector <- ifelse(is.na(res_DESeq2$log2FoldChange), FALSE, res_DESeq2$log2FoldChange > fold_change_threshold)
#res_DESeq2 <- res_DESeq2[logical_vector,]

# there are NaNas in log2FoldChange. Make these 0s
res_DESeq2$log2FoldChange[is.na(res_DESeq2$log2FoldChange)] <- 0

# Identify up-regulated genes
upregulated_genes <- rownames(res_DESeq2[res_DESeq2$log2FoldChange > fold_change_threshold, ])

# Identify down-regulated genes
downregulated_genes <- rownames(res_DESeq2[res_DESeq2$log2FoldChange < -fold_change_threshold, ])


# Optional: Create a column in res_DESeq2 to mark up-regulated and down-regulated genes
res_DESeq2$regulation <- ifelse(rownames(res_DESeq2) %in% upregulated_genes, "Up-regulated",
                                ifelse(rownames(res_DESeq2) %in% downregulated_genes, "Down-regulated", "Not-regulated"))

# Summary of the number of up-regulated and down-regulated genes
summary(res_DESeq2$regulation)
print(res_DESeq2$regulation)
# Optional: Plot log2 fold changes
#making regulation a factor so that I can use it as a colour
res_DESeq2$regulation <- as.factor(res_DESeq2$regulation)
plot(res_DESeq2$log2FoldChange, pch = 16, col = res_DESeq2$regulation)
abline(h = c(-fold_change_threshold, fold_change_threshold), lty = 2, col = "red")
legend("topright", legend = c("Up-regulated", "Down-regulated"), col = c("green", "black"), pch = 16)

#Grouping
#Case-control
#Up-/Downregulated

#For multiple comparison
group1_condition <- (meta$Organ == 'Brain') & (meta$condition == 'steady_state')
group2_condition <- (meta$Organ == 'Brain') & (meta$condition == 'LCMVd8')
group3_condition <- (meta$Organ == 'Liver') & (meta$condition == 'steady_state')
group4_condition <- (meta$Organ == 'Liver') & (meta$condition == 'LCMVd8')

DEG <- res_DE$gene[res_DE$DE]
avg_expr <- sapply(sort(unique(meta$Organ)), 
                   function(Organ) rowMeans(expr[,which(meta$Organ == Organ)])
                   )

#All are maximally expressed in the brain
# max_layer_DEG <- setNames(colnames(avg_expr)[apply(avg_expr[DEG,], 1, which.max)], DEG)
# avg_expr_DEG_list <- tapply(names(max_layer_DEG), max_layer_DEG, function(x) avg_expr[x,])
# 
# scaled_expr_DEG_list <- lapply(avg_expr_DEG_list, function(x) t(scale(t(x))))
# layout(matrix(1, nrow = 2))
# par(mar=c(3,3,3,3))
# Organ=c('Brain','Liver')
# for(i in 1:2) boxplot(avg_expr[,i], main = paste0(Organ[i], " (", nrow(scaled_expr_DEG_list), ")"))
# boxplot(avg_expr[,1], main = paste0(Organ[1], " (", nrow(scaled_expr_DEG_list), ")"))
# boxplot(avg_expr[,2], main = paste0(Organ[2], " (", nrow(scaled_expr_DEG_list), ")"))


# hcl_spearman <- hclust(as.dist(1 - corr_spearman))
# avg_expr <- sapply(sort(unique(meta$Organ)), function(Organ) rowMeans(expr[,which(meta$Organ == Organ)]))
# corr_DEG <- cor(avg_expr[res_DE$gene[res_DE$DE],], method = "spearman")
# mgi_DEG <- hclust(as.dist(1 - corr_DEG), method = "complete")
# #Doesn't work
# plot(hcl_spearman, labels = FALSE)
# 
# 
# library(gplots)
# heatmap.2(corr_DEG, Rowv = as.dendrogram(mgi_DEG), Colv = as.dendrogram(mgi_DEG), trace = "none", scale = "none", labRow = NA, labCol = NA)
# 
# install.packages("viridis")
# library(viridis)
# heatmap.2(corr_DEG, Rowv = as.dendrogram(mgi_DEG), Colv = as.dendrogram(mgi_DEG), trace = "none", scale = "none", labRow = NA, labCol = NA, col = viridis)
# 
# cl_DEG <- cutree(mgi_DEG, k = 15)
# heatmap.2(corr_DEG, Rowv = as.dendrogram(mgi_DEG), Colv = as.dendrogram(mgi_DEG), trace = "none", scale = "none", labRow = NA, labCol = NA, col = viridis, ColSideColors = rainbow(15)[cl_DEG])
# 
# avg_expr <- sapply(sort(unique(meta$Organ)), function(Organ) rowMeans(expr[,which(meta$Organ == Organ)]))
# avg_expr_DEG_list <- tapply(names(cl_DEG), cl_DEG, function(x) avg_expr[x,])
# scaled_expr_DEG_list <- lapply(avg_expr_DEG_list, function(x) t(scale(t(x))))
# layout(matrix(1:15, nrow = 3, byrow = T))
# par(mar=c(3,3,3,3))
# for(cl in 1:15) 
#   boxplot(scaled_expr_DEG_list[[cl]], main = paste0(cl, " (", nrow(scaled_expr_DEG_list[[cl]]), ")"))
# 
# unique(cl_DEG[hcl_DEG$order])
# 
# layout(matrix(1:15, nrow = 3, byrow = T))
# par(mar=c(3,3,3,3))
# for(layer in unique(cl_DEG[hcl_DEG$order])) 
#  boxplot(scaled_expr_DEG_list[[Organ]], main = paste0(layer, " (", nrow(scaled_expr_DEG_list[[Organ]]), ")"))

# need < 3000 rows to work in DAVID so should said the thresholding above higher
combined_genes <- c(upregulated_genes, downregulated_genes)
DE_genes_data <- meta_genes[meta_genes$ensembl_gene_id_version %in% combined_genes, "ensembl_gene_id"]

write.table(DE_genes_data, file = "genes_C2.txt", quote = F, row.names = F, col.names = F)

write.table(meta_genes[meta_genes$expressed, "ensembl_gene_id"], file = "genes_expressed.txt", quote = F, row.names = F, col.names = F)


scores <- setNames(sign(log(res_DE$fc)) * (-log10(res_DE$pval)), setNames(meta_genes$ensembl_gene_id, meta_genes$ensembl_gene_id_version)[res_DE$gene])
scores_ordered <- sort(scores, decreasing=T)

#instead calcuating with the DESeq2 model
scores <- setNames(sign(res_DESeq2$log2FoldChange) * (-log10(res_DESeq2$pvalue)), setNames(meta_genes$ensembl_gene_id, meta_genes$ensembl_gene_id_version)[rownames(res_DESeq2)])
scores_ordered <- sort(scores, decreasing=T)


#install.packages("msigdbr")
library(msigdbr)
# Using C7 since it is genes in immune pathways. Could also use H for hallmark groupings
genesets_celltype <- msigdbr(species = "Mus musculus", category = "H")
genesets_celltype_list <- tapply(genesets_celltype$ensembl_gene, genesets_celltype$gs_name, list)

#BiocManager::install("fgsea")
library(fgsea)
fgsea_kegg <- fgsea(pathways = genesets_celltype_list, stats = scores_ordered,
minSize = 15, maxSize = 500)

fgsea_kegg[order(NES,decreasing=T),][1:10,1:7] 

#pick a gene from above to go here
plotEnrichment(genesets_celltype_list[["ZHONG_PFC_C2_THY1_POS_OPC"]], scores_ordered) + labs(title="ZHONG pathway")

fgsea_kegg[order(NES,decreasing=F),][1:10,1:7] 


# lets do a bunch of this again but look for DE differences in liver vs brain (so we now consider a model with an interaction term)
dds <- DESeqDataSetFromTximport(txi, colData = meta, 
                                design = ~ Organ + condition + Organ:condition)

dds_filtered <- dds[intersect(rownames(expr)[meta_genes$expressed], rownames(dds)),]
dds_filtered <- DESeq(dds_filtered, test="LRT", reduced= ~ Organ + condition)
res_DESeq2 <- results(dds_filtered)
summary(res_DESeq2)


fold_change_threshold <- 1  # Adjust this threshold as needed
scores <- setNames(sign(res_DESeq2$log2FoldChange) * (-log10(res_DESeq2$pvalue)), setNames(meta_genes$ensembl_gene_id, meta_genes$ensembl_gene_id_version)[rownames(res_DESeq2)])
scores_ordered <- sort(scores, decreasing=T)
genesets_celltype <- msigdbr(species = "Mus musculus", category = "H")
genesets_celltype_list <- tapply(genesets_celltype$ensembl_gene, genesets_celltype$gs_name, list)
fgsea_kegg <- fgsea(pathways = genesets_celltype_list, stats = scores_ordered,
                    minSize = 15, maxSize = 500)

fgsea_kegg[order(NES,decreasing=T),][1:10,1:7] 
fgsea_kegg[order(NES,decreasing=F),][1:10,1:7] 
#Some of the pathways we get here eg HALLMARK_COAGULATION make sense

