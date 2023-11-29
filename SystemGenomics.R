setwd("Documents")
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
  select(c("Run",
           "Age",
           "Sample.Name",
           "Organ",
           "condition"))
expr <- expr[,meta$Run]

library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "mmusculus_gene_ensembl")
meta_genes <- getBM(attributes = c("ensembl_gene_id",
                                   "ensembl_gene_id_version",
                                   "hgnc_symbol",
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
avg_expr <- rowMeans(expr)
layout(matrix(1:2, nrow=1))
hist(avg_expr)
hist(log10(avg_expr + 1))
library(ggplot2)
ggplot(data.frame(avg_expr), aes(x=avg_expr)) + geom_histogram(bins = 50) + scale_x_continuous(breaks = c(0,1,10,100,1000,10000,20000), trans="log1p"
           , expand=c(0,0)) + scale_y_continuous(breaks = c(0,1), expand=c(0,0), trans="log1p") +
  theme_minimal()
num_det <- rowSums(expr > 0)
hist(num_det)
expressed <- rowMeans(expr > 0) >= 0.5 | rowMeans(expr) >= 1
#first way
#expr <- expr[which(expressed),]
#meta_genes <- meta_genes[which(expressed),]
#second way
meta_genes$expressed <- expressed
corr_pearson <- cor(log1p(expr[meta_genes$expressed,]))
corr_spearman <- cor(expr[meta_genes$expressed,], method = "spearman")
hcl_pearson <- hclust(as.dist(1 - corr_pearson))
hcl_spearman <- hclust(as.dist(1 - corr_spearman))
layout(matrix(1:2,nrow=1))
plot(hcl_pearson)
plot(hcl_spearman)
layout(matrix(1:2,nrow=1))
plot(hcl_spearman, labels = meta$Organ)
plot(hcl_spearman, labels = meta$condition)
pca <- prcomp(log1p(t(expr[meta_genes$expressed,])), center = TRUE, scale. = TRUE)
eigs <- pca$sdev^2
plot(1:length(eigs), eigs)
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
layout(matrix(1:2,nrow=1))
plot(hcl_spearman_highvar, labels = meta$Organ)
plot(hcl_spearman_highvar, labels = meta$condition)

pca_highvar <- prcomp(log1p(t(expr[meta_genes$highvar,])), center = TRUE, scale. = TRUE)
ggplot(data.frame(pca_highvar$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = Organ, shape = condition), size = 5)
library(sva)
expr_combat <- ComBat_seq(counts = expr,
                          batch = meta$Organ)

corr_spearman_combat <- cor(expr_combat[meta_genes$expressed,], method = "spearman")
hcl_spearman_combat <- hclust(as.dist(1 - corr_spearman_combat))
layout(matrix(1:2,nrow=1))
plot(hcl_spearman_combat, labels = meta$Organ)
plot(hcl_spearman_combat, labels = meta$condition)

pca_combat <- prcomp(log1p(t(expr_combat[meta_genes$expressed,])), center = TRUE, scale. = TRUE)
ggplot(data.frame(pca_combat$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = Organ, shape = condition), size = 5)
dat <- data.frame(y = log1p(as.numeric(expr["ENSMUSG00000000001.5",])), meta)
m1 <- lm(y ~ Organ + condition, data = dat)
m0 <- lm(y ~ Organ, data = dat)
test <- anova(m1, m0)
pval <- test$Pr[2]
library(pbapply)
pvals <- pbapply(expr[meta_genes$expressed,], 1, function(e){
  dat <- data.frame(y = log1p(e),
                    meta)
  m1 <- lm(y ~ Organ + condition, data = dat)
  m0 <- lm(y ~ condition, data = dat)
  test <- anova(m1, m0)
  pval <- test$Pr[2]
  return(unname(pval))
})
padj <- p.adjust(pvals, method = "bonferroni")
fc <- pbapply(expr[meta_genes$expressed,], 1, function(e){
  avg_layers <- tapply(log1p(e), meta$Organ, mean)
  return(exp(max(avg_layers) - min(avg_layers)))
})
DE_test <- function(expr,
                    cond,
                    ctrl = NULL,
                    covar = NULL,
                    padj_method = p.adjust.methods){
  pval_fc <- data.frame(t(pbapply(expr, 1, function(e){
    dat <- data.frame(y = log1p(e),
                      cond = cond)
    if (! is.null(covar))
      dat <- data.frame(dat, covar)
    
    m1 <- lm(y ~ ., data = dat)
    m0 <- lm(y ~ . - cond, data = dat)
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

res_DE <- DE_test(expr = expr[meta_genes$expressed,],
                  cond = meta$Organ,
                  covar = meta %>% dplyr::select(condition)) %>%
  tibble::rownames_to_column("gene")
res_DE <- res_DE %>%
  mutate(DE = padj < 0.1 & fc > 2)

library(ggrepel)
ggplot(res_DE, aes(x = log(fc), y = -log10(padj), col=DE, label=DEG)) +
  geom_point() +
  geom_text_repel() +
  geom_vline(xintercept=c(log(2), 0), col="#303030", linetype="dotted") +
  geom_hline(yintercept=-log10(0.1), col="#303030", linetype="dotted") +
  scale_color_manual(values=c("#909090", "red")) +
  theme_minimal()
