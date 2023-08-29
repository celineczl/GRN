## class I mombf (non-local prior variable selection) model to infer regulatory effect
## assume TF activity = expression level

## ======== mombf variable selection model ====================
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
expr <- expr_tissue_median_gtex$data
target <- "CDKN1A"
tf.candidates <- readRDS("data/tf.candidates.rds")

# https://cran.r-project.org/web/packages/mombf/mombf.pdf
# https://cran.r-project.org/web/packages/mombf/vignettes/mombf.pdf
# https://github.com/davidrusi/mombf

## sessionInfo()
# install.packages("rlang")
library(rlang) # rlang_1.1.1

# install.packages("mombf")
library(mombf)

# input matrix x is TF expression data
x <- t(expr[tf.candidates,])
# response variable is the target gene expression
y <- expr[target,]

msfit <- modelSelection(y, x, family="normal", priorCoef=momprior())
sampl <- rnlp(msfit=msfit)
beta <- colMeans(sampl)[-ncol(sampl)]
# remove phi in the last column

hist(sampl[, "AR"], breaks=100)
hist(sampl[, "GLI2"], breaks=100)
hist(sampl[, "TP73"], breaks=100)

mombf_I_coef <- as.matrix(beta)
saveRDS(mombf_I_coef, "results/mombf_I_coef.rds")
