# Infer TF activity based on mean expression of target genes
setwd("D:/GRN")

# install.packages("io", lib="D:/R-4.1.2/library")
library(io)

# input expression data matrix (rownames: gene symbols)
expr_tissue_median_gtex <- readRDS("C:/Users/Celin/Downloads/expr_tissue-median_gtex.rds")
expr <- expr_tissue_median_gtex$data

# input TF-target list of targets
# TP53
#   - CDKN1A
#   - ...
# RUNX1
#   - ...

tf.targets <- TF_list_DORO_cdkn1a; # 35 TFs and their target list

# Infer TF activity matrix (all samples by all TFs)
# target is from the list
tf.activity_list <- lapply(tf.targets,
    function(target)
      {
      target_name <- unlist(target)
        colMeans(as.matrix(expr[target_name, ]))
      }
)
# convert list of TF-target_mean to matrix
tf.activities <- do.call(rbind, tf.activity_list)

# @param target.expr    expression vector of target across N samples
# @param tf.activities  TF activity matrix (N samples by K TFs)
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

target <- "CDKN1A"
tf.candidates <- unique(TF_matrix_DORO_cdkn1a$source)

## == Fit a univariate linear model to infer regulatory effects =====================

ulinear_model <- function(expr, tf.activities, target){
  apply(tf.activities, MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
  coef(fit)
  })
}
ulinear_model(expr, tf.activities, target)

#To check whether the function output the same value
fit <- lm(expr[target,] ~ tf.activities["TFAP2A",]);
coef(fit) # NA???


## == Fit a multivariate linear model to infer regulatory effects ====================

# check before use: rownames(tf.activities) %in% rownames(expr)
mlinear_model <- function(expr, tf.activities, target){
  fit <- lm(expr[target, ] ~ t(tf.activities));
    coef(fit)
}
# call the mlinear function
mlinear_model(expr, tf.activities, target)


## == Fit a lasso model to infer regulatory effects ==================================

lasso_model <- function(expr, tf.activities, target){
  library(glmnet,lib="D:/R-4.1.2/library")
    x <- t(tf.activities)
    y <- expr[target,]
      fit <- glmnet(x, y, alpha = 1)
      coef(fit)
}
lasso_model(expr, tf.activities, target)
# In a Lasso model, the coefficients of some variables can be represented as a dot . 
# when the L1 penalty parameter (lambda) is set high enough such that the coefficients 
# of these variables are shrunk to zero by the L1 penalty.


