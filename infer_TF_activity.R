## Infer TF activity based on mean expression of target genes
## univariate, multivariate linear regression, lasso and mombf
##================================================================================

# install.packages("io")
library(io)

# input expression data matrix (rownames: gene symbols)
expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
expr <- expr_tissue_median_gtex$data
TF <- read.csv("data/TF.csv", header = FALSE)

# input TF-target list of targets
# TP53
#   - CDKN1A
#   - ...
# RUNX1
#   - ...

tf.targets <- tf_list_dorothea_cdkn1a; # 35 TFs and their target list

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

# @param expr           GTEX expression matrix across N samples
# @param tf.activities  TF activity matrix (N samples by K TFs)
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

target <- "CDKN1A"
tf.candidates <- unique(TF_matrix_DORO_cdkn1a$source)

## == Fit a univariate linear model to infer regulatory effects =====================

ulinear_model_II <- function(expr, tf.activities, target){
  apply(tf.activities, MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
  coef(fit)
  })
}
ulinear_model_II(expr, tf.activities, target)

#To check whether the function output the same value
fit <- lm(expr[target,] ~ tf.activities["TFAP2A",]);
coef(fit) # NA???


## == Fit a multivariate linear model to infer regulatory effects ====================

# check before use: rownames(tf.activities) %in% rownames(expr)
mlinear_model_II <- function(expr, tf.activities, target){
  fit <- lm(expr[target, ] ~ t(tf.activities));
    coef(fit)
}
# call the mlinear function
mlinear_model_II(expr, tf.activities, target)


## == Fit a lasso model to infer regulatory effects ==================================


# TODO: install the package again and into the central location
# restart and rerun the code when finish writing
lasso_model_II <- function(expr, tf.activities, target){
  library(glmnet,lib="D:/R-4.1.2/library")
    x <- t(tf.activities)
    y <- expr[target,]
      fit <- glmnet(x, y, alpha = 1)
      coef(fit)
}
lasso_model_II(expr, tf.activities, target)
# In a Lasso model, the coefficients of some variables can be represented as a dot . 
# when the L1 penalty parameter (lambda) is set high enough such that the coefficients 
# of these variables are shrunk to zero by the L1 penalty.

## == Fit a mombf model to infer regulatory effects ==================================

mombf_model_II <- function(expr, tf.activities, target){
  library(mombf,lib="D:/R-4.1.2/library")
  x <- t(tf.activities)
  y <- expr[target,]
    fit <-bestBIC(y ~ x)
    summary(fit)
}
mombf_model_II(expr, tf.activities, target)
# Error in modelSelection(..., priorCoef = bic(), priorDelta = modelunifprior(), :
# There are >1 constant columns in x (e.g. two intercepts)

x1 <- t(tf.activities)
y1 <- expr[target,]

# debug(function_name)
# TODO write 
# don't indent for no reason, only when have a curly bracket