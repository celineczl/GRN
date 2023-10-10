## Infer TF activity based on mean expression of target genes
## univariate, multivariate linear regression, lasso and mombf
##================================================================================

# install.packages("io")
library(io)

# input TF-target list of targets
# TP53
#   - CDKN1A
#   - ...
# RUNX1
#   - ...

tf.targets <- readRDS("data/tf.list.rds"); # 156 TFs and their target list

# Infer TF activity matrix (all samples by all TFs)
# target is from the list
tf.activity_list <- lapply(tf.targets,
                           function(value)
                           {
                             target_name <- unlist(value)
                             # set drop dimension = FALSE, 
                             # or colMeans() will average the row expr of CDKN1A when no other targets
                             colMeans(as.matrix(expr[target_name, , drop=FALSE]))
                           }
)
# convert list of TF-target_mean to matrix
tf.activities_matrix <- do.call(rbind, tf.activity_list)
saveRDS(tf.activities_matrix, "data/tf.activities_matrix.rds")

## need to remove three TFs with only target "CDKN1A"
# identify all rows that have exactly same expression as CDKN1A
tf_cdkn1a_only <- which(apply(tf.activities_matrix, 1, function(row) all(row == expr["CDKN1A", ])))

# remove those rows to avoid assigning coefficient ~1 to those TFs in the models  
tf.activities_filtered <- tf.activities_matrix[-tf_cdkn1a_only,]

## get the same result, strings of target names
unlist(tf.list[["DMAP1"]])
tf.list$AR

## test the calculation of activities
target_name <- unlist(tf.list[["DMAP1"]])
colMeans(as.matrix(expr[target_name, ])) == tf.activities["DMAP1",]

## another method to solve colMean() bug =======================================
## add if else loop
# tf.activity_list <- lapply(tf.targets,
#                             function(value)
#                             {
#                               target_name <- unlist(value)
#                               if (length(target_name)>1){
#                                 colMeans(as.matrix(expr[target_name, ]))
#                               } else (expr[target_name, ])
#                             }
# )
# ==============================================================================

# @param expr           GTEX expression matrix across N samples
# @param tf.activities  TF activity matrix (N samples by K TFs)
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

# input expression data matrix (rownames: gene symbols)
expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
expr <- expr_tissue_median_gtex$data

tf.activities <- tf.activities_filtered
target <- "CDKN1A"
tf.candidates <- readRDS("data/tf.candidates.rds")

## == Fit a univariate linear model to infer regulatory effects =====================

ulinear_model_II <- function(expr, tf.activities, target){
  apply(tf.activities, MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
    coef(fit)
  })
}
ulinear_II_coef <- t(ulinear_model_II(expr, tf.activities, target))
qwrite(ulinear_II_coef, "results/ulinear_II_coef.rds");

## == Fit a multivariate linear model to infer regulatory effects ====================

# check before use: rownames(tf.activities) %in% rownames(expr)
mlinear_model_II <- function(expr, tf.activities, target){
  fit <- lm(expr[target, ] ~ t(tf.activities));
  coef(fit)
}
# call the mlinear function
mlinear_II_coef <- as.matrix(mlinear_model_II(expr, tf.activities, target))

# Remove 't(tf.activities)' from each row name
rownames(mlinear_II_coef) <- gsub("t\\(tf.activities\\)", "", rownames(mlinear_II_coef))
qwrite(mlinear_II_coef, "results/mlinear_II_coef.rds");

## CXXC1, EBF3 and ZBTB5 only have one target 'CDKN1A'
## If didn't use tf.activities_filtered [remove the 3 rows], coef will be NA
## CXXC1 appear first, EBF3 activity exactly the same as CXXC1 => co-linearity

mlinear_II_summary <- function(expr, tf.activities, target){
  fit <- lm(expr[target, ] ~ t(tf.activities));
  summary(fit)
}
mlinear_II_summary(expr, tf.activities, target)

## == Fit a lasso model to infer regulatory effects ==================================

lasso_model_II <- function(expr, tf.activities, target){
  library(glmnet)
  x <- t(tf.activities)
  y <- expr[target,]
  fit <- glmnet(x, y, alpha = 1)
  plot(fit)
  cv_model <- cv.glmnet(x, y, alpha = 1)
  plot(cv_model)
  coef(cv_model, s = "lambda.min")
}

## since CXXC1 activity is exactly same as CDKN1A expression
## must use tf.activities_filtered OR
## lasso just chose CXXC1 and gave a coefficient ~ 0.97

## save as matrix to convert all the dots to 0
lasso_II_coef <- as.matrix(lasso_model_II(expr, tf.activities, target))

# In a Lasso model, the coefficients of some variables can be represented as a dot . 
# when the L1 penalty parameter (lambda) is set high enough such that the coefficients 
# of these variables are shrunk to zero by the L1 penalty.

qwrite(lasso_II_coef, "results/lasso_II_coef.rds")

## When running the Lasso regression with the same data, 
## it is expected to obtain different coefficients at s = "lambda.min" each time.

## == Fit a mombf model to infer regulatory effects ==================================

mombf_model_II <- function(expr, tf.activities, target){
  library(mombf)
  x <- t(tf.activities)
  y <- expr[target,]
  msfit <- modelSelection(y, x, family="normal", priorCoef=momprior())
  sampl <- rnlp(msfit=msfit)
  colMeans(sampl)[-ncol(sampl)]
}

mombf_II_coef <- as.matrix(mombf_model_II(expr, tf.activities, target))
# Remove 'x' from each row name
rownames(mombf_II_coef) <- sub("x", "", rownames(mombf_II_coef))

saveRDS(mombf_II_coef, "results/mombf_II_coef.rds")
