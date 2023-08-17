## class I univariate linear regression model to infer regulatory effect
## assume TF activity = expression level

library(io)

# GTEX median expression of each gene across different tissues:
expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)

expr <- expr_tissue_median_gtex$data
# qwrite(expr, "expr.rds");

TF <- read.csv("data/TF.csv", header = FALSE)

## == Fit a univariate linear model to infer regulatory effects =====================
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

target <- 'CDKN1A'
# tf.candidates <- tf_cdkn1a$TF
# saveRDS(tf.candidates, "tf.candidates.rds")
tf.candidates <- readRDS("tf.candidates.rds")

ulinear_model_I <- function(expr, tf.candidates, target){
  apply(expr[tf.candidates,], MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
    coef(fit)
  })
}
ulinear_I_coef <- t(ulinear_model_I(expr, tf.candidates, target))
saveRDS(ulinear_I_coef,"ulinear_I_coef.rds")

## scatter plot of target expression vs each TF expression
## TODO: y-axis label is always "row"
ulinear_model_plot <- function(expr, tf.candidates, target){
  apply(expr[tf.candidates,], MARGIN = 1, function(row){
    plot(expr[target,],row)
  })
}
ulinear_model_plot(expr, tf.candidates, target)

