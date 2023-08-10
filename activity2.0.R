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

tf.targets <- tf.list; # 35 TFs and their target list

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
tf.candidates <- tf_cdkn1a$TF

## == Fit a univariate linear model to infer regulatory effects =====================

ulinear_model_II <- function(expr, tf.activities, target){
  apply(tf.activities, MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
    coef(fit)
  })
}
ulinear_coef <- t(ulinear_model_II(expr, tf.activities, target))
qwrite(ulinear_coef, "ulinear_coef.rds");

## == Fit a multivariate linear model to infer regulatory effects ====================

# check before use: rownames(tf.activities) %in% rownames(expr)
mlinear_model_II <- function(expr, tf.activities, target){
  fit <- lm(expr[target, ] ~ t(tf.activities));
  coef(fit)
}
# call the mlinear function
mlinear_coef <- as.matrix(mlinear_model_II(expr, tf.activities, target))

# Remove 't(tf.activities)' from each row name
rownames(mlinear_coef) <- gsub("t\\(tf.activities\\)", "", rownames(mlinear_coef))
qwrite(mlinear_coef, "mlinear_coef.rds");

## == Fit a lasso model to infer regulatory effects ==================================

lasso_model_II <- function(expr, tf.activities, target){
  library(glmnet)
  x <- t(tf.activities)
  y <- expr[target,]
  fit <- glmnet(x, y, alpha = 1)
  coef(fit)
}
lasso_model_II(expr, tf.activities, target)
lasso_coef <- as.matrix(lasso_model_II(expr, tf.activities, target))
# In a Lasso model, the coefficients of some variables can be represented as a dot . 
# when the L1 penalty parameter (lambda) is set high enough such that the coefficients 
# of these variables are shrunk to zero by the L1 penalty.
qwrite(lasso_coef, "lasso_coef.rds");

## == Fit a mombf model to infer regulatory effects ==================================

mombf_model_II <- function(expr, tf.activities, target){
  library(mombf)
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

plot(ulinear_coef)
plot(mlinear_coef)
plot(lasso_coef)

# Combine the coefficients into a data frame
coefficients <- data.frame(Model = c("ulinear", "mlinear", "lasso"),
                           Coefficient = c(ulinear_coef["x1"], mlinear_coef["x2"], lasso_coef["x3"]))

# Create a bar plot
library(ggplot2)

plot <- ggplot(coefficients, aes(x = Model, y = Coefficient)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Model", y = "Coefficient") +
  ggtitle("Comparison of Coefficients") +
  theme_minimal()

print(plot)