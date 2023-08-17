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

tf.targets <- readRDS("tf.list.rds"); # 156 TFs and their target list

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
saveRDS(tf.activities, "tf.activities.rds")

# @param expr           GTEX expression matrix across N samples
# @param tf.activities  TF activity matrix (N samples by K TFs)
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

# input expression data matrix (rownames: gene symbols)
expr <- readRDS("expr.rds")
tf.activities <- readRDS("tf.activities.rds")
target <- "CDKN1A"
tf.candidates <- readRDS("tf.candidates.rds")

## == Fit a univariate linear model to infer regulatory effects =====================

ulinear_model_II <- function(expr, tf.activities, target){
  apply(tf.activities, MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
    coef(fit)
  })
}
ulinear_II_coef <- t(ulinear_model_II(expr, tf.activities, target))
qwrite(ulinear_II_coef, "ulinear_II_coef.rds");

fit <- lm(expr[target,] ~ tf.activities["ZBTB5",]);
coef(fit) # NA
tf.activities["ZBTB5",]
# TODO: NO missing values in data, why NA?

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
qwrite(mlinear_II_coef, "mlinear_II_coef.rds");

## == Fit a lasso model to infer regulatory effects ==================================

lasso_model_II <- function(expr, tf.activities, target){
  library(glmnet)
  x <- t(tf.activities)
  y <- expr[target,]
  fit <- glmnet(x, y, alpha = 1)
  coef(fit)
}
lasso_model_II(expr, tf.activities, target)

lasso_II_coef <- lasso_model_II(expr, tf.activities, target)
# In a Lasso model, the coefficients of some variables can be represented as a dot . 
# when the L1 penalty parameter (lambda) is set high enough such that the coefficients 
# of these variables are shrunk to zero by the L1 penalty.
qwrite(lasso_II_coef, "lasso_II_coef.rds");

## == Fit a mombf model to infer regulatory effects ==================================

mombf_model_II <- function(expr, tf.activities, target){
  library(mombf)
  x <- t(tf.activities)
  y <- expr[target,]
  fit <-bestBIC(y ~ x)
  coef(fit)
}

mombf_II_coef <- as.matrix(mombf_model_II(expr, tf.activities_filtered, target))
# Remove 'x' from each row name
rownames(mombf_II_coef) <- sub("x", "", rownames(mombf_II_coef))

saveRDS(mombf_II_coef, "mombf_II_coef.rds")

## mombf_model_II function debug:
# mombf_model_II(expr, tf.activities, target)
# Error in modelSelection(..., priorCoef = bic(), priorDelta = modelunifprior(), :
# There are >1 constant columns in x (e.g. two intercepts)

x <- t(tf.activities)
y <- expr[target,]

# Check for constant columns
constant_cols <- apply(x, 2, function(col) var(col) == 0)
# Get the indices of constant columns
constant_col_indices <- which(constant_cols) # 14 23 144

x1 <- x[,-constant_col_indices]
fit <-bestBIC(y ~ x1)

tf.activities_filtered <- t(x1)
# coincidence or poor data?



plot(ulinear_II_coef)
plot(mlinear_coef)
plot(lasso_coef)

# Combine the coefficients into a data frame
coefficients <- data.frame(Model = c("ulinear", "mlinear", "lasso"),
                           Coefficient = c(ulinear_coef["x1"], mlinear_coef["x2"], lasso_coef["x3"]))

# Create a bar plot
library(ggplot2)

# TODO: plot tf coefficients on y axis and TF names in the same order on x axis.

plot <- ggplot(coefficients, aes(x = Model, y = Coefficient)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Model", y = "Coefficient") +
  ggtitle("Comparison of Coefficients") +
  theme_minimal()

print(plot)