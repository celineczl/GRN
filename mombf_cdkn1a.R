## class I mombf (non-local prior variable selection) model to infer regulatory effect
## assume TF activity = expression level

## ======== mombf variable selection model ====================
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

expr <- readRDS("expr.rds")
target <- "CDKN1A"
tf.candidates <- readRDS("tf.candidates.rds")

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
fit <-bestBIC(y ~ x)
summary(fit) # usual GLM summary
coef(fit)

mombf_I_coef <- as.matrix(coef(fit)) # MLE under top model
## TODO: NA coefficient started from HES5

# Remove 'x' from each row name
rownames(mombf_I_coef) <- sub("x", "", rownames(mombf_I_coef))

saveRDS(mombf_I_coef, "mombf_I_coef.rds")

confint(fit) # conf int under top model

selected_models <- fit$topmodel.fit
# best model selected, no need to perform mombf?
