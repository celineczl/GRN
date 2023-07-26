setwd("D:/GRN")
expr_tissue_median_gtex <- readRDS("C:/Users/Celin/Downloads/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)
expr <- expr_tissue_median_gtex$data
TF <- c("TP53","TFAP4","E2F1","E2F3","SP1","SP3","TCF3","TFAP2A","TFAP2C","TFAP2E","STAT1")

.libPaths("D:/R-4.1.2/library")
# .libPaths()
# install.packages("glmnet",lib="D:/R-4.1.2/library")
# need to specify library path everytime, defualt usually not work
library(glmnet,lib="D:/R-4.1.2/library")
# glmnet require latest R version 4.3.1 and Rtools for install

## ======== LASSO variable selection model ====================
## https://glmnet.stanford.edu/articles/glmnet.html#quick-start

# The default model used in the package is the Guassian linear model or “least squares”
# input matrix x is TF expression data
x <- t(expr[TF,])
# response variable is the target gene expression
y <- expr["CDKN1A",]

# data distribution
hist(expr)
hist(expr[TF,])

# Fit Lasso model
lasso_model <- glmnet(x, y, alpha = 1)
# default, Gaussian linear regression
# alpha = 1, lasso
# alpha = 0, ridge

# Perform cross-validation
cv_model <- cv.glmnet(x, y, alpha = 1)

plot(lasso_model)
print(lasso_model)
plot(cv_model)

# min lambda, max model complexity
cv_model$lambda.min

coef(cv_model, s = "lambda.min")
predict(cv_model, newx = x[1:5,], s = "lambda.min")

expr["CUX",]

coef(cv_model, s = "lambda.1se")

fit <- glmnet(x, y)
# the fitted values for the observations at 𝜆=0.05
predict(fit, newx = x, type = "response", s = 0.05)

# log lambda
plot(fit, xvar = "lambda", label = TRUE)
# fraction deviance explained
plot(fit, xvar = "dev", label = TRUE)

## ============== lasso inferetial model ====================
# choose 
  


# TRUNCATED normal distribution, truncated at 0


