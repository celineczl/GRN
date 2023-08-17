## class I lasso model to infer regulatory effect
## assume TF activity = expression level

library(glmnet)
# glmnet require latest R version 4.3.1 and Rtools for install

## ======== LASSO variable selection model ====================
## https://glmnet.stanford.edu/articles/glmnet.html#quick-start
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

expr <- readRDS("expr.rds")
target <- "CDKN1A"
tf.candidates <- readRDS("tf.candidates.rds")

# The default model used in the package is the Guassian linear model or “least squares”
# input matrix x is TF expression data
x <- t(expr[tf.candidates,])
# response variable is the target gene expression
y <- expr[target,]

# data distribution
hist(expr)
hist(expr[tf.candidates,])

# Fit Lasso model
lasso_model <- glmnet(x, y, alpha = 1)
# default, Gaussian linear regression
# alpha = 1, lasso
# alpha = 0, ridge

plot(lasso_model)

# summary information for the Lasso model:
print(lasso_model)

# store lasso coefficients, format: dgCMatrix
lasso_I_coef <- coef(lasso_model)
print(lasso_I_coef)

##============== Perform cross-validation ==========================

cv_model <- cv.glmnet(x, y, alpha = 1)

plot(cv_model)

# min lambda, max model complexity
cv_model$lambda.min

lasso_I_coef <- as.matrix(coef(cv_model, s = "lambda.min"))


saveRDS(lasso_I_coef, "lasso_I_coef.rds")

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


