## class I multivariate linear regression model to infer regulatory effect
## assume TF activity = expression level

library(io)

# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
expr <- expr_tissue_median_gtex$data

target <- "CDKN1A"
tf.candidates <- readRDS("data/tf.candidates.rds")
# check before use: tf.candidates %in% rownames(expr)

# fit <- lm(expr[target, ] ~ t(expr[tf.candidates,]))
# summary(fit)

mlinear_model_I <- function(expr, tf.candidates, target){
  fit <- lm(expr[target, ] ~ t(expr[tf.candidates,]));
  coef(fit)
}

mlinear_I_coef <- as.matrix(mlinear_model_I(expr, tf.candidates, target))
## from HES1, all coefficients turned to NA
## explain why: the number of parameters to be estimated is greater than the number of sample tissues.
## N*1 = N*D D*1, D is greater than N. more unknowns than equations, 
## only first 53 genes' coefficients and the intercept can be obtained the multivariate model.
## solutions: regularization/variable selection: lasso and non-local prior
## lasso: designed to do shrinkage, the variable selection function is not accurate

# Remove 't(expr[tf.candidates, ])' from each row name
rownames(mlinear_I_coef) <- sub("t\\(expr\\[tf\\.candidates, \\]\\)", "", rownames(mlinear_I_coef))

saveRDS(mlinear_I_coef,"results/mlinear_I_coef.rds")

## ===================== model testing =============================================
## NOT OUR FOCUS SO FAR

# testing data is stored in x_test and y_test variables
TF_test <- c("STAT3","CEBPA","CEBPB","VDR","RARA")
x_test <- as.data.frame(expr[TF_test,])
y_test <- expr["CDKN1A",]

predictions <- predict(model, newdata = x_test)

plot(y_test, predictions, xlab = "Actual Values", ylab = "Predicted Values", main = "Predicted vs Actual")

# normalized root mean squared error for regression model evaluation
nrmse <- sqrt(mean((predictions - y_test)^2)) / (max(y_test) - min(y_test))
summary(nrmse)
print(nrmse)

expr["STAT3",]
expr["STAT5",]
# no STAT5
expr["CEBPA",]
expr["CEBPB",]
expr["VDR",]
expr["RARA",]

residual <- predictions - y_test
hist(residual)
## check assumption, residual is normally distributed