## class I multivariate linear regression model to infer regulatory effect
## assume TF activity = expression level

library(io)

# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

expr <- readRDS("expr.rds")
tf.candidates <- readRDS("tf.candidates.rds")
# check before use: tf.candidates %in% rownames(expr)

# fit <- lm(expr[target, ] ~ t(expr[tf.candidates,]))
# summary(fit)

mlinear_model_I <- function(expr, tf.candidates, target){
  fit <- lm(expr[target, ] ~ t(expr[tf.candidates,]));
  coef(fit)
}

mlinear_I_coef <- as.matrix(mlinear_model_I(expr, tf.candidates, target))
## from HES1, all coefficients turned to NA
## possibly reasons:
## multicollinearity
## missing values [but no NA in dataset]
## insufficient variability

# TODO: examine the correlation matrix or calculate the variance inflation factor (VIF) for each predictor

# Remove 't(expr[tf.candidates, ])' from each row name
rownames(mlinear_I_coef) <- sub("t\\(expr\\[tf\\.candidates, \\]\\)", "", rownames(mlinear_I_coef))

saveRDS(mlinear_I_coef,"mlinear_I_coef.rds")

## ================== testing ==============================================

# remove the 152:155
tf.candidates1 <- setdiff(tf.candidates, c("ZNF263", "ZNF300", "ZNF331", "ZNF76"))
# tf.candidates1 <- tf.candidates[c(1:151,156)]

fit <- lm(expr[target, ] ~ t(expr[tf.candidates1,]));
summary(fit)

tf.candidates2 <- tf.candidates[1:50] # GOOD
tf.candidates3 <- tf.candidates[51:100] # GOOD
# remove term 152:155
tf.candidates4 <- tf.candidates[c(101:151,156)] # GOOD

# checked the data matrix, no NA
expr1 <- expr[tf.candidates[101:156],]

fit <- lm(expr[target, ] ~ t(expr[tf.candidates1,]));
summary(fit)
## still NA

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