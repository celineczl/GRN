## plot barplots with tf coefficients on y axis and TF names in the same order on x axis

# @param coef           coefficient values for each TF in certain model
# @param title          string specifying title of the plot, e.g., "Class I univariate linear regression model"

model_coef_plot <-function(coef,title){
  gene_names <- rownames(coef)
  par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
  barplot(coef, xlab = "Gene Names", ylab = "Coefficients", 
  names.arg = gene_names, las = 2, cex.names = 0.6, 
  main = title)
}

## ============== univariate linear regression model coefficient plot ==================

coef <- ulinear_I_coef[, "row"]
title <- "Class I univariate linear regression model"
model_coef_plot(coef,title)

coef <- ulinear_II_coef[, "row"]
title <- "Class II univariate linear regression model"
model_coef_plot(coef,title)

## ============== multivariate linear regression model coefficient plot ==================

## remove the intercept in the first row
coef <- mlinear_I_coef[-1,1]
title <- "Class I multivariate linear regression model"
model_coef_plot(coef,title)

coef <- mlinear_II_coef[-1,1]
title <- "Class II multivariate linear regression model"
model_coef_plot(coef,title)

## ============== lasso model coefficient plot ==================

## remove the intercept in the first row
coef <- lasso_I_coef[-1,1]
title <- "Class I lasso model"
model_coef_plot(coef,title)

coef <- lasso_II_coef[-1,1]
title <- "Class II lasso model"
model_coef_plot(coef,title)

## ============== mombf model coefficient plot ==================

coef <- mombf_I_coef[-1,1]
title <- "Class I mombf model"
model_coef_plot(coef,title)

coef <- mombf_II_coef[-1,1]
title <- "Class II mombf model"
model_coef_plot(coef,title)

plot(rowMeans(expr[tf.candidates, ]),rowMeans(tf.activities_matrix),
     xlab = "TF expression", ylab = "TF activities", 
     xlim = c(0, max(rowMeans(expr[tf.candidates, ]))),
     ylim = c(0, max(rowMeans(expr[tf.candidates, ]))),
     cex = 0.8)
abline(a=0, b=1, col = "red")

qqplot(rowMeans(expr[tf.candidates, ]),rowMeans(tf.activities_matrix),
       xlab = "TF expression", ylab = "TF activities")

library(io)
qdraw(
  {
    plot(rowMeans(expr[tf.candidates, ]),rowMeans(tf.activities_matrix),
         xlab = "TF expression", ylab = "TF activities", 
         xlim = c(0, max(rowMeans(expr[tf.candidates, ]))),
         ylim = c(0, max(rowMeans(expr[tf.candidates, ]))))
  },
  width = 3, height = 3, res = 1000,
  file = "expr-activity.jpg"
)
