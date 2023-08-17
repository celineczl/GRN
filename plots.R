## plot tf coefficients on y axis and TF names in the same order on x axis
## TODO: add line y=0 to the plot

gene_names <- rownames(ulinear_I_coef)
coefficients <- ulinear_I_coef[, "row"]
par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
plot(coefficients, xlab = "Gene Names", ylab = "coefficients",xaxt = "n", main = "class I univariate linear regression model")
axis(1,at=1:length(gene_names),labels=gene_names,las=2,cex.axis=0.6)

gene_names <- rownames(ulinear_II_coef)
coefficients <- ulinear_II_coef[, "row"]
par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
plot(coefficients, xlab = "Gene Names", ylab = "coefficients",xaxt = "n", main = "class II univariate linear regression model")
axis(1,at=1:length(gene_names),labels=gene_names,las=2,cex.axis=0.6)

## ============================================================

gene_names <- rownames(mlinear_I_coef)
par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
plot(mlinear_I_coef, xlab = "Gene Names", ylab = "coefficients",xaxt = "n", main = "class I multivariate linear regression model")
axis(1,at=1:length(gene_names),labels=gene_names,las=2,cex.axis=0.6)

gene_names <- rownames(mlinear_II_coef)
par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
plot(mlinear_II_coef, xlab = "Gene Names", ylab = "coefficients",xaxt = "n", main = "class II multivariate linear regression model")
axis(1,at=1:length(gene_names),labels=gene_names,las=2,cex.axis=0.6)

## ===========================================================

gene_names <- rownames(mombf_I_coef)
par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
plot(mombf_I_coef, xlab = "Gene Names", ylab = "coefficients",xaxt = "n", main = "class I mombf model")
axis(1,at=1:length(gene_names),labels=gene_names,las=2,cex.axis=0.6)

gene_names <- rownames(mombf_II_coef)
par(mar = c(5, 4, 4, 2) + 0.1,cex.axis=0.8)
plot(mombf_II_coef, xlab = "Gene Names", ylab = "coefficients",xaxt = "n", main = "class II mombf model")
axis(1,at=1:length(gene_names),labels=gene_names,las=2,cex.axis=0.6)

## =========================================================
## TODO: lasso
