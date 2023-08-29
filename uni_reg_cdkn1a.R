## class I univariate linear regression model to infer regulatory effect
## assume TF activity = expression level

library(io)

# GTEX median expression of each gene across different tissues:
expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)

expr <- expr_tissue_median_gtex$data

# TF <- read.csv("data/TF.csv", header = FALSE)

## == Fit a univariate linear model to infer regulatory effects =====================
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

target <- 'CDKN1A'
tf.candidates <- readRDS("data/tf.candidates.rds")

ulinear_model_I <- function(expr, tf.candidates, target){
  apply(expr[tf.candidates,], MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
    coef(fit)
  })
}
ulinear_I_coef <- t(ulinear_model_I(expr, tf.candidates, target))
saveRDS(ulinear_I_coef,"results/ulinear_I_coef.rds")

ulinear_model_plot <- function(expr, tf.candidates, target, output_file){
  
  # Create a PDF device to save the plots
  pdf(output_file)
  
  # scatter plot of target expression vs each TF expression
  for (candidate in tf.candidates) {
    plot(expr[candidate, ], expr[target, ],
      xlab = paste0(target, " expression"),
      ylab = paste0(candidate, " expression")
    )
    
    # add ulinear model line [y = intercept + coefficient * x]
    abline(
      a = ulinear_I_coef[candidate, 1],
      b = ulinear_I_coef[candidate, 2]
    )
  }
  
  # Save the current plot to the PDF file
  dev.off()
}

# Specify the name of the output PDF file
output_file <- "results/ulinear_output_plots.pdf"

ulinear_model_plot(expr, tf.candidates, target, output_file)
