## Function Documentation

#===== class I models to infer regulatory effects (coef) [assume TF activity = expression level] ====================
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

# class I univariate linear regression model
# output the coefficients of variables
ulinear_model_I <- function(expr, tf.candidates, target){
  apply(expr[tf.candidates,], MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
    coef(fit) 
  })
}

# class I multivariate linear regression model
mlinear_model_I <- function(expr, tf.candidates, target){
  fit <- lm(expr[target, ] ~ t(expr[tf.candidates,]));
  coef(fit)
}

# class I lasso variable selection model
# TODO: function to be completed

# class I mombf variable selection model
# TODO: function to be completed

#===== class I models to infer regulatory effects [assume TF activity = expression level] ====================
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs

# class II univariate linear regression model
# output the coefficients of variables
ulinear_model_II <- function(expr, tf.activities, target){
  apply(tf.activities, MARGIN = 1, function(row) {
    fit <- lm(expr[target,] ~ row);
    coef(fit)
  })
}

# class II multivariate linear regression model
mlinear_model_II <- function(expr, tf.activities, target){
  fit <- lm(expr[target, ] ~ t(tf.activities));
  coef(fit)
}

# class II lasso variable selection model
lasso_model_II <- function(expr, tf.activities, target){
  library(glmnet)
  x <- t(tf.activities)
  y <- expr[target,]
  fit <- glmnet(x, y, alpha = 1)
  plot(fit)
  cv_model <- cv.glmnet(x, y, alpha = 1)
  plot(cv_model)
  coef(cv_model, s = "lambda.min")
}

# class II mombf variable selection model
mombf_model_II <- function(expr, tf.activities, target){
  library(mombf)
  x <- t(tf.activities)
  y <- expr[target,]
  msfit <- modelSelection(y, x, family="normal", priorCoef=momprior())
  sampl <- rnlp(msfit=msfit)
  colMeans(sampl)[-ncol(sampl)]
}

#========= functions for plotting and summary of coefficients from models ==================================
# @param expr           GTEX expression matrix across N samples
# @param target         gene name of target
# @param tf.candidates  gene names of candidate TFs
# @param output_file    file path of output pdf file "results/ulinear_output_plots.pdf"
# @param ulinear_coef   coefficient matrix obtained from univariate linear regression models

# plot scatter plots for class I univariate linear regression to visualize
# target expression vs TF expression and output to a pdf file
ulinear_model_plot <- function(expr, tf.candidates, target, ulinear_coef, output_file){
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
      a = ulinear_coef[candidate, 1], 
      # intercept is the first column of univariate model coef matrix
      b = ulinear_coef[candidate, 2]  
      # gradient = coef is the second column of univariate model coef matrix
    )
  }
  # Save the current plot to the PDF file
  dev.off()
}

# get the summary of class II multivariate linear regression model
mlinear_II_summary <- function(expr, tf.activities, target){
  fit <- lm(expr[target, ] ~ t(tf.activities));
  summary(fit)
}

# ================== 
# Infer TF activity matrix (all samples by all TFs)
# target is from the list
tf.activity_list <- lapply(tf.targets,
                           function(value)
                           {
                             target_name <- unlist(value)
                             # set drop dimension = FALSE, 
                             # or colMeans() will average the row expr of CDKN1A when no other targets
                             colMeans(as.matrix(expr[target_name, , drop=FALSE]))
                           }
)