## class I univariate linear regression model
## assume TF activity = expression level

# GTEX median expression of each gene across different tissues:
expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)

expr <- expr_tissue_median_gtex$data
TF <- read.csv("data/TF.csv", header = FALSE)

## ========== CLASS I UNIVARIATE LINEAR REGRESSION ==============

uni_regr_CDKN1A_TF <- list()
for (i in TF){
  tf_exp <- expr[i,]
  print(as.matrix(tf_exp))
  # for each TF in the list, perform univariate linear regression
  # output into pdf
  
  uni_regr_CDKN1A_TF[[i]] <- lm(expr["CDKN1A",] ~ tf_exp)
  summary(uni_regr_CDKN1A_TF[[i]])
  ## [[i]] are used to directly access the value of a specific element within a list
  
  pdf(file = paste0("D:/R-Workshop/output1_summary/",i,".pdf"),width=7,height=7)
  print(summary(uni_regr_CDKN1A_TF[[i]]))
  print(plot(uni_regr_CDKN1A_TF[[i]]))
  ## [i] are used to subset and return a sub-list or sub-data frame
  dev.off()
}

# Generate model summary and capture output
# model_summary <- capture.output(summary(uni_regr_CDKN1A_TF[[i]]))

# Print model summary within the PDF
# cat(model_summary, sep = "\n")
