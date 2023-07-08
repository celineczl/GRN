setwd("D:/GRN")
expr_tissue_median_gtex <- readRDS("C:/Users/Celin/Downloads/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)
expr <- expr_tissue_median_gtex$data
TF <- c("TP53","TFAP4","E2F1","E2F3","SP1","SP3","TCF3","TFAP2A","TFAP2C","TFAP2E","STAT1")

## ======== mombf variable selection model ====================

## sessionInfo()
# install.packages("rlang")
library(rlang) # rlang_1.1.1

# install.packages("mombf",lib="D:/R-4.1.2/library")
library(mombf,lib="D:/R-4.1.2/library")

# y (target gene) and x (transcription factors)
x <- t(expr[TF,])
y <- expr["CDKN1A",]

fit= bestBIC(y ~ x)
summary(fit) # usual GLM summary
coef(fit) # MLE under top model
confint(fit) # conf int under top model

selected_models <- fit$topmodel.fit
# best model selected, no need to perform mombf?
