setwd("D:/GRN")


expr_tissue_median_gtex <- readRDS("C:/Users/Celin/Downloads/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)
expr <- expr_tissue_median_gtex$data

row.names(expr)

# Get the row index where the row name matches the search word
row_index <- which(row.names(expr) == "CDKN1A")

# Print the row index
print(row_index)
# 17656
CDKN1A <- expr[17656,]
# double check the row name
row.names(expr)[17656]

TF <- c("TP53","TFAP4","E2F1","E2F3","SP1","SP3","TCF3","TFAP2A","TFAP2C","TFAP2E","STAT1")
# row_index <- which(row.names(expr) == TF)
print(row_index)

row_index <- which(row.names(expr) == "TP53")
print(row_index) # 42272

# set multiple uni-variate linear regression across tissues
model <- lm(expr[17656,] ~ expr[42272,])
summary(model)

row_index <- which(row.names(expr) == "AP1")
print(row_index)
# no AP1

row_index <- which(row.names(expr) == "TFAP2A")
print(row_index) #16838
model <- lm(expr[17656,] ~ expr[16838,])
summary(model)

model <- lm(expr[17656,] ~ expr[16838,]+expr[42272,])
summary(model)

row_index <- which(row.names(expr) == "TFAP2C")
print(row_index) #49604
model <- lm(expr[17656,] ~ expr[49604,])
summary(model)
# smaller p-value than TFAP2A
# focus on coeficient but not p-values

# TFAP2 A B C D E
# TF2P2A and TFAP2C and TFAP2E important to CDKN1A

## ========== UNIVARIATE REGRESSION ==============
expr_tissue_median_gtex <- readRDS("C:/Users/Celin/Downloads/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)
expr <- expr_tissue_median_gtex$data

TF <- c("TP53","TFAP4","E2F1","E2F3","SP1","SP3","TCF3","TFAP2A","TFAP2C","TFAP2E","STAT1")
# can directly call row name

uni_regr_CDKN1A_TF = list()

for (i in TF){
  tf_exp <- expr[i,]
  print(as.matrix(tf_exp))
  
  uni_regr_CDKN1A_TF[[i]] = lm(expr["CDKN1A",] ~ tf_exp)
  summary(uni_regr_CDKN1A_TF[[i]])
  ## [[i]] are used to directly access the value of a specific element within a list
  
  pdf(file = paste0("D:/R-Workshop/output1_summary/",i,".pdf"),width=7,height=7)
  print(summary(uni_regr_CDKN1A_TF[[i]]))
  print(plot(uni_regr_CDKN1A_TF[[i]]))
  ## [i] are used to subset and return a sub-list or sub-data frame
  dev.off()
}

# Generate model summary and capture output
model_summary <- capture.output(summary(uni_regr_CDKN1A_TF[[i]]))

# Print model summary within the PDF
cat(model_summary, sep = "\n")
