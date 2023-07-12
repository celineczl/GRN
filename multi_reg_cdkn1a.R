expr_tissue_median_gtex <- readRDS("C:/Users/Celin/Downloads/expr_tissue-median_gtex.rds")
head(expr_tissue_median_gtex)
summary(expr_tissue_median_gtex)
expr <- expr_tissue_median_gtex$data

## ========= multivariate linear regression ================
model <- lm(expr["CDKN1A",] ~ expr["TP53",]+expr["TFAP4",]+expr["E2F1",]+expr["E2F3",]+expr["SP1",]
            +expr["SP3",]+expr["TCF3",]+expr["TFAP2A",]+expr["TFAP2C",]+expr["TFAP2E",]+expr["STAT1",])
summary(model)

TF <- c("TP53","TFAP4","E2F1","E2F3","SP1","SP3","TCF3","TFAP2A","TFAP2C","TFAP2E","STAT1")

# to be updated with new TFs...

# testing data is stored in x_test and y_test variables
TF_test <- c("STAT3","STAT5","CEBPA","CEBPB","VDR")
x_test <- 

predictions <- predict(model, newdata = TF_test)

plot(y_test, predictions, xlab = "Actual Values", ylab = "Predicted Values", main = "Predicted vs Actual")

nrmse <- sqrt(mean((predictions - y_test)^2)) / (max(y_test) - min(y_test))

expr["STAT3",]
expr["STAT5",]
expr["CEBPA",]
expr["CEBPB",]
expr["VDR",]
