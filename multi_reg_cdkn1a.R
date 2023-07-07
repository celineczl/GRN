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
