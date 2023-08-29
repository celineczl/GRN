## retrieve TF-target list from DoRothEA raw data table
## remove entries only inferred from co-expression

# install.packages("readxl")
library(readxl)
library(dplyr)

# GTEX median expression of each gene across different tissues:
expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
# expression matrix
expr <- expr_tissue_median_gtex$data

# TF list of CDKN1A from literature_curated table
TF <- read.csv("data/TF.csv", header = FALSE)

## revised table from https://genome.cshlp.org/content/29/8/1363/suppl/DC2
dorothea_data <- read.csv("data/Revised_Supplemental_Table_S3.csv", header = TRUE)

# filter out the TF-target interaction ONLY inferred from co-expression
coexpr_only <- subset(dorothea_data,
  is_evidence_curateddatabase == FALSE & 
  is_evidence_chipSeq == FALSE & 
  is_evidence_TFbindingMotif == FALSE & 
  is_evidence_coexpression == TRUE
)

# remove 'coexpr_only' from 'dorothea_data'
filtered_doro <- anti_join(dorothea_data, coexpr_only)
qwrite(filtered_doro, "data/filtered_doro.rds");

# obtain all the TF-target entries (target = CDKN1A in this case)
doro_cdkn1a <- subset(filtered_doro, target == "CDKN1A") # 156 obs

## haven't integrate the doro_cdkn1a with literature_curated table yet

# exclude TF not included in expr
in_target <- doro_cdkn1a$TF %in% rownames(expr) # 156 obs
tf_cdkn1a <- doro_cdkn1a[in_target,]
tf.candidates <- tf_cdkn1a$TF
saveRDS(tf.candidates, "data/tf.candidates.rds")

# retrieve all the targets of all the TFs of CDKN1A
tf.matrix <- subset(filtered_doro, TF %in% tf_cdkn1a$TF) # 270731 obs

# remove targets not included in GTEX data
tf.matrix_gtex <- subset(tf.matrix, target %in% rownames(expr))
                    
# split the matrix to a list
tf.list <- split(tf.matrix_gtex[, 2], tf.matrix_gtex[, 1])
saveRDS(tf.list, "data/tf.list.rds")
