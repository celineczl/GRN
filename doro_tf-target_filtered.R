# install.packages("readxl")
library(readxl)

# Table S3. Scored TF-target interactions for normal tissues. Columns indicate: 
# 1) TF, 2) Target, 3) Effect, 4) Confidence score; 5) logical if from curated database, 
# 6) logical if from ChIPseq, 7) logical if from TF binding motif predictions 8) logical if inferred from GTEx expression data, 
# 9) curated database name, 10) ChIPseq resource name, 10) TF binding motif name, 11) inference strategy name, 
# 12) pubmedID -from curated resources-, and 13) kegg pathway id. 
# For the regulons inferred from GTEx we only include the consensus ones (present in more than 2 tissues), 
# the rest are available at https://github.com/saezlab/TFbenchmark.

## truncated table from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/
## dorothea_data <- read_excel("data/GarciaAlonso_supplemental_table_S3_regulonsNormal.xlsx", sheet = "Sheet 1 - Table S3. Scored TF-t")

## =====================================================================================================================================

# GTEX median expression of each gene across different tissues:
expr_tissue_median_gtex <- readRDS("data/expr_tissue-median_gtex.rds")
# expression matrix
expr <- expr_tissue_median_gtex$data

# TF list of CDKN1A from literature_curated table
TF <- read.csv("data/TF.csv", header = FALSE)

## revised table from https://genome.cshlp.org/content/29/8/1363/suppl/DC2
dorothea_data <- read.csv("data/Revised_Supplemental_Table_S3.csv", header = TRUE)

# remove all the TF-target interaction inferred from co-expression
filtered_doro <- subset(dorothea_data, is_evidence_coexpression == FALSE)

# obtain all the TF-target entries (target = CDKN1A in this case)
doro_cdkn1a <- subset(filtered_doro, target == "CDKN1A") # 144 obs

## haven't integrate the doro_cdkn1a with literature_curated table yet

# exclude TF not included in expr
in_target <- doro_cdkn1a$TF %in% rownames(expr) # 144 obs
tf_cdkn1a <- doro_cdkn1a[in_target,]

# retrieve all the targets of all the TFs of CDKN1A
tf.matrix <- subset(filtered_doro, TF %in% tf_cdkn1a$TF) # 249289 obs

# remove targets not included in GTEX data
tf.matrix_gtex <- subset(tf.matrix, target %in% rownames(expr))
                    
# split the matrix to a list
tf.list <- split(tf.matrix_gtex[, 2], tf.matrix_gtex[, 1])




