## Code by djhshih. Title: cshkg. 
## Retrieved from https://github.com/djhshih/cshkg/blob/main/pancan/data/preprocess.R. 
## Last updated: March 1, 2023.

## Data source: RNA (Final) - EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv
## Retrieved from https://gdc.cancer.gov/about-data/publications/pancanatlas [Supplemental Data]
## The Cancer Genome Atlas Pan-Cancer analysis project [https://doi-org.eproxy.lib.hku.hk/10.1038/ng.2764]

library(io)
library(data.table)
library(org.Hs.eg.db)

# human genome build in PanCanAtlas: hg19

in.fn <- as.filename("tcga-pancan-expr.tsv");
out.fn <- set_fext(in.fn, "rds"); # change file extension from tsv to rds

expr0 <- fread(tag(in.fn)); # data.frame

features <- expr0[[1]]; # $gene_id
features.ss <- strsplit(features, "|", fixed=TRUE);
# The `features` vector is split by the "|" character using the `strsplit` function 

# do not use gene symbols from original data:
# genes <- unlist(lapply(features.ss, function(x) x[1]));
# genes[genes == "?"] <- NA; 

entrezs <- unlist(lapply(features.ss, function(x) x[2]));
# second element of each split string is extracted into a new vector called `entrezs`

# This is done because the gene IDs in the input file are in the format "gene_symbol|entrez_id", 
# and we want to extract only the entrez IDs for further analysis.

expr <- expr0;
expr[, gene_id:=NULL]; # remove gene_id column
expr <- as.matrix(expr);

sum(duplicated(entrezs)) # check whether entrez IDs are unique

# re-obtain gene symbols from database
syms.d <- select(org.Hs.eg.db, keys=entrezs, columns=c("SYMBOL"),keytype="ENTREZID");
genes <- syms.d[match(entrezs, syms.d[,1]), 2];

idx <- duplicated(genes);
print(genes[idx])

valid <- !idx & !is.na(genes); # not duplicated AND is not `NA`
expr <- expr[valid, ];
genes <- genes[valid];
entrezs <- entrezs[valid];
rownames(expr) <- genes;

min.expr <- min(as.numeric(expr), na.rm=TRUE); # negative expression

annot <- data.frame(gene = genes, entrez_id = entrezs);

# log transform
# re-scale by minus min.expr
eset <- list(
  features = annot,
  data = log2(expr - min.expr + 1)
);

min(as.numeric(eset$data), na.rm=TRUE) # 0

qwrite(eset, out.fn); # save as tcga-pancan-expr.rds in same directory
