## retrieve TF-target list from database
## compare different databases
## 1) OmnipathR
## 2) DoRothEA and dorothea_hs
## 3) collecTRI

## To install the package from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dorothea")

## install.packages("devtools")
## To install the development version from the Github repo:
# devtools::install_github("saezlab/dorothea")

# install.packages("remotes")
# library(remotes)
# remotes::install_github('saezlab/OmnipathR')
# remotes::install_github('saezlab/decoupleR')

library(dorothea)
library(OmnipathR)
library(decoupleR)
detach("package:stats")
library(tidyverse) # detach("package:stats") or conflict with tidyverse


##=========================== OmnipathR ===============================================
## protein-protein interaction database
# OmnipathR_interaction_res: all the resources in OmnipathR
OmnipathR_interaction_res <- get_interaction_resources(dataset = NULL) # 142 resources entry

# no matter what datasets specified, under 'interactions' query_type, output the same 38 resources
interaction_res <- get_resources(query_type = 'interactions', datasets = 'dorothea', generic_categories = NULL)
get_resources(query_type = 'interactions', datasets = 'tf_target', generic_categories = NULL)
get_resources(query_type = 'interactions', datasets = 'collectri', generic_categories = NULL)

# raw regulons
OmnipathR::collectri(organism=9606L, genesymbols=TRUE, loops=TRUE)
# 9606 human (default), 10116 rat and 10090 Mouse
# https://www.bioconductor.org/packages/devel/bioc/manuals/OmnipathR/man/OmnipathR.pdf

TRI <- collectri(
  resources = 'CollecTRI',
  organism = 9606L,
  references_by_resource = TRUE,
  exclude = NULL,
  strict_evidences = TRUE
)

DORO <- dorothea(
  resources = 'DoRothEA',
  organism = 9606,
  dorothea_levels = c("A", "B", "C"),
  fields = NULL,
  default_fields = TRUE,
  references_by_resource = TRUE,
  exclude = NULL,
  strict_evidences = TRUE
)


##===================== processed regulons from collectri ==========================================
# reference: https://www.biorxiv.org/content/10.1101/2023.03.30.534849v1
# reference: collectri firstly introduced in 2021, ExTRI https://pubmed.ncbi.nlm.nih.gov/34875418/

net_TRI <- decoupleR::get_collectri(organism='human', split_complexes=FALSE) # 43178 obs
head(net_TRI) # no confidence level column like dorothea

# DoRothEA was downloaded using the function 'get_dorothea' from the decoupler 2.4.0 bioconductor package
#   and filtered by its confidence levels A, B and C
net_TRI_cdkn1a <- subset(net_TRI, target == "CDKN1A") # 254 obs


##===================== 'dorothea_hs' database in dorothea package: human TF-target interactions ================
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/
# 1395 TFs targeting 20,244 genes with 486,676 unique interactions
# http://bioconductor.org/packages/release/data/experiment/manuals/dorothea/man/dorothea.pdf

library(ggplot2)
library(dplyr)

net_DORO_hs <- subset(dorothea_hs, confidence %in% c("A", "B", "C", "D")) # 29086 obs

# visualize edge confidence level of target CDKN1A 
net_DORO_hs_cdkn1a_all <- subset(net_DORO_hs, target == "CDKN1A")
n_edges <- net_DORO_cdkn1a_all %>%
  group_by(confidence) %>%
  summarize(n = n())
ggplot(data=n_edges, aes(x=confidence, y=log10(n), color=confidence, fill=confidence)) +
  geom_bar(stat="identity") +
  theme(text = element_text(size=12)) +
  xlab('Confidence') +
  ylab('log10(Number of edges)') +
  theme_bw() +
  theme(legend.position = "none")

net_DORO_hs_cdkn1a <- subset(net_DORO_hs, target == "CDKN1A" & confidence %in% c("A", "B", "C")) # 42 obs


##====================== dorothea v1.7.3 'get_dorothea' ==========================================================
# https://saezlab.github.io/dorothea/articles/dorothea.html
# http://bioconductor.org/packages/release/data/experiment/manuals/dorothea/man/dorothea.pdf

library(ggplot2)
library(dplyr)

net_DORO <- decoupleR::get_dorothea(organism = "human", levels = c('A', 'B', 'C', 'D')) # 276731 obs

# visualize edge confidence level of target CDKN1A 
net_DORO_cdkn1a_all <- subset(net_DORO, target == "CDKN1A")
n_edges <- net_DORO_cdkn1a_all %>%
  group_by(confidence) %>%
  summarize(n = n())
ggplot(data=n_edges, aes(x=confidence, y=log10(n), color=confidence, fill=confidence)) +
  geom_bar(stat="identity") +
  theme(text = element_text(size=12)) +
  xlab('Confidence') +
  ylab('log10(Number of edges)') +
  theme_bw() +
  theme(legend.position = "none")

net_DORO_cdkn1a <- subset(net_DORO, target == "CDKN1A" & confidence %in% c("A", "B", "C")) # 47 obs
# To better estimate TF activities, recommend to select regulons from the confidence levels A, B and C

read.csv("data/GarciaAlonso_supplemental_table_S3_regulonsNormal.csv")

## ========================= output TF-target list =======================================================

# don't use the abosolute path, use relative path
# put the script in the same file with data (in a subdirectory)
# organize project

TF <- read.csv("TF.csv", header = FALSE)
target_query <- "CDKN1A"

literature_found <- tibble(
  source = TF, 
  target = target_query,
  literature_curated = "found")

literature_found <- as.data.frame(literature_found)

## ========================= collecTRI ==================================================

net_TRI_cdkn1a # TF of CDKN1A from collecTRI

net_TRI_cdkn1a_rec <- net_TRI_cdkn1a |> 
  full_join(literature_found, join_by(source == source, target == target))

net_TRI_cdkn1a_rec <- as.data.frame(net_TRI_cdkn1a_rec)

net_TRI_cdkn1a_rec <- net_TRI_cdkn1a_rec |> 
  mutate(literature_in_collecTRI = ifelse(is.na(mor), FALSE, ifelse(is.na(net_TRI_cdkn1a_rec$literature_curated), NA, TRUE))) |> 
  arrange(literature_curated, desc(literature_in_collecTRI))

view(net_TRI_cdkn1a_rec)

TF_TRI_cdkn1a <- net_TRI_cdkn1a_rec$source # 258 TFs

# get target list for each TF of CDKN1A from collecTRI
TF_net_TRI_cdkn1a <- subset(net_TRI, source == TF_TRI_cdkn1a)

# Use split() to create a list of targets for each TF
TF_list_TRI_cdkn1a <- split(TF_net_TRI_cdkn1a[, 2], TF_net_TRI_cdkn1a[, 1]) # 78 TFs


## ======================== DoRothEA ==========================================================

net_DORO_cdkn1a # TF of CDKN1A from DoRothEA

net_DORO_cdkn1a_rec <- net_DORO_cdkn1a |> 
  full_join(literature_found, join_by(source == source, target == target))

net_DORO_cdkn1a_rec <- as.data.frame(net_DORO_cdkn1a_rec)

net_DORO_cdkn1a_rec <- net_DORO_cdkn1a_rec |> 
  mutate(literature_in_DoRothEA = ifelse(is.na(mor), FALSE, ifelse(is.na(net_DORO_cdkn1a_rec$literature_curated), NA, TRUE))) |> 
  arrange(literature_curated, desc(literature_in_DoRothEA))

view(net_DORO_cdkn1a_rec)

TF_DORO_cdkn1a <- net_DORO_cdkn1a_rec$source # 65 TFs

# get target list for each TF of CDKN1A from DoRothEA 
# ONLY 35 TFs remain, since evidence level "D" "E" excluded
# matrix # 219 rows
TF_net_DORO_cdkn1a <- subset(net_DORO, source %in% TF_DORO_cdkn1a & confidence %in% c("A", "B", "C"))

# exclude genes not included in expr
in_target <- TF_net_DORO_cdkn1a$target %in% rownames(expr)
# 212 rows
TF_matrix_DORO_cdkn1a <- TF_net_DORO_cdkn1a[in_target, ]

# splits the data based on the values in column 3 (target) of tf.matrix and creates a list 
# where each element corresponds to a unique value in column 1 (TF). 
# The function inside lapply() extracts the target values from each split group.
TF_list_DORO_cdkn1a <- lapply(
  split(TF_matrix_DORO_cdkn1a[, 3], TF_matrix_DORO_cdkn1a[, 1]),
  function (d) d$target
)
