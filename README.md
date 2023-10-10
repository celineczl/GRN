# Gene Regulatory Network Inference

## Project Objectives
To develop and validate a computational method for inferring gene regulatory networks from genomic data using a non-local prior Bayesian model as a foundation for high-dimensional variable selection. Other statistical models and variable selection methods (lasso and mombf) are tested and compared.

## Set-up
### Github
`git clone git@github.com:Celine-ZL-Chen/GRN.git`
### R Packages
* io: read from, write to, plot to in a unified manner
* ggplot2
* glmnet
* mombf
* dplyr
* to be completed

## Databases
1. GTEx
2. CGTA
3. DoRothEA
4. to be completed

## Data Preprocessing
GTEx ...
doro_tf-target_filtered.R

## Model Description
### Class 1: assume TF expression = TF activity

Model 0 - univariate linear model

`uni_reg_cdkn1a.R`

Model 1 - multivariate linear model
'multi_reg_cdkn1a.R'

Model 2 - multivariate lasso model (glmnet)
'lasso_cdkn1a.R'

Model 3 - multivariate non-local prior model (mombf)
'mombf_cdkn1a.R'

‌Class 2: infer TF activity from target genes' expression
* mean expression of target genes (DoRothEA) as activity of TF [no normalization performed]
* Summerized GTEX data [pick out all activities of the TF]
