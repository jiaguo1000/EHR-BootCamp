library(tidyverse)
library(ranger)
library(pROC)
library(gplots)
library(SNFtool)

source("cv.RF.R")
source("PsDF.R")

# simulation --------------------------------------------------------------
set.seed(2021)
n = 200
p = 50
case_X = matrix(rbinom(n/2*p, size = 1, prob = 0.6), nrow = n/2)
ctrl_X = matrix(rbinom(n/2*p, size = 1, prob = 0.2), nrow = n/2)
data_X = rbind(case_X, ctrl_X)

std_X = as.matrix(standardNormalization(data_X))
dist_X = as.matrix(dist(std_X, method = "euclidean"))
smlt_X = affinityMatrix(dist_X)
diag(smlt_X) = NA

hm_col = colorRampPalette(c("white","dodgerblue3","dodgerblue4"))
sample_col = c(rep("orange",n/2), rep("purple",n/2))

heatmap.2(smlt_X, srtCol=15, cexCol=1, Colv=F, Rowv=F,
          labRow=c(rep(NA, 200)), labCol=c(rep(NA, 200)),
          RowSideColors=sample_col, col=hm_col, dendrogram = "none", trace = "none",
          main="Similarity between 100 cases and 100 controls")
par(lend = 1)
legend("top", legend=c("cases", "controls"), col=c("orange", "purple"), lty= 1, lwd=10)

# read data ---------------------------------------------------------------
outcome = readRDS("Data/outcome.rds")
predictor = readRDS("Data/predictor.rds")
colnames(predictor) = str_replace_all(colnames(predictor), "-", "_")
colnames(predictor) = str_replace_all(colnames(predictor), "/", "_")
colnames(predictor) = str_replace_all(colnames(predictor), "%", "_")

# rank by prev ------------------------------------------------------------
prev = apply(outcome, 2, sum)
prev = sort(prev, decreasing = T)
outcome = outcome[,match(names(prev), colnames(outcome))]
identical(rownames(outcome), rownames(predictor))

# dict for disease --------------------------------------------------------
dict_ICD = read_csv("Data/D_ICD_DIAGNOSES.csv")
dict_outcome = tibble(ICD9_CODE = str_remove_all(colnames(outcome), "ICD_")) %>% 
  left_join(dict_ICD) %>% 
  select(-ROW_ID)

# predict hypertensive encephalopathy -------------------------------------
set.seed(2021)
i = 986
colnames(outcome)[i] == "ICD_4372"
y = outcome[,i]
case_idx = which(y==1)
ctrl_idx = which(y==0)
ctrl_idx = sample(ctrl_idx, 200)

# split training data (for cv) and new data (for test)
train_case_idx = sample(case_idx, length(case_idx)/2)
train_ctrl_idx = sample(ctrl_idx, length(ctrl_idx)/2)
new_case_idx = setdiff(case_idx, train_case_idx)
new_ctrl_idx = setdiff(ctrl_idx, train_ctrl_idx)

train_idx = c(train_case_idx, train_ctrl_idx)
new_idx = c(new_case_idx, new_ctrl_idx)
train_X = predictor[train_idx,]
new_X = predictor[new_idx,]
train_y = y[train_idx]
new_y = y[new_idx]

# RF cv
p = ncol(predictor)
m_list = c(sqrt(p)-20, sqrt(p)-10, sqrt(p), sqrt(p)+10, sqrt(p)+20)
m_list = c(sqrt(p)-20, sqrt(p))
nfold = 5

RF_cv = cv.RF(train_X=train_X, train_y=train_y, m_list=m_list, nfold=nfold)
RF_cv$m_best
RF_pred = predict(RF_cv$model, data = new_X)$predictions[,2]
RF_AUC = auc(new_y~RF_pred, quiet=T)
RF_AUC

# PsDF prediction ---------------------------------------------------------
data_X = rbind(train_X, new_X)
data_y = c(train_y, new_y)
test_idx = (length(train_y)+1):length(data_y)

PsDF_pred = PsDF(predictor = data_X, y = data_y, test_idx = test_idx)
PsDF_AUC = auc(new_y~PsDF_pred, quiet=T)
PsDF_AUC

