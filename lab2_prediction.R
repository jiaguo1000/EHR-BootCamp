library(tidyverse)
library(glmnet)
library(ranger)
# library(xgboost)
library(pROC)

source("cv.R")

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

# predict hypertension ----------------------------------------------------
i = 1
colnames(outcome)[i] == "ICD_4019"
y = outcome[,i]
case_idx = which(y==1)
ctrl_idx = which(y==0)

# split training data (for cv) and new data (for test)
set.seed(2021)
train_case_idx = sample(case_idx, length(case_idx)/2)
train_ctrl_idx = sample(ctrl_idx, length(ctrl_idx)/2)
new_case_idx = setdiff(case_idx, train_case_idx)
new_ctrl_idx = setdiff(ctrl_idx, train_ctrl_idx)

train_idx = c(train_case_idx, train_ctrl_idx)
new_idx = c(new_case_idx, new_ctrl_idx)
train_X = Matrix(predictor[train_idx,], sparse = TRUE)
new_X = Matrix(predictor[new_idx,], sparse = TRUE)
train_y = y[train_idx]
new_y = y[new_idx]

# LASSO cv
lasso_cv = cv.glmnet(x = train_X, y = train_y, family = "binomial", alpha = 1, nfolds = 5)
lasso_coef = as.matrix(coef(lasso_cv, s = lasso_cv$lambda.min))
lasso_pred = predict(lasso_cv, newx = new_X, s = lasso_cv$lambda.min, type = "response")
lasso_AUC = auc(new_y~lasso_pred, quiet=T)
lasso_AUC

# RF cv
p = ncol(predictor)
m_list = c(sqrt(p)-20, sqrt(p)-10, sqrt(p), sqrt(p)+10, sqrt(p)+20)
m_list = c(sqrt(p)-20, sqrt(p))
nfold = 5
train_case_idx = which(train_y==1)
train_ctrl_idx = which(train_y==0)
case_cv = cross_validation(train_case_idx, nfold)
ctrl_cv = cross_validation(train_ctrl_idx, nfold)

res = NULL
for (j in 1:nfold) {
  cv_valid_idx = c(case_cv[[j]], ctrl_cv[[j]])
  cv_train_idx = setdiff(1:length(train_y), cv_valid_idx)
  
  cv_train_X = Matrix(train_X[cv_train_idx,], sparse = TRUE)
  cv_valid_X = Matrix(train_X[cv_valid_idx,], sparse = TRUE)
  cv_train_y = train_y[cv_train_idx]
  cv_valid_y = train_y[cv_valid_idx]
  
  for (m in m_list) {
    RF_model = ranger(x = cv_train_X, y = cv_train_y, num.trees = 100, mtry = m, probability = TRUE)
    cv_pred = predict(RF_model, data = cv_valid_X)$predictions[,2]
    cv_AUC = auc(cv_valid_y~cv_pred, quiet=T)
    tmp = tibble(cv=j, m=m, AUC=cv_AUC)
    res = rbind(res, tmp)
  }
}
m_best = res %>% 
  group_by(m) %>% 
  mutate(AUC=mean(AUC)) %>% 
  ungroup() %>% 
  select(-cv) %>% 
  distinct()

m_best = m_best$m[m_best$AUC==min(m_best$AUC)]

RF_model = ranger(x = train_X, y = train_y, num.trees = 100, mtry = m_best, probability = TRUE)
RF_pred = predict(RF_model, data = new_X)$predictions[,2]
RF_AUC = auc(new_y~RF_pred, quiet=T)
RF_AUC


