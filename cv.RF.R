library(ranger)
library(pROC)
source("cv.R")

cv.RF = function(train_X, train_y, m_list, nfold){
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
  m_best = m_best$m[m_best$AUC==max(m_best$AUC)]
  RF_model = ranger(x = train_X, y = train_y, num.trees = 100, mtry = m_best, probability = TRUE)
  return(list(model=RF_model, m_best=m_best))
}