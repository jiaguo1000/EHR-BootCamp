source("fastSNF.R")

PsDF = function(predictor, y, test_idx){
  n_all = nrow(predictor)
  knn = round(n_all/2)
  iterN = 10
  
  # random_order = sample(1:(n_valid+n_train), n_valid+n_train, replace = F)
  M = unique(str_extract(colnames(predictor), "[^_]+"))
  orig_data = list()
  for (i in 1:length(M)) {
    orig_data[[i]] = predictor[,str_detect(colnames(predictor), M[i])]
    orig_data[[i]] = as.matrix(standardNormalization(orig_data[[i]]))
  }
  
  similarity = list()
  for (i in 1:length(M)) {
    similarity[[i]] = (dist2(orig_data[[i]], orig_data[[i]]))^0.5
    similarity[[i]] = affinityMatrix(similarity[[i]], K = knn, sigma = 0.5)
  }
  
  # SNF_result = fastSNF(similarity, K = knn, t = iterN, test_idx = test_idx, R_version = "open")
  SNF_result = fastSNF(similarity, K = knn, t = iterN, test_idx = test_idx, R_version = "base")
  
  # test --------------------------------------------------------------------
  tran_idx = setdiff(1:n_all, test_idx)
  tran_case_idx = intersect(tran_idx, which(y==1))
  tran_ctrl_idx = intersect(tran_idx, which(y==0))
  Y_test_pred_PsDF_t = rep(NA, length(test_idx))
  for (i in 1:length(test_idx)) {
    x = SNF_result[test_idx[i],]
    caseIndex = x[tran_case_idx]
    ctrlIndex = x[tran_ctrl_idx]
    Y_test_pred_PsDF_t[i] = t.test(caseIndex, ctrlIndex)$statistic
  }
  
  return(Y_test_pred_PsDF_t)
}




