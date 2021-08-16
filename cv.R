cross_validation = function(full_list, fold){
  cv = list()
  size = length(full_list)
  remaining = full_list
  for(i in 1:(fold-1)){
    cv[[i]] = sample(remaining, size = size/fold, replace = F)
    remaining = setdiff(remaining, cv[[i]])
  }
  cv[[fold]] = remaining
  return(cv)
}