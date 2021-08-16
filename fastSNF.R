library(SNFtool)
library(SMUT)

fast_dominateset = function(xx, KK=20, test_idx) {
  ###This function outputs the top KK neighbors.	
  
  zero = function(x) {
    x[test_idx] = min(x)
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize = function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
  }
  
  return(normalize(A))
}

fastSNF = function (Wall, K = 20, t = 20, test_idx, R_version) {
  # message("start fastSNF...")
  check_wall_names = function(Wall) {
    name_match = function(names_A, names_B) {
      return(identical(dimnames(names_A), dimnames(names_B)))
    }
    return(all(unlist(lapply(Wall, FUN = name_match, Wall[[1]]))))
  }
  wall.name.check = check_wall_names(Wall)
  wall.names = dimnames(Wall[[1]])
  if (!wall.name.check) {
    warning("Dim names not consistent across all matrices in Wall.\n            Returned matrix will have no dim names.")
  }
  LW = length(Wall)
  normalize = function(X, test_idx) {
    temp = X
    temp[,test_idx] = 0
    diag(temp) = diag(X)
    row.sum.mdiag = rowSums(temp) - diag(X)
    row.sum.mdiag[row.sum.mdiag == 0] = 1
    X = X/(2 * (row.sum.mdiag))
    diag(X) = 0.5
    return(X)
  }
  newW = vector("list", LW)
  nextW = vector("list", LW)
  for (i in 1:LW) {
    Wall[[i]] = normalize(Wall[[i]], test_idx)
    Wall[[i]] = (Wall[[i]] + t(Wall[[i]]))/2
  }
  for (i in 1:LW) {
    newW[[i]] = (fast_dominateset(Wall[[i]], K, test_idx))
  }
  for (i in 1:t) {
    for (j in 1:LW) {
      sumWJ = matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
      for (k in 1:LW) {
        if (k != j) {
          sumWJ = sumWJ + Wall[[k]]
        }
      }
      if(R_version=="base"){nextW[[j]] = eigenMapMatMult(eigenMapMatMult(newW[[j]], (sumWJ/(LW - 1))), t(newW[[j]]))}
      if(R_version=="open"){nextW[[j]] = newW[[j]] %*% (sumWJ/(LW - 1)) %*% t(newW[[j]])}
    }
    for (j in 1:LW) {
      Wall[[j]] = normalize(nextW[[j]], test_idx)
      Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2
    }
  }
  W = matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
  for (i in 1:LW) {
    W = W + Wall[[i]]
  }
  W = W/LW
  W = normalize(W, test_idx)
  W = (W + t(W))/2
  if (wall.name.check) {
    dimnames(W) = wall.names
  }
  return(W)
}

