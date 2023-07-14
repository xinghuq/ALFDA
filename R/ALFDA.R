### define functions
## The input dataset x is a N*L genotype matrix coding with 0,1, 2. It can be transformed to allele frequency by divided by 2 or just using the copy number of reference alleles, where N is number of individuals and L is the number of locus 
## The calculation is analogous to other distance measures(EU distance) that use the data matrix as input 
## @ x is input dataset with N*L
##@ para whether using parallel cluster to compute the parawise distance, default is to use ncore-2 cores




pwDis=function(x,para){
## The distance measures here is analogous to species in a community, locus is a site  
### transpose the dataset to a L*N matrix
  
  x=t(as.matrix(x))
# define the basic individual/pop distance function
  DeltaD = function(abun, struc) {
    # get the total number of copy of reference alleles in all individuals
    n = sum(abun,na.rm = TRUE)
    ## Number of individuals
    N = ncol(abun)
    # total number of copy of reference alleles in an individual at a locus
    ga = rowSums(abun,na.rm = TRUE)
    #get the sample average allele frequency across all individuals at a locus
    gp = ga[ga > 0]/n
    ## Entropy of a locus 
    G = sum(-gp * log(gp))
    ## If these individuals are structured in hierarchies, H is hierarchy levels defined below in "str" data
    H = nrow(struc) ## str is structured into three levels, individuals, pops, ecosystem, for individual pairwise distance, pops=ecosystem, so there are only individual and ecosystem
    A = numeric(H - 1) ##
    W = numeric(H - 1)
    Diff = numeric(H - 1) # distance/differentation
    
    wi = colSums(abun,na.rm = TRUE)/n  ## the weight/allele frequency/proportion of variants / in an individual i across all individuals. 
    W[H - 1] = -sum(wi[wi > 0] * log(wi[wi > 0]),na.rm = TRUE) # get the weighted-entropy of an individual i in the pooled population
    pi = sapply(1:N, function(k) abun[, k]/sum(abun[, k], na.rm = TRUE)) ## the allele frequency of at locus l in indivdual i
    Ai = sapply(1:N, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0]),na.rm = TRUE)) ## the entropy of an individual i
    A[H - 1] = sum(wi * Ai)  ## The entropy of the pooled population/(here is the whole individuals)
    
    ## now enter a conditional statement, if the herarchy has more than 2 levels (individual, pops, ecosystems, compuute W(H-1), unless W(1)=1
    ## individual, pops, ecosystem
    ## the calculation is the same as above
    if (H > 2) {
      for (i in 2:(H - 1)) {
        I = unique(struc[i, ]) 
        NN = length(I)
        ai = matrix(0, ncol = NN, nrow = nrow(abun))
        c
        for (j in 1:NN) {
          II = which(struc[i, ] == I[j])
          if (length(II) == 1) {
            ai[, j] = abun[, II]
          }
          else {
            ai[, j] = rowSums(abun[, II],na.rm = TRUE)
          }
        }
        pi = sapply(1:NN, function(k) ai[, k]/sum(ai[,k],na.rm = TRUE))
        wi = colSums(ai)/sum(ai) ## the proportion of variants in a individual across all variants  
        W[i - 1] = -sum(wi * log(wi))
        Ai = sapply(1:NN, function(k) -sum(pi[, k][pi[,k] > 0] * log(pi[, k][pi[, k] > 0])))
        A[i - 1] = sum(wi * Ai)
      }
    }
    #### distace or the differentiation between individuals 
    if(W[1]==0) {
      Diff[1] = (G - A[1])  
      } else{
      Diff[1] = (G - A[1])/W[1]
      }
    
    #### if H >2 calculate other levels of differentation 
    
    if (H > 2) {
      for (i in 2:(H - 1)) {
        Diff[i] = (A[i - 1] - A[i])/(W[i] - W[i - 1])
      }
    }
    Diff = Diff
    out = matrix(c(Diff), ncol = 1)
    return(out)
  }
  
  v1 = c("ecosystem", "region1", "pop1") ## This can vary but here is all levels are the same 
  v2 = c("ecosystem", "region1", "pop2")##
  str = data.frame(v1, v2)
  str = as.matrix(str)
  nvar = ncol(x)
  napops=colnames(x)
  # pairwise matrix index
  pw <- combn(nvar, 2,simplify = FALSE) #,simplify = FALSE
  
  # set up parallel cluster
  #if(para){
  
  #  ncor <- parallel::detectCores()
  #}
  #cl <- parallel::makeCluster(ncor)
  
  # Dmat = matrix(data = 0, nrow = nvar, ncol = nvar,dimnames = list(napops, napops))
  # default is using parallel computation
  if(para){
    ncor <- parallel::detectCores()
    print(paste("Using parallel",ncor-2,"cores"))
    Dp=parallel::mclapply(pw,function(i, abun) DeltaD(abun[,i],struc=str)[2], abun = x,  mc.cores = ncor-2)
  }
  else {
    print(paste("Using normal computer resource"))
    Dp=lapply(pw,function(i, abun) DeltaD(abun[,i],struc=str)[2], abun = x)
  }
  ## set a distance pairwise matrix with diagonal is 0
  pwD <- diag(nvar)
  diag(pwD) <- 0
  pwD[lower.tri(pwD)] <- unlist(Dp)
  # this ugly looking bit is due to the matrix filling up by columns instead of rows
  pwD[upper.tri(pwD)] <- t(pwD)[upper.tri(t(pwD))]
  colnames(pwD)  <- colnames(x)
  row.names(pwD)=colnames(x)
  
  
  #  Dmat=lapply(x,DeltaD()[2] i,j)
  #  b=parallel::parLapply(cl,pw,function(i,abun) HierDpart::IDIP(abun[,i],struc=str), abun=dat)
  #x[i]
  
  # parallel::stopCluster(cl)
  
  ###  use this loop only for small data
  
  #   for (i in 2:nvar) {
  #   for (j in 1:(i-1)) {
  #    Dmat[i, j] =Dmat[j, i]=DeltaD((x[, c(i, j)]),str)[2]
  #}
  #  }
  #  Dm=parallel::parLapply(x,DeltaD(x,str)[2], i=2:nvar,j=1:(i-1))
  #parallel::stopCluster(cl)
  
  # pairwiseDav=as.matrix(Dmat)
  #colnames(pairwiseDav) = colnames(x)
  #rownames(pairwiseDav) = colnames(x)
  ##DeltaDmat = as.dist(pairwiseDav, diag = FALSE, upper = FALSE)
  # or may try sapply(1:N, function(x) DeltaD((x[, c(i, j)]), str)[2],)
  
  return(list(PairwiseD = pwD))
}

##
### compute the affinity matrix that represents the neighborhood graph of the individuals, Wang et al, 2012

## pwd here is 1-pwdistance 

### search the K nearest neighborhood and assign the weight to each individuals

Af=function (pwD, K = 10, sigma = 0.5) 
{
  N <- nrow(pwD)
  pwD <- (pwD + t(pwD))/2
  diag(pwD) <- 0
  sortedColumns <- as.matrix(t(apply(pwD, 2, sort)))
  finiteMean <- function(x) {
    return(mean(x[is.finite(x)]))
  }
  means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
    .Machine$double.eps
  avg <- function(x, y) {
    return((x + y)/2)
  }
  Sig <- outer(means, means, avg)/3 * 2 + pwD/3 + .Machine$double.eps
  Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
  densities <- dnorm(pwD, 0, sigma * Sig, log = FALSE)
  W <- (densities + t(densities))/2
  return(W)
}



require("lfda")

### Here is KLFDA function, now using affinity matrix as input 
## genmat : input matrix is a genotype matrix N*L

ALFDA=function (genmat, y, r, kaf=10,sigma=0.5, metric = c("weighted", "orthonormalized",
                                    "plain"), knn = 6, reg = 0.001)
{
  
  ShannonDis=pwDis(as.matrix(genmat),para=T)
  k=Af(as.matrix(1-ShannonDis$PairwiseD),K = kaf, sigma = sigma)
  
  metric <- match.arg(metric)
  #D1=pwDeltaD(x,para)$PairwiseDeltaD
  
  # require("SNFtool")
  #k=Af(1-D1,K = K, sigma = sigma)
  y <- t(as.matrix(y))
  n <- nrow(k)
  if (is.null(r))
    r <- n
  tSb <- mat.or.vec(n, n)
  tSw <- mat.or.vec(n, n)
  for (i in unique(as.vector(t(y)))) {
    Kcc <- k[y == i, y == i]
    Kc <- k[, y == i]
    nc <- nrow(Kcc)
    Kccdiag <- diag(Kcc)
    distance2 <- repmat(Kccdiag, 1, nc) + repmat(t(Kccdiag),
                                                 nc, 1) - 2 * Kcc
    A <- getAf(distance2, knn, nc)
    Kc1 <- as.matrix(rowSums(Kc))
    Z <- Kc %*% (repmat(as.matrix(colSums(A)), 1, n) * t(Kc)) -
      Kc %*% A %*% t(Kc)
    tSb <- tSb + (Z/n) + Kc %*% t(Kc) * (1 - nc/n) + Kc1 %*%
      (t(Kc1)/n)
    tSw <- tSw + Z/nc
  }
  K1 <- as.matrix(rowSums(k))
  tSb <- tSb - K1 %*% t(K1)/n - tSw
  tSb <- (tSb + t(tSb))/2
  tSw <- (tSw + t(tSw))/2
  F=tSb/tSw
  eigTmp <- suppressWarnings(rARPACK::eigs(A = solve(tSw +
                                                       reg * diag(1, nrow(tSw), ncol(tSw))) %*% tSb, k = r,
                                           which = "LM"))
  eigVec <- Re(eigTmp$vectors)
  eigVal <- as.matrix(Re(eigTmp$values))
  Tr <- getMetricOfType(metric, eigVec, eigVal, n)
  Z <- t(t(Tr) %*% k)
  out <- list(T = Tr, Z = Z,eigVec=eigVec,eigVal=eigVal,nc=nc,A=k,F=F,distance2=distance2)
  class(out) <- "alfda"
  return(out)
}




