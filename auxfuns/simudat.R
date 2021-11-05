## simulate data
simudat <- function(n,coefb,type,eta,tau,W,seqBank,seqhom,seed){
    set.seed(seed)

    ## convariates
    X = cbind(1,rbinom(n,1,0.5),rnorm(n))
    mu = X%*%coefb
    
    ## simu TCR
    sTCR = sample(2:25,size = n,TRUE) # size of TCR
    idTCR = lapply(sTCR,function(x) sample(1:10^4,x)) # sample TCR id
    abundance = lapply(sTCR,function(x) sample(1:5,x,TRUE))

    ## extract feature
    fR = t(mapply(AA_freq,idTCR,abundance,MoreArgs = list(seqBank=seqBank))) # f(R)
    fRW = fR%*%W

    ## homology matrix
    S  <- matrix(NA,n,n)
    for(i in 1:(n-1))
        for(j in (i+1):n){
            S[i,j] <- S[j,i] <-  Sij(idTCR[[i]],idTCR[[j]],abundance[[i]],abundance[[j]],seqhom)
        }
    diag(S)  = 1
    if(!is.positive.semi.definite(S)) S <- psd(S)
    
    if(tau!=0){
        hR <- mvrnorm(1,mu=rep(0,n),Sigma = tau^2*S)
        mu = mu + fRW*eta  + hR
        if(type == 'bin'){
            Y = rbinom(n,1,exp(mu)/(exp(mu)+1))
        }else {
            Y = rnorm(n,mu,1)
        }
        Y = cbind(Y,Y_p)
    } else {
        mu = mu +fRW*eta
        if(type == 'bin') Y = rbinom(n,1,exp(mu)/(exp(mu)+1)) else Y = rnorm(n,mu,1)
        Y = matrix(Y,ncol=1)
    }
    
    ## output
    dat = list(Y=Y,X=X,S=S,fR=fR,W=W,idTCR=idTCR,abundance=abundance,seed=seed)
    return(dat)
    
}


AA_freq <- function(id,abund,seqBank){ ## extract features
    sAA = strsplit(paste(rep(seqBank[id],abund),collapse=""),"")[[1]]
    freq = table(factor(sAA,levels=AA_STANDARD))/length(sAA)
    return(freq)
}


Sij <- function(id1,id2,abund1,abund2,seqhom){
    submat = seqhom[id1,id2]
    s1 = sum(apply(submat,1,max)*abund1)
    s2 = sum(apply(submat,2,max)*abund2)
    (s1+s2)/sum(c(abund1,abund2))
}


psd <- function(S){
  if (isSymmetric(S)==FALSE){
    stop("S should be symmetric.")
  }
  es = eigen(S)
  S <- as.matrix(S)
  Q <- es$values
  m <- sum(Q>=0)
  P <- es$vectors
  n <- dim(S)[1]
  my.psd <- matrix(0, n, n)
  rownames(my.psd) <- rownames(S)
  colnames(my.psd) <- colnames(S)
  for (j in 1: m){
    my.psd <- my.psd + Q[j]*P[,j]%*%t(P[,j])
  }
  return(my.psd)
}

