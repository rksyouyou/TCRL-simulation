## build TCR bank
rm(list=ls())
source('../auxfuns/pkg.R')
 
fAA <- function(len)
    unlist(lapply(len,function(x) paste(sample(AA_STANDARD,x,TRUE),collapse = "")))


seqhom <- function(seq1, seq2,type){
    s12 <- Biostrings::pairwiseAlignment(seq1, seq2, substitutionMatrix = type)
    return(s12@score)    
}

###################################
######## generate sequence ########
###################################
hinfo= readRDS('../dat/simu/head.rds')
tinfo = readRDS('../dat/simu/tail.rds')

set.seed(510)
n = 10000 # bank size

hseg <- sample(hinfo$seq,n,TRUE,prob=hinfo$prob) # head segment
tseg <- sample(tinfo$seq,n,TRUE,prob=tinfo$prob) # tail segment

## V, D, J segments in middle seg
mlen <- sample(4:10,n,TRUE)
vlen <- ifelse(mlen==4,1,ifelse(mlen>7,3,2))
dlen <- ifelse(mlen<6,1,ifelse(mlen>8,3,2))
jlen <- ifelse(mlen <7,2,ifelse(mlen>9,4,3))

Vseg <- fAA(vlen)
Dseg <- fAA(dlen)
Jseg <- fAA(jlen)
mstatus <- sample(0:3,n,TRUE,c(0.96,0.01,0.01,0.02))
VD <- unlist(lapply(mstatus,function(x) ifelse(x==1|x==3,sample(AA_STANDARD,1),"")))
DJ <- unlist(lapply(mstatus,function(x) ifelse(x==2|x==3,sample(AA_STANDARD,1),"")))

## combine 
TCR <- paste0(hseg,Vseg,VD,Dseg,DJ,Jseg,tseg)
TCRlen <- as.numeric(lapply(TCR,nchar))

saveRDS(TCR,'../dat/tcrBank/TCR_AA_seq_10k.rds')

####################################
######## compute homology  #########
####################################
b = 1 
type = 'BLOSUM62'
TCR <- readRDS('../dat/tcrBank/TCR_AA_seq_10k.rds') # seq data    
n = length(TCR)
omat = matrix(NA,n,n)
for(j in 1:n){
    for(i in j:n)
        omat[i,jp] = seqhom(TCR[i],TCR[j],type)
}
saveRDS(omat,'../dat/tcrBank/scoreMat_BLOSUM62.rds') 

## use score matrix to compute homology
scoreMat <- readRDS('../dat/tcrBank/scoreMat_BLOSUM62.rds')
smat1 <- scoreMat/sqrt(diag(scoreMat))
smat2 <- t(t(smat1)/sqrt(diag(scoreMat)))
smat2[is.na(smat2)] = 0
smat3 = smat2+t(smat2)
diag(smat3) = 1
saveRDS(smat3,'../dat/tcrBank/seqhom_BLOSUM62.rds')





