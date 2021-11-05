## load data and pkg
require(tcrl)
require(matrixcalc)
source('../auxfuns/simudat.R')
seqBank = readRDS('../dat/tcrBank/TCR_AA_seq_10k.rds')
seqhom_B = readRDS('../dat/tcrBank/seqhom_BLOSUM62.rds')
seeds = readRDS('../dat/seed1M.rds')
W = readRDS('../dat/KDH_Wmat.rds')

## parameter 
coefb = matrix(c(0.1,0.5,-0.4),ncol=1)
eta = 1.6
tau = 0
n = 350
rtype = 'bin'
B = 1000

out = matrix(NA,B,3)
for(i in 1:B){
    dat = simudat(n,coefb,type=rtype,eta=eta,tau=tau,W=W,seqBank = seqBank,seqhom=seqhom_B,seed=seeds[i])
    out1 = TCRseq_bin(dat$Y,dat$X,dat$S)
    out2 = TCRL_bin(dat$Y,dat$X,dat$fR,dat$W,dat$S)
    out[i,] = c(out2[1],out1,out2[3])
}

