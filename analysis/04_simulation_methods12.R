# Reference implementation for binomial resampling Methods 1 and 2.
# This RECONSTRUCTS the current methods and DOES NOT claim to reproduce the
# original unarchived Table 2/3 random-number streams or the sequential competitor.
# Default runs a small diagnostic; use N_OUTER=5000 N_INNER=10000 for the very slow
# full-scale author check, preferably with parallel computing and reproducible seeds.
source('R/nma_fe.R'); source('R/closed_form.R'); source('R/resampling.R')
N_OUTER <- as.integer(Sys.getenv('N_OUTER','20'))
N_INNER <- as.integer(Sys.getenv('N_INNER','1000'))
TOTAL_N <- as.integer(Sys.getenv('TOTAL_N','50'))
RULE <- Sys.getenv('RULE','lt')       # 'lt' or 'band'
SCENARIO <- Sys.getenv('SCENARIO','null') # 'null' or 'alternative'
set.seed(20261029)
d <- read_brd('data/BRD_Arm_Network.csv')
f <- nma_fe(d);risk <- nma_risks(d,f)
A <- 'Ceftiofur hydrochloride'; B <- 'Trimethoprim'
nt <- attr(risk,'n_by_trt')
p0 <- (risk[A]*nt[A]+risk[B]*nt[B])/(nt[A]+nt[B])
if (SCENARIO=='null') {pA <- p0; pB <- p0} else {pA <- .35; pB <- .55}
newrisk <- risk;newrisk[A] <- pA;newrisk[B] <- pB
selection_ok <- function(p) if (RULE=='lt') p < 0.10 else p>0.05 & p<0.10
out <- matrix(NA_real_,N_OUTER,6)
i <- 0; trials <- 0
while(i<N_OUTER){
  trials <- trials+1
  d_old <- resample_events(d,newrisk)
  old <- nma_fe(d_old)
  if (!selection_ok(old$p)) next
  i <- i+1
  rA <- rbinom(1,TOTAL_N/2,pA);rB <- rbinom(1,TOTAL_N/2,pB)
  d_new <- add_trial(d_old,A,B,rA,TOTAL_N/2,rB,TOTAL_N/2)
  rA <- max(.5,min(TOTAL_N/2-.5,rA));rB <- max(.5,min(TOTAL_N/2-.5,rB))
  phA <- rA/(TOTAL_N/2);phB <- rB/(TOTAL_N/2)
  yy <- qlogis(phA)-qlogis(phB)
  vv <- 1/((TOTAL_N/2)*phA*(1-phA))+1/((TOTAL_N/2)*phB*(1-phB))
  cf <- nma_update_test(old$z,old$se^2,yy,vv,RULE)
  a <- resample_method1(d_new,A=A,B=B,criterion=RULE,n_iter=N_INNER)['p']
  b <- resample_method2(d_new,A=A,B=B,n_iter=N_INNER)['p']
  out[i,] <- c(cf['p_method1'],cf['p_method2'],cf['p_alone'],cf['p_naive'],a,b)
  if (i%%5==0) cat('selected',i,'/',N_OUTER,'attempted',trials,'\n')
}
colnames(out) <- c('closed1','closed2','alone','naive','resample1','resample2')
print(colMeans(out<0.05,na.rm=TRUE)); print(colMeans(out,na.rm=TRUE))
