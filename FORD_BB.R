## mu_pos and mu_neg are true treatment effect of biomarker positive and negative groups
## mu is the effect size used to power the study
## sigma is the common variance assumed for the treatment effect across patients
## N is the total sample size per arm
## prev is the prevalence of biomarker positive subgroup
## t is the information fraction at the interim analysis
## c1n and c1p is the interim analysis futility cutoffs
## eff_bound takes in a vector of length 2 specifying the efficacy bound of users choice
FORD_BB = function(mu_pos,mu_neg,sigma,N,prev,t,c1n,c1p,eff_bound){
  #browser()
  PET = PET.n = PETE =  SS= rejpos = rejneg =0
  req = 0
  Nn = round(N * (1-prev))
  Np = N - Nn
  ## First stage sample size
  N1n = round(Nn * t)
  N1p = round(Np * t)
  pos1 = rnorm(N1p,mu_pos,sqrt(2)*sigma) 
  neg1 = rnorm(N1n,mu_neg,sqrt(2)*sigma)
  posbar1=mean(pos1)
  negbar1 = mean(neg1)
  pool1 = mean(c(pos1,neg1))
  sigmapos1 = sd(pos1)
  sigmaneg1 = sd(neg1)
  sigmapool1 = mean(c(sigmapos1,sigmaneg1))
  ## Second stage sample size and statistics  
  N2p = Np - N1p
  N2n = Nn - N1n
  N2e = N - N1n-N1p
  pos2 = rnorm(N2p,mu_pos,sqrt(2)*sigma) 
  neg2 = rnorm(N2n,mu_neg,sqrt(2)*sigma)
  pool2 = c(pos2,neg2)
  pose = rnorm(N2e,mu_pos,sqrt(2)*sigma) 
  posbare = mean(pose)
  sigmapose = sd(pose)
  nege = rnorm(N2e,mu_neg,sqrt(2)*sigma) 
  negbare = mean(nege)
  sigmanege = sd(nege)
  test.stat.neg2 = negbare/(sqrt(2/N2e))
  test.stat.pos.ia = (posbar1-0.3)/(sqrt(2/N1p))
  test.stat.neg.ia = (negbar1-0.3)/(sqrt(2/N1n))
  test.stat.pool.ia = (pool1-0.3)/(sqrt(2/(N1n+N1p)))
  test.stat.pos.ia1 = (posbar1-0)/(sqrt(2/N1p))
  test.stat.neg.ia1 = (negbar1-0)/(sqrt(2/N1n))
  test.stat.pool.ia1 = (pool1-0)/(sqrt(2/(N1n+N1p)))
  ## Enrichment
  test.stat.pos2 = posbare/(sqrt(2/N2e))
  #Go both
  test.stat.pos = mean(pos2)/(sqrt(2/N2p))
  test.stat.pool = mean(pool2)/(sqrt(2/(N2p+N2n)))
  test.stat.neg = mean(neg2)/sqrt(2/N2n)
  ts.ia = min(test.stat.neg.ia,test.stat.pool.ia)
  ts.ia1 = max(test.stat.pos.ia1,test.stat.pool.ia1)
  ts = max(test.stat.pos,test.stat.pool)
  pneg.ia = 1 - pnorm(test.stat.neg.ia1)
  pall.ia = 1-psn(ts.ia1,alpha = (1-sqrt(prev))/sqrt(1-prev))
  pall.se = 1-psn(ts, alpha = (1-sqrt(prev))/sqrt(1-prev) )
  #pe.ia = 1- psn(test.stat.pos.ia1,alpha = (1-sqrt(prev))/sqrt(1-prev))
  pneg.se = 1 - pnorm(test.stat.neg)
  pe.ia = pall.ia
  pe.se = 1- pnorm(test.stat.pos2)
  ne.se = 1 - pnorm(test.stat.neg2)
  w = t
  ne = 1 - pnorm(sqrt(w)*qnorm(1-pneg.ia)+sqrt(1-w)*qnorm(1-ne.se))
  pe = 1 - pnorm(sqrt(w)*qnorm(1-pe.ia)+sqrt(1-w)*qnorm(1-pe.se))
  pall = 1 - pnorm(sqrt(w)*qnorm(1-pall.ia)+sqrt(1-w)*qnorm(1-pall.se))
  pneg = 1 - pnorm(sqrt(w)*qnorm(1-pneg.ia)+sqrt(1-w)*qnorm(1-pneg.se))
  if(ts.ia < c1n){ # go A plus only
    PET.n = 1
    if(test.stat.pos.ia>c1p){
      if(pe.ia < eff_bound[1]){
        SS = N1n + N1p
        PETE =1
        rejpos = 1
      }else{
        SS = N
        if(pe <= eff_bound[2]){
          rejpos=1
        } 
      }
    }
    else{ # stop all for futility
      SS=N1n + N1p
      PET =1
    }
  }
  else{ 
    if(pe.ia < eff_bound[1]){ 
      PETE =1
      rejpos =  1
      if(pneg.ia < eff_bound[1]){
        rejneg = 1
        SS = N1n + N1p
      }else{
        SS = N
        if(ne <= eff_bound[2]){ #enrich A-
          rejneg=1
        }
      }
    }
    else{ #go both
      SS= N
      if(pall <= eff_bound[2]){
        rejpos =  1
        if(pneg <= eff_bound[2]){
          rejneg = 1 
        }
      }
    }
    
  }
  # returns probability of early stopping of both subgroups (PET); stopping of biomarker negative group only (PET.n)
  # PETE is the probability of any interim efficacy stopping
  # SS is the total sample size of the study, rejpos and rejneg are decisions to reject the positive and negative subgroup.
  return(c(PET,PET.n,PETE,SS,rejpos,rejneg))
}
