library(sn)

# N2e is the added patients

simu_or = function(mu_pos,mu_neg,delta,N,N2en,N2ep,prev,t,c1n,c1p,eff_bound,TruePrev){
  #browser()
  PET = PET.n = PETE = SS = screen = rejpos = rejneg =0
  Nn = round(N * (1-prev))
  Np = N - Nn
  ## First stage sample size
  N1n = round(Nn * t)
  N1p = round(Np * t)
  pos1 = rnorm(N1p,mu_pos,sqrt(2)) 
  neg1 = rnorm(N1n,mu_neg,sqrt(2))
  posbar1=mean(pos1)
  negbar1 = mean(neg1)
  pool1 = mean(c(pos1,neg1))
  ## Second stage sample size and statistics  
  N2p = Np - N1p
  N2n = Nn - N1n
  pos2 = rnorm(N2p,mu_pos,sqrt(2)) 
  neg2 = rnorm(N2n,mu_neg,sqrt(2))
  pool2 = c(pos2,neg2)
  pose = rnorm(N2ep,mu_pos,sqrt(2)) 
  posbare = mean(pose)
  nege = rnorm(N2en,mu_neg,sqrt(2)) 
  negbare = mean(nege)
  test.stat.neg2 = negbare/(sqrt(2/N2en))
  test.stat.pos.ia = (posbar1-delta)/(sqrt(2/N1p))
  test.stat.neg.ia = (negbar1-delta)/(sqrt(2/N1n))
  test.stat.pool.ia = (pool1-delta)/(sqrt(2/(N1n+N1p)))
  test.stat.pos.ia1 = (posbar1-0)/(sqrt(2/N1p))
  test.stat.neg.ia1 = (negbar1-0)/(sqrt(2/N1n))
  test.stat.pool.ia1 = (pool1-0)/(sqrt(2/(N1n+N1p)))
  ## Enrichment
  test.stat.pos2 = posbare/(sqrt(2/N2ep))
  #Go both
  test.stat.pos = mean(pos2)/(sqrt(2/N2p))
  test.stat.pool = mean(pool2)/(sqrt(2/(N2p+N2n)))
  test.stat.neg = mean(neg2)/sqrt(2/N2n)
  ts.ia = min(test.stat.neg.ia,test.stat.pool.ia)
  ts.ia1 = max(test.stat.pos.ia1,test.stat.pool.ia1)
  ts = max(test.stat.pos,test.stat.pool)
  pall.ia = 1-psn(ts.ia1,alpha = (1-sqrt(prev))/sqrt(1-prev))
  pneg.ia = 1 - pnorm(test.stat.neg.ia1)
  pall.se = 1-psn(ts, alpha = (1-sqrt(prev))/sqrt(1-prev))
  pe.ia = pall.ia
  pneg.se = 1 - pnorm(test.stat.neg)
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
        screen = max(round(N1n/(1-TruePrev)),round(N1p/TruePrev))
        PETE =1
        rejpos = 1
      }else{
        SS = N1n +N1p + N2ep
        screen = max(N1n/(1-TruePrev) +round(N2ep/TruePrev) ,N1p/TruePrev + round(N2ep/TruePrev))
        if(pe <= eff_bound[2]){
          rejpos=1
        } 
      }
    }
    else{ # stop all for futility
      SS=N1n + N1p
      screen =max(round(N1n/(1-TruePrev)),round(N1p/TruePrev))
      PET =1
    }
  }
  else{ 
    if(pe.ia < eff_bound[1]){ 
      rejpos =  1
      if(pneg.ia < eff_bound[1]){
        PETE =1
        rejneg = 1
        SS = N1n + N1p
        screen = max(round(N1n/(1-TruePrev)),round(N1p/TruePrev))
      }else{
        SS = N1p + N1n + N2en
        screen = max(round(N1n/(1-TruePrev) +N2en/(1-TruePrev)) ,round(N1p/TruePrev + N2en/(1-TruePrev)))
        if(ne <= eff_bound[2]){ #enrich A-
          rejneg=1
        }
      }
    }
    else{ #go both
      SS= N
      screen = max(round(Nn/(1-TruePrev)),round(Np/TruePrev))
      if(pall <= eff_bound[2]){
        rejpos =  1
        if(pneg <= eff_bound[2]){
          rejneg = 1 
        }
      }
    }
    
  }
  return(list(PET = PET,PET_BN = PET.n,PET_Eff = PETE,N = SS,M = screen,rej_BP = rejpos,rej_BN = rejneg))
}


simu_order_restricted = function(stoprateN,stoprateP,eff_bound,delta,N,N2en,N2ep,prev,t,TruePrev,n.sim){
  c1n = qsn(stoprateN,alpha = -(1-sqrt(1-prev))/sqrt(prev))
  c1p = qnorm(stoprateP)
  mu00 = matrix(ncol=7,nrow = n.sim)
  mu01 = matrix(ncol=7,nrow = n.sim)
  mu11 = matrix(ncol=7,nrow = n.sim)
  for(i in 1:n.sim){
    mu00[i,] = unlist(simu_or(0,0,delta,N,N2en,N2ep,prev,t,c1n,c1p,eff_bound = eff_bound,TruePrev))
    mu01[i,] = unlist(simu_or(delta,0,delta,N,N2en,N2ep,prev,t,c1n,c1p,eff_bound = eff_bound,TruePrev))
    mu11[i,] = unlist(simu_or(delta,delta,delta,N,N2en,N2ep,prev,t,c1n,c1p,eff_bound = eff_bound,TruePrev))
  }
  mu00 = colMeans(mu00)
  mu01 = colMeans(mu01)
  mu11 = colMeans(mu11)
  return(list(mu00=mu00,mu01=mu01,mu11=mu11))
}

