# `simu_or`: Two-stage order-restrcited B+/B- design

## Description:

The `simu_or` function is designed to simulate a order-restricted B+/B- design under given design parameters. It generates simulated data based on the provided parameters and returns the outcome of the simulated trial. It generates one instance of the trial. Requires package `sn`.

## Usage:

simu_or(mu_pos, mu_neg, delta, N, N2en, N2ep, prev, t, c1n, c1p, eff_bound, truePrev)

## Function Parameters:

1. `mu_pos`: Mean treatment effect for biomarker positive (B+) subgroup.
2. `mu_neg`: Mean treatment effect for biomarker positive (B-) subgroup.
3. `delta`: Effect size assuming a standard deviation of 1.
4. `N`: Integer, total number of patients per arm if both B+ and B- are enrolled in the second stage.
5. `N2en`: Integer, number of B- patients enriched per arm in the second stage.
6. `N2ep`: Integer, number of B+ patients enriched per arm in the second stage.
7. `prev`: Sampling prevalence of the biomarker positive subgroup.
8. `t`: Information fraction expected at the interim analysis.
9. `c1n`: Futility cutoff for the B- subgroup.
10. `c1p`: Futility cutoff for the B+ subgroup.
11. `eff_bound`: Numeric nonnegative vector length of 2, the efficacy bound used in terms of significance level. User should verify if the efficacy bound controls the type I error rate.
12. `truePrev`: True prevalence of the biomarker positive subgroup.

## Value:

1. `PET`: Probability of early termination for both subgroups.
2. `PET_BN`: Probability of early termination for B- subgroup only.
3. `PET_Eff`: Probability of early termination for both subgroups due to efficacy stopping.
4. `N`: Total sample size enrolled.
5. `M`: Total number of patients sampled.
6. `rej_BP`: Indicator of whether hypothesis regarding B+ is rejected.
7. `rej_BN`: Indicator of whether hypothesis regarding B- is rejected.

## Examples:
library(sn)

simu_or(mu_pos = 0.3, mu_neg= 0.3, delta= 0.3, N=300, N2en=150, N2ep=150, prev=0.5, t=0.5, c1n=1.2, c1p=1.4, eff_bound=c(0.0025, 0.024), TruePrev = 0.6)






# `simu_order_restricted` : Repeated Simulation of Two-stage order-restrcited B+/B- design based on three scenario considered

## Description:

The `simu_order_restricted` function is designed to simulate a order-restricted B+/B- design under given design parameters for multiple times. This function calls 'simu_or' function simulating the three scenarios considered in the paper. This function requires package `sn`.

## Usage:

simu_order_restricted(stoprateN,stoprateP,eff_bound,delta,N,N2en,N2ep,prev,t,truePrev,n.sim)

## Function Parameters:

1. `stoprateN`: Ideal probability of stopping under the alternative for biomarker negative (B-) subgroup.
2. `stoprateP`: Ideal probability of stopping under the alternative for biomarker positive (B+) subgroup.
3. `eff_bound`: Numeric nonnegative vector length of 2, the efficacy bound used in terms of significance level. User should verify if the efficacy bound controls the type I error rate.
4. `delta`: Effect size assuming a standard deviation of 1.
5. `N`: Integer, total number of patients per arm if both B+ and B- are enrolled in the second stage.
6. `N2en`: Integer, number of B- patients enriched per arm in the second stage.
7. `N2ep`: Integer, number of B+ patients enriched per arm in the second stage.
8. `prev`: Sampling prevalence of the biomarker positive subgroup.
9. `t`: Information fraction expected at the interim analysis.
10. `truePrev`: True prevalence of the biomarker positive subgroup.
11. `n.sim`: Numbers of copies of trials desired. 

## Value:

Three vectors length of 7 will be returned. It is the average of the 7 returned value from `simu_or`. 

## Examples:
library(sn)

simu_order_restricted(stoprateN = 0.1, stoprateP = 0.1, eff_bound=c(0.0025, 0.024),delta= 0.3, N=300, N2en=150, N2ep=150, prev=0.5, t=0.5,TruePrev = 0.6,n.sim = 100000)
