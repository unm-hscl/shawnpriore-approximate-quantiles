# Approximate Quantiles
Code for the paper, "Approximate Quantiles for Stochastic Optimal Control of LTI Systems with Arbitrary Disturbances," ACC 2022.

## Requirements
* SReachTools [https://unm-hscl.github.io/SReachTools/](https://unm-hscl.github.io/SReachTools/)
* CVX [http://cvxr.com/cvx/](http://cvxr.com/cvx/)

## Examples
### 4d CWH with Cauchy Disturbance
We use a 4d CWH system with additive Cauchy noise to model the dynamics of three satellites performing rendezvous operations while maintaining minimum distances from eachother. This example compares the efficacy of the numerical quantile with an alaytical quantile. The results show the numerical quantiles result in nearly identical trajectories. Thus, implying the error of the numerical quantile is very small.

### 6d CWH with Gaussian Disturbance
We use a 6d CWH system with additive Gaussian noise to model the dynamics of three satellites performing rendezvous operations while maintaining minimum distances from eachother. This example compares the proposed method with particle control. Our method performed significantly faster and satisfied the required probabilistic thresholds.
