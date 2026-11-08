# Mathematical Derivation: Student-t Log-Posterior & Jacobian Correction

**Course:** MATH 3428 Bayesian Statistics  
**Lead:** CHEUNG YUK YING 21063768 (Framework & Sampler)  
**Deliverable:** Model Formulation & Derivation Draft  

---

## 1. Model Likelihood (Student-$t$)

We assume the daily financial log-returns $y = (y_1, y_2, \dots, y_n)^T$ follow an independent and identically distributed (i.i.d.) three-parameter Student-$t$ distribution:

$$y_i \overset{\text{iid}}{\sim} t_\nu(\mu, \sigma^2), \quad i = 1, \dots, n$$

The probability density function (PDF) of a single observation $y_i$ is given by:

$$p(y_i \mid \mu, \sigma, \nu) = \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\pi \nu}\,\sigma\,\Gamma\left(\frac{\nu}{2}\right)} \left[1 + \frac{1}{\nu}\left(\frac{y_i - \mu}{\sigma}\right)^2\right]^{-\frac{\nu+1}{2}}$$

Taking the product across all $n$ observations and computing the natural logarithm yields the log-likelihood function:

$$\ln p(y \mid \mu, \sigma, \nu) = n \left[ \ln\Gamma\left(\frac{\nu+1}{2}\right) - \ln\Gamma\left(\frac{\nu}{2}\right) - \frac{1}{2}\ln(\pi\nu) - \ln\sigma \right] - \frac{\nu+1}{2} \sum_{i=1}^n \ln\left[1 + \frac{1}{\nu}\left(\frac{y_i - \mu}{\sigma}\right)^2\right]$$

---

## 2. Parameter Transformation & Jacobian Adjustment

To ensure numerical stability and prevent illegal parameter proposals during Random Walk Metropolis sampling (specifically enforcing scale $\sigma > 0$ and degrees of freedom $\nu > 2$ for finite variance), we transform the constrained parameter vector $(\mu, \sigma, \nu)$ into an unconstrained space $\mathbb{R}^3$:

$$\theta = (\theta_1, \theta_2, \theta_3)^T \in \mathbb{R}^3$$

The forward and inverse transformations are defined as:
* $\theta_1 = \mu \iff \mu = \theta_1$
* $\theta_2 = \ln\sigma \iff \sigma = \exp(\theta_2)$
* $\theta_3 = \ln(\nu - 2) \iff \nu = 2 + \exp(\theta_3)$
