# Catalytic-Cox-Regrerssion

In this repo, we implement the methods we used in the paper (reference). 

## Point estimation
By sourceing the file "point_estimation_cat.R", we have function for computing following two estimators:
### Estimation with Posterior mode
```math
\hat{\boldsymbol{\beta}}_{CR, \tau}=\arg \max _{\boldsymbol{\beta}}\left\{\log P L(\boldsymbol{\beta})+\frac{\tau}{M} \log L\left(\boldsymbol{\beta}, h_0^{+} \mid\left\{\left(\boldsymbol{X}_i^*, Y_i^*, \delta_i^*\right)\right\}_{i=1}^M\right)\right\}        
```
The function **betahat_RL_estimation(...)** will return the $\hat{\boldsymbol{\beta}}_{R L, \tau}$
### Estimation utilizing the synthetic data

```math
P \tilde{L}(\boldsymbol{\beta})=\left[\prod_{i \in \boldsymbol{I}_D}\left(\frac{\exp \left(\boldsymbol{X}_i^{\top} \boldsymbol{\beta}\right)}{\sum_{j \in \tilde{\mathcal{R}}_i \cap \boldsymbol{I}_D} \theta_j+\frac{\tau}{M} \sum_{j \in \tilde{\mathcal{R}}_i \cap \boldsymbol{I}_{D^*}} \theta_j^*}\right)^{\delta_i}\right]\left[\prod_{i \in \boldsymbol{I}_{D^*}} \frac{\exp \left(\boldsymbol{X}_{i-n}^*{ }^{\top} \boldsymbol{\beta}\right)}{\sum_{j \in \tilde{\mathcal{R}}_i \cap \boldsymbol{I}_D} \theta_j+\frac{\tau}{M} \sum_{j \in \tilde{\mathcal{R}}_i \cap \boldsymbol{I}_{D^*}} \theta_j^*}\right]^{\tau / M}
```
$\hat{\boldsymbol{\beta}}_{\text {WM, } \tau}$ that maximize above $P \tilde{L}(\boldsymbol{\beta})$ can be computed via **betahat_mix_estimation(...)**.

## Bayesian analysis
The Bayesian sampling can be performed via function **Bayesian_cox_cat_stan(...)**, which is defined in "Bayesian_cox_cat.R".

## Demonstration
See file "Demo.R"
