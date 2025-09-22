# Wiener Process Degradation Models (WPDM)

This repository runs Wiener process degradation models (WPDMs) using first difference estimators. These functions were used in the manuscript "Assessment of Wiener Process Degradation Models with Application to Amyotrophic Lateral Sclerosis Decline". The user can run this class of models with:

1.  Optimization of the full joint likelihood

2.  Profile likelihood methods with maximum likelihood or empirically unbiased estimators

As an implementation example, please refer to *sample code.R*. For additional details, review the documentation in the primary function WPDM.R.

Note, simulation data from the manuscript is organized as:

1.  sim1: $\mu$=-1.00, $\sigma^2 = 0.25$, $\sigma_{\mu}^2=0.25$, $\sigma_{\varepsilon}^2 = 0.25$
2.  sim2: $\mu$=-1.00, $\sigma^2 = 0.25$, $\sigma_{\mu}^2=0.05$, $\sigma_{\varepsilon}^2 = 0.49$
3.  param_checks**:** additional simulations in the Supporting Information
4.  nfu_checks:simulations relating to Figure S2 in the Supporting Information
5.  followup_time_checks: simulations relating to Figure 1 focusing on clustering of follow-up times.
