# Measles - GHA

Simulations examining RDT impact using Ghana as an example context. Work presented as part of Feb 2022 review.

Parameters in the baseline model were adjusted to fit observed timeseries of measles incidence. Poisson-based likelihood function is maximized over one free parameter that scales total incidence and is interpreted as a reporting rate.

Important features include:

- Multi-node network used to represent the country of Ghana at the 10km scale.
- Network infectivity contagion transfer between nodes.
- Spatially varying rates of routine immunization.
- Acquisition-transmission covariance.
- Pre-specified calendar of supplemental immunization activities (SIAs).
- Maternally derived immunity.
- Overdispersion of the infection rate.
- Maximum simulation duration based on elapsed time.
- Weighted agents to represent multiple individuals per agent.
- Age-based immunity initialization to approximate endemic transmission.
- Infectivity reservoir to represent persistent exogeneous importation.

Test scenario projects incidence forward and implements outbreak response SIAs based on observed incidence.

Important features include:

- SQL-based event reporting of symptomatic incidence.
- Spatial varying reporting rates sampled from a beta distribution.
- Response interventions created dynamically using python in-process.
