# Polio

Demonstration simulations for poliovirus type-2 transmission in Nigeria. The model structure was used to examine the impact Sabin vaccine reversion dynamics in contrast to nOPV vaccines.

Important features include:

- Multi-node network used to represent the country of Nigeria at either the
  admin-2 level or the 10km length scale.
- Network infectivity contagion transfer between nodes.
- Acquisition-transmission covariance.
- HINT structure to create non-participatory populations.
- Genome mutation dynamics and mutation labeling.
- Genome based variable infectivity.
- Clade labeling to create variable genome reversion trajectories.
- Immunity initialization from external data (IDM's Polio Immunity Mapper).
- Pre-specificed calendar of supplemental immunization activities (SIAs).
- Maternally derived immunity.
- Node-level variance in R0 value.
- Overdispersion of the infection rate.
- Maximum simulation duration based on elapsed time.
- Weighted agents to represent multiple individuals per agent.
- Population serosurveys implemented using custom reporters.
