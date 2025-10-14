# Measles - GHA

This model explores outbreak response strategies for measles, with Ghana as an example near-elimination context. Routine coverage of MCV1 in Ghana has been around 90% for more than a decade. Coverage at this level results in occasional interruption but is not high enough to prevent periodic resurgences or outbreaks from re-importation.

UN WPP estimates for total population and population age structure were used to create a representation of the demography of Ghana.

![Figure 1: Population pyramids](figures/ref_pyr_measles_gha01.png)

WorldPop estimates of the spatial distribution of the population were aggregated from constrained 1km approximations to the 10km scale, with individual nodes respecting district boundaries (second administrative level). Districts smaller than 100km^2^ were not subdivided.

![Figure 2: Population density and administrative boundaries](figures/ref_gha_map_pop01.png)

Measles infectivity moved between nodes based on a gravity-type network model. A spatially uniform rate of MCV1 was assumed for each region (first administrative level) based on time varying estimates of MCV1 coverage from IHME.

Each agent represented 25 individuals. Prevalence could fall to zero during the simulation, but exogeneous importation ensured long-term elimination did not occur. Importation probability was proportional to population. Simulations started in 2008. Outputs prior to 2011 were ignored because initial measles immunity is only roughly approximate and a few years of simulated time are required to establish internal consistency. Each simulation time step represents 5 days.

Parameters in the baseline model were adjusted to fit observed timeseries of measles incidence. A Poisson-based likelihood function was maximized over one free parameter that scales total incidence and was interpreted as a reporting rate.

![Figure 3: Baseline timeseries and reporting rate](figures/ref_gha_timeseries.png)

Posterior estimates suggest around 4% of measles infections are captured in case reports.

Test scenarios implemented hypothetical regional outbreak response SIAs starting in 2022. These responses were triggered based on varying thresholds of observed cases. Spatially heterogeneous reporting rates were assumed with a mean reporting rate of 4%

![Figure 4: Baseline timeseries and reporting rate](figures/ref_gha_obr_contours.png)

Increasing the surveillance rate in districts with below average reporting reduced the expected number of cumulative infections over the period Jan 2022 - Dec 2024. This effect was present at all response thresholds but was most notable for higher thresholds when outbreak response was less common.
