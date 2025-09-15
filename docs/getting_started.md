# Getting Started

1. Set up a virtual python environment.

2. Install requirements via `pip` using the IDM artifactory:

    `pip install . --extra-index-url=https://packages.idmod.org/api/pypi/pypi-production/simple`

3. Run an experiment:

    `cd EMOD-Generic-Scripts/model_covariance01/experiment_covariance01`

    `python make01_param_dict.py`

    `python make02_lauch_sims.py`

    `python make03_pool_brick.py`

4. Generate figures:

    `cd EMOD-Generic-Scripts/model_covariance01/figure_attackfrac01`

    `python make_fig_attackrate.py`
