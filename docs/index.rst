=========================
Welcome to |EMOD_Generic|
=========================

|EMOD_Generic| is a collection of Python scripts and demonstration models for use with the generic simtype of |EMOD_s| (current with the Generic-Ongoing branch), |emod_api|, and |IT_s|. All of these models are for use on |COMPS_l| at |IDM_s| but can be easily adapted for use in other Slurm systems with support for SIF containers.

Each model is independent and implemented using the generic simtype of |EMOD_s|. For each model, one or more simulation ensembles (called experiments on |COMPS_s|) are provided that highlight interesting features and/or phenomena. |EMOD_s| simulates stochastic processes; individual outcomes may not be representative, so an ensemble of simulations is the natural basis for consideration.

Input files for each simulation are constructed using |emod_api| as a Python pre-processing step during execution on |COMPS_s|. Values to be varied between simulations or between experiments are considered simulation parameters and are conceptually distinct from software parameters. Simulation parameters are specified as arguments to the Python scripts that construct input files; those scripts are responsible for setting all necessary software parameters. Simulation parameters are typically a (small) subset of the software parameters required and there may not be a one-to-one mapping between simulation parameters and software parameters (e.g., a single simulation parameter named R0 may set several software parameters in order to ensure the model produces the desired R0 value). Model construction is the process of determining which software parameters are static and which software parameters are to be varied between simulations. The software parameter "Run_Number" is used to set the random number seed for a simulation and should always be included as a simulation parameter.

Additional information about how to use |IT_s| can be found in the :doc:`documentation <idmtools:index>`. Additional information about software parameters for the generic simtype of |EMOD_s| can be found in the :doc:`parameter overview <emod-generic:parameter-overview>`.

:doc:`Click here <idmtools:index>` for a diagram showing how |IT_s| and each of the related packages are used in an end-to-end workflow using |EMOD_s| as the disease transmission model.

Contents
========

.. toctree::
   :maxdepth: 3
   :titlesonly:

   get-started
   workflow
   examples
