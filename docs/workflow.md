# Workflow Overview

The workflow used here separates defining parameters and uploading files (client-side operations) from writing inputs and running simulations (server-side).

## Client-side

On the client side, a Python environment with emodpy is used to communicate with COMPS.

1. Create a parameter dictionary that defines the set of simulations (experiment) by specifying values for all simulation parameters.

2. Upload the parameter dictionary along with the Python and data files used to construct model inputs. Run the simulations.

3. Collect outputs from COMPS.

## Server-side

The COMPS server supports a variety of environments via SIF containers. By default, EMOD will run using Debian and a Python environment with `emod-api`. 

On the COMPS server, each EMOD simulation automatically runs functions named `application` from three separate Python modules: `dtk_pre_process`, `dtk_in_process`, and `dtk_post_process`.

All EMOD input file creation is done as part of the function named `application` from `dtk_pre_process`. That function will find the ID assigned to the simulation, open the parameter dictionary and retrieve parameter values for the simulation, then use emod-api to create input files.

Every model directory containas an `Assets` folder with `python` and `data` subfolders. These subfolders are copied to COMPS and available for all siumulations in the experiment.

* `python`: Contains user-defined Python files used by EMOD. Some filenames must remain unchanged: `dtk_pre_process.py`, `dtk_in_process.py`, and `dtk_post_process.py` are assumed to be present and must contain a function named `application`. Everything else can be user-defined.
* `data`: Contains csv files, json files, etc. for use by the Python scripts when running EMOD. There are no required data files, but models without data files include an unused empty file (`no_data.json`) to ensure the directory structure is retained.
