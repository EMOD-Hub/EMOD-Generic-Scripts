# *****************************************************************************
#
# *****************************************************************************

import os

from idmtools.core.platform_factory import Platform
from idmtools.entities.command_task import CommandTask
from idmtools.entities.experiment import Experiment
from idmtools.assets import AssetCollection
from idmtools_platform_comps.utils.singularity_build \
                                         import SingularityBuildWorkItem
from idmtools_platform_comps.utils.assetize_output.assetize_output \
                                         import AssetizeOutput

from py_assets_common.emod_constants import ID_OS, ID_EXE, ID_ENV, ID_SCHEMA

# ****************************************************************************


# Retrieve base image from Docker hub; save locally to avoid repeated pulls
def make_OS_asset():

    # Name for the OS
    OS_NAME = os.getcwd().split('_')[1]

    # Prepare the platform
    plat_obj = Platform(type='COMPS',
                        endpoint='https://comps.idmod.org',
                        environment='Calculon')

    # Creates a work item to build image
    def_file_name = 'EMOD_OS_'+OS_NAME+'.def'
    sbwi_obj = SingularityBuildWorkItem(name='Build_EMOD_OS_'+OS_NAME,
                                        definition_file=def_file_name,
                                        force=True)

    # Wait until the image is built
    ac_obj = sbwi_obj.run(wait_until_done=True, platform=plat_obj)

    # Save asset id for sif to file
    ac_obj.to_id_file(ID_OS)

    return None

# ****************************************************************************


# Build executable and reporters; generate schema
def make_EXE_asset():

    # Name for the OS
    OS_NAME = os.getcwd().split('_')[1]

    # Prepare the platform
    plat_obj = Platform(type='COMPS',
                        endpoint='https://comps.idmod.org',
                        environment='Calculon',
                        priority='Highest',
                        simulation_root='$COMPS_PATH(USER)',
                        node_group='idm_48cores',
                        num_cores='1',
                        num_retries='0',
                        exclusive='False')

    # Add image for base OS
    os_image = AssetCollection.from_id_file(ID_OS)

    # Creates a work item to create the build image
    def_file_name = 'EMOD_EXE_'+OS_NAME+'.def'
    sbwi_obj = SingularityBuildWorkItem(name='Build_EMOD_EXE_'+OS_NAME,
                                        definition_file=def_file_name,
                                        force=True)
    sbwi_obj.assets.add_assets(os_image)

    # Wait until the build image is finished
    ac_obj = sbwi_obj.run(wait_until_done=True, platform=plat_obj)

    # Magic words
    s_exe = 'singularity exec Assets/EMOD_EXE_'+OS_NAME+'.sif '

    # Command to generate schema
    cmd_opts = ' --get-schema --schema-path schema.json -P /PyScripts'
    cmd_line = s_exe + '/outputs/Eradication' + cmd_opts

    # Create CommandTask
    task_obj = CommandTask(command=cmd_line)
    task_obj.common_assets.add_assets(ac_obj)

    # Wait until asset collection has been built; save to file
    ao_obj01 = AssetizeOutput(no_simulation_prefix=True)
    ao_obj01.from_items(Experiment.from_task(task_obj))
    ao_obj01.run(wait_until_done=True, platform=plat_obj)
    ao_obj01.asset_collection.to_id_file(ID_EXE)
    print(ao_obj01.asset_collection.uid.hex)

    # Command to copy executable and reporters
    cmd_line = s_exe + 'cp -r /outputs/. .'

    # Create CommandTask
    task_obj = CommandTask(command=cmd_line)
    task_obj.common_assets.add_assets(ac_obj)

    # Wait until asset collection has been built
    ao_obj02 = AssetizeOutput(no_simulation_prefix=True)
    ao_obj02.from_items(Experiment.from_task(task_obj))
    ao_obj02.run(wait_until_done=True, platform=plat_obj)
    ao_obj02.asset_collection.to_id_file(ID_SCHEMA)
    print(ao_obj02.asset_collection.uid.hex)

    return None

# ****************************************************************************


# Run work items on COMPS
def make_ENV_asset():

    # Name for the OS
    OS_NAME = os.getcwd().split('_')[1]

    # Prepare the platform
    plat_obj = Platform(type='COMPS',
                        endpoint='https://comps.idmod.org',
                        environment='Calculon')

    # Add image for base OS
    os_image = AssetCollection.from_id_file(ID_OS)

    # Creates a single work item to create the image
    def_file_name = 'EMOD_ENV_'+OS_NAME+'.def'
    sbwi_obj = SingularityBuildWorkItem(name='Build_EMOD_ENV_'+OS_NAME,
                                        definition_file=def_file_name,
                                        force=True)
    sbwi_obj.assets.add_assets(os_image)

    # Wait until the image is built
    ac_obj = sbwi_obj.run(wait_until_done=True, platform=plat_obj)

    # Save asset id for sif to file
    ac_obj.to_id_file(ID_ENV)

    return None

# ****************************************************************************
