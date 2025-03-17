# *****************************************************************************
#
# *****************************************************************************

import json
import os

from multiprocessing import Pool

from idmtools.core.id_file import read_id_file
from idmtools.core.platform_factory import Platform
from idmtools.assets import Asset
from idmtools_models.python.python_task import PythonTask
from idmtools_platform_comps.ssmt_work_items.comps_workitems \
                                        import SSMTWorkItem

from COMPS import Client
from COMPS.Data import Experiment
from COMPS.Data.Simulation import SimulationState

from py_assets_common.emod_constants import DOCK_PACK, COMPS_URL, D_FILE, \
                                            COMPS_ID_FILE, LOCAL_EXP_DIR, \
                                            O_FILE

# *****************************************************************************


def getter_worker(sim):

    Client.login(COMPS_URL)

    sim_obj = sim.retrieve_output_file_info(None)
    ret_val = None

    for file_info in sim_obj:
        if (file_info.friendly_name == O_FILE):
            file_bytes = sim.retrieve_output_files_from_info([file_info])
            ret_val = file_bytes[0].decode()

    return ret_val

# *****************************************************************************


# Runs on cluster
def pool_manager(exp_id=None):

    Client.login(COMPS_URL)

    if (not exp_id):
        (exp_id, _, _, _) = read_id_file(os.path.join('Assets', COMPS_ID_FILE))

    exp_obj = Experiment.get(exp_id)
    sims_all = exp_obj.get_simulations()
    sims_valid = [s for s in sims_all if s.state.value >=
                  SimulationState.Commissioned.value]

    with Pool() as pool_obj:
        resp_list = pool_obj.map(getter_worker, sims_valid)

    merged_dict = dict()
    for resp in resp_list:
        if (resp):
            merged_dict.update(json.loads(resp))

    with open(D_FILE, 'w') as fid01:
        json.dump(merged_dict, fid01)

    return None

# *****************************************************************************


# Runs locally
def get_sim_files(exp_id=''):

    # Connect to COMPS
    plat = Platform(type='COMPS', endpoint=COMPS_URL, environment='Calculon')

    # Create python task for SSMT work item
    f_path = os.path.abspath(__file__)
    f_dir = os.path.dirname(f_path)
    f_name = os.path.basename(f_path)
    task_obj = PythonTask(python_path='python3', script_path=f_name)

    # Add script for python task and exp id file to assets
    asset01 = Asset(filename=COMPS_ID_FILE, content=exp_id)
    asset02 = Asset(filename=f_path)
    task_obj.common_assets.add_asset(asset01)
    task_obj.common_assets.add_asset(asset02)

    # Add everything in the common python scripts directory as assets
    PATH_COMMON = os.path.join(f_dir, 'py_assets_common')
    task_obj.common_assets.add_directory(PATH_COMMON,
                                         relative_path='py_assets_common')

    # Reduce experiment output to single file
    wi_obj = SSMTWorkItem(name='ReduceExpOutput', task=task_obj,
                          docker_image=DOCK_PACK)
    wi_obj.run(wait_until_done=True)

    # Download reduced output and delete work item
    resp_dict = plat.get_files(wi_obj, [D_FILE])
    ret_val = json.loads(resp_dict[D_FILE].decode())
    plat_obj = wi_obj.get_platform_object()
    plat_obj.delete()

    return ret_val

# *****************************************************************************


def get_data_brick(run_local=False):

    if (run_local):
        with open(LOCAL_EXP_DIR) as fid01:
            exp_dir = fid01.readlines()[0].strip()

        data_brick = dict()
        for dir_item in os.listdir(exp_dir):
            fpath = os.path.join(exp_dir, dir_item, O_FILE)
            if (os.path.isfile(fpath)):
                with open(fpath) as fid01:
                    data_bit = json.load(fid01)
                data_brick.update(data_bit)
    else:
        # Get Experiment ID
        (exp_id, _, _, _) = read_id_file(COMPS_ID_FILE)

        # Reduce output and write data brick
        data_brick = get_sim_files(exp_id)

    with open(D_FILE, 'w') as fid01:
        json.dump(data_brick, fid01)

    return None


# *****************************************************************************


if (__name__ == "__main__"):

    pool_manager()

# *****************************************************************************
