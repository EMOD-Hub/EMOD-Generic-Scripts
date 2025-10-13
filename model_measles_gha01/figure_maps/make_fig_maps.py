# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))

from py_assets_common.emod_local_proc import shape_patch, shape_line

# *****************************************************************************


def make_fig():

    dat_root = os.path.join('..', 'Assets', 'data')

    tpath = os.path.join(dat_root, 'shapes_gha_adm00.json')
    with open(tpath) as fid01:
        shp_adm00 = json.load(fid01)

    tpath = os.path.join(dat_root, 'shapes_gha_adm01.json')
    with open(tpath) as fid01:
        shp_adm01 = json.load(fid01)

    tpath = os.path.join(dat_root, 'shapes_gha_adm02.json')
    with open(tpath) as fid01:
        shp_adm02 = json.load(fid01)

    tpath = os.path.join(dat_root, 'shapes_gha_adm02_sub_100km.json')
    with open(tpath) as fid01:
        shp_adm02sub = json.load(fid01)

    tpath = os.path.join(dat_root, 'demog_data_ADM02_sub_100km.csv')
    with open(tpath) as fid01:
        pop_dat = [lval.strip().split(',') for lval in fid01.readlines()[1:]]
        dict_pop = {lval[0]: int(lval[3]) for lval in pop_dat}
        dict_area = {lval[0]: float(lval[4]) for lval in pop_dat}

    # Figure dimensions
    llpad = 0.2
    (xmin, xmax, ymin, ymax) = (180, -180, 90, -90)
    for cname in shp_adm00:
        pts = np.array(shp_adm00[cname]['points'])
        xmin = min(xmin, np.min(pts[:, 0])) - llpad
        xmax = max(xmax, np.max(pts[:, 0])) + llpad
        ymin = min(ymin, np.min(pts[:, 1])) - llpad
        ymax = max(ymax, np.max(pts[:, 1])) + llpad
        xyrat = (xmax-xmin)/(ymax-ymin)

    # Figure - Admin
    base_dim = 10
    cbar_add = 0.75
    fig01 = plt.figure(figsize=(base_dim, base_dim/xyrat))
    axs01 = fig01.add_subplot(1, 1, 1, label=None)
    plt.sca(axs01)

    axs01.axis('off')
    axs01.set_xlim(xmin, xmax)
    axs01.set_ylim(ymin, ymax)

    for cname in shp_adm00:
        adm_prt = shp_adm00[cname]['parts']
        adm_pts = shp_adm00[cname]['points']
        shape_patch(axs01, adm_pts, adm_prt)
        shape_line(axs01, adm_pts, adm_prt, wid=2.0)

    for cname in shp_adm01:
        adm_prt = shp_adm01[cname]['parts']
        adm_pts = shp_adm01[cname]['points']
        shape_line(axs01, adm_pts, adm_prt, wid=2.0)

    for cname in shp_adm02:
        adm_prt = shp_adm02[cname]['parts']
        adm_pts = shp_adm02[cname]['points']
        shape_line(axs01, adm_pts, adm_prt, wid=0.2)

    plt.tight_layout()
    plt.savefig('fig_admin01.png', transparent=False)
    plt.close()

    # Figure - Pop
    fig01 = plt.figure(figsize=(base_dim/cbar_add, base_dim/xyrat))
    axs01 = fig01.add_subplot(1, 1, 1, label=None)
    plt.sca(axs01)

    axs01.axis('off')
    axs01.set_xlim(xmin, xmax)
    axs01.set_ylim(ymin, ymax)

    for cname in shp_adm00:
        adm_prt = shp_adm00[cname]['parts']
        adm_pts = shp_adm00[cname]['points']
        shape_line(axs01, adm_pts, adm_prt, wid=2.0)

    for cname in shp_adm01:
        adm_prt = shp_adm01[cname]['parts']
        adm_pts = shp_adm01[cname]['points']
        shape_line(axs01, adm_pts, adm_prt, wid=2.0)

    for cname in shp_adm02:
        adm_prt = shp_adm02[cname]['parts']
        adm_pts = shp_adm02[cname]['points']
        shape_line(axs01, adm_pts, adm_prt, wid=0.2)

    virmap = mpl.colormaps['viridis']
    clr_map = mpl.cm.ScalarMappable(cmap=virmap)

    for cname in shp_adm02sub:
        adm_prt = shp_adm02sub[cname]['parts']
        adm_pts = shp_adm02sub[cname]['points']
        area_val = dict_area[cname]
        pop_val = dict_pop[cname]
        dense_val = pop_val/area_val + 1.0e-5

        lden_val = min(np.log(dense_val)/np.log(5000), 1.0)
        cval = virmap(lden_val)[:3]

        shape_patch(axs01, adm_pts, adm_prt, clr=cval)

    cbar_handle = plt.colorbar(mappable=clr_map, ax=axs01,
                               fraction=(1-cbar_add), shrink=cbar_add)

    ticloc = [0.1889, 0.4593, 0.72965, 1.0]
    ticlab = ['5', '50', '500', '5000']
    cbar_handle.set_ticks(ticks=ticloc)
    cbar_handle.set_ticklabels(ticlab)
    cbar_handle.set_label(r'Population Density (km$^{-2}$)',
                          fontsize=20, labelpad=10)
    cbar_handle.ax.tick_params(labelsize=16)

    plt.tight_layout()
    plt.savefig('fig_pop01.png', transparent=False)
    plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
