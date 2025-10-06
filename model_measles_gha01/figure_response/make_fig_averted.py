# *****************************************************************************

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as cm

# Ought to go in emodpy
sys.path.append(os.path.abspath(os.path.join('..', '..', 'local_python')))
sys.path.append(os.path.abspath(os.path.join('..', 'Assets', 'python')))

from py_assets_common.emod_constants import NUM_SIMS, P_FILE, D_FILE, \
                                            EXP_C, EXP_V, BASE_YEAR

# *****************************************************************************

DIRNAMES = ['experiment_gha_test']

# *****************************************************************************


def make_fig():

    CCUT = -4e3

    for dirname in DIRNAMES:

        targfile = os.path.join('..', dirname, P_FILE)
        with open(targfile) as fid01:
            param_dict = json.load(fid01)

        targfile = os.path.join('..', dirname, D_FILE)
        with open(targfile) as fid01:
            dbrick = json.load(fid01)

        nsims = int(param_dict[NUM_SIMS])
        log10_rep = np.array(param_dict[EXP_V]['log10_min_reporting'])
        num_cases = np.array(param_dict[EXP_V]['adm01_case_threshold'])

        tvals = dbrick.pop('tstamps')

        inf_data = np.zeros((nsims, len(tvals)))
        scale_vec = np.zeros((nsims, 1))-1
        cal_vec = np.zeros(nsims)

        for sim_idx_str in dbrick:
            if (not sim_idx_str.isdigit()):
                continue
            sim_idx = int(sim_idx_str)
            sim_obj = dbrick[sim_idx_str]

            inf_data[sim_idx, :] = np.array(sim_obj['timeseries'])
            scale_vec[sim_idx, 0] = sim_obj['rep_rate']
            cal_vec[sim_idx] = sim_obj['cal_val']

        gidx = (scale_vec[:,0] >= 0)
        gidx = gidx & (cal_vec > CCUT)

        ntval = np.array(tvals)/365.0 + BASE_YEAR
        test_idx = np.argwhere(ntval>2020)[0][0]

        inf_data = inf_data[gidx, :]
        log10_rep = log10_rep[gidx]
        num_cases = num_cases[gidx]
        inf_data[:, :test_idx] = 0
        inf_data = np.cumsum(inf_data, axis=1)

        # Figure
        fig01 = plt.figure(figsize=(8, 5))
        axs01 = fig01.add_subplot(1, 1, 1, label=None)
        plt.sca(axs01)

        axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
        axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
        axs01.set_axisbelow(True)

        axs01.set_xlabel('Minimum Surveillance',fontsize=14)
        axs01.set_ylabel('Case Threshold for Response', fontsize=14)

        axs01.set_xscale('log')

        xmin = 1e-3
        xmax = 1e-1
        lxmin = np.log10(xmin)
        lxmax = np.log10(xmax)
        ymin = 0
        ymax = 1000
        axs01.set_ylim(ymin, ymax)

        xval = log10_rep
        yval = num_cases
        cval = inf_data[:, -1]
        gdx2 = np.argsort(cval)[::-1]

        # Binning
        nbins = 7
        dyval = (ymax-ymin)/(nbins-1)
        dxval = (lxmax-lxmin)/(nbins-1)
        (xvec,yvec) = np.meshgrid(np.linspace(lxmin-dxval/2, lxmax+dxval/2, nbins+1),
                                  np.linspace(ymin-dyval/2, ymax+dyval/2, nbins+1))
        zmat = np.zeros((nbins,nbins))

        for k1 in range(nbins):
            for k2 in range(nbins):
                gidx = (xval>xvec[k1,k2]) & \
                       (xval<xvec[k1,k2+1]) & \
                       (yval>yvec[k1,k2]) & \
                       (yval<yvec[k1+1,k2])
                zmat[k1,k2] = np.mean(cval[gidx])

        plt.contour(np.power(10.0,xvec[:-1,:-1]+dxval/2),
                                  yvec[:-1,:-1]+dyval/2,
                                  zmat, levels=[25e3,40e3,60e3,80e3], linewidths=3, vmin=0, vmax=100e3)

        axs01.scatter(np.power(10.0, xval[gdx2]), yval[gdx2], c=cval[gdx2], s=2, vmin=0, vmax=100e3, alpha=0.7)

        ticloc = [0.001, 0.01, 0.1]
        ticlab = ['0.1%', '1.0%', '10%']
        axs01.set_xticks(ticks=ticloc)
        axs01.set_xticklabels(ticlab)

        virmap      = plt.get_cmap('viridis')
        cbar_handle = plt.colorbar(cm.ScalarMappable(cmap=virmap), ax=axs01, shrink=0.75)

        axs01.plot([0.04, 0.04],[0,1000],'k:')

        ticloc = [0.0,0.2,0.4,0.6,0.8,1.0]
        ticlab = ['0','20k','40k','60k','80k','100k']

        cbar_handle.set_ticks(ticks=ticloc)
        cbar_handle.set_ticklabels(ticlab)
        cbar_handle.set_label('Cumulative Infections',fontsize=14,labelpad=10)
        cbar_handle.ax.tick_params(labelsize=14)

        plt.tight_layout()
        plt.savefig('fig_averted01.png')
        plt.close()

    return None

# *****************************************************************************


if (__name__ == "__main__"):

    make_fig()

# *****************************************************************************
