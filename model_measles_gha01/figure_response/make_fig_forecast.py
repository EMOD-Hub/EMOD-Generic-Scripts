#*******************************************************************************

import os, sys, json

import numpy              as np
import matplotlib.pyplot  as plt
import matplotlib.patches as patch
import matplotlib         as mpl

#*******************************************************************************

DIRNAME = 'experiment_meas_gha_base01'
YMAX    =  350e3

targfile = os.path.join('..',DIRNAME,'param_dict.json')
with open(targfile) as fid01:
  param_dict = json.load(fid01)

targfile = os.path.join('..',DIRNAME,'data_brick.json')
with open(targfile) as fid01:
  data_brick = json.load(fid01)

nsims    = int(param_dict['NUM_SIMS'])
tvals    = data_brick.pop('tstamps')
ntstp    = len(tvals)
infBlock = np.zeros((nsims,ntstp))
sclVal   = np.zeros((nsims,1))-1
cdata    = np.zeros(nsims)

for sim_idx_str in data_brick:
  try:
    sim_idx = int(sim_idx_str)
  except:
    continue
  infBlock[sim_idx,:] = np.array(data_brick[sim_idx_str]['timeseries'])
  sclVal[sim_idx,0]   = data_brick[sim_idx_str]['rep_rate']
  cdata[sim_idx]      = data_brick[sim_idx_str]['cal_val']

fidx = (sclVal[:,0]>=0)

# Figure
fig01 = plt.figure(figsize=(8,6))

axs01  = fig01.add_subplot(111)
plt.sca(axs01)

axs01.grid(visible=True, which='major', ls='-', lw=0.5, label='')
axs01.grid(visible=True, which='minor', ls=':', lw=0.1)
axs01.set_axisbelow(True)

xval     = np.array(tvals)/365+1900

gidx     = (cdata > -4000) & fidx
test_idx = np.argwhere(xval>2020)[0][0]

infBlock              = infBlock[gidx,:]
infBlock[:,:test_idx] = 0
infBlock              = np.cumsum(infBlock,axis=1)

yval = np.mean(infBlock,axis=0)


infDatSetSort = np.sort(infBlock,axis=0)
infDatSetSort = infDatSetSort


for patwid in [0.475,0.375,0.25]:
  xydat = np.zeros((2*infDatSetSort.shape[1],2))
  xydat[:,0] = np.hstack((xval,xval[::-1]))
  tidx = int((0.5-patwid)*infDatSetSort.shape[0])
  xydat[:,1] = np.hstack((infDatSetSort[tidx,:],infDatSetSort[-tidx,::-1]))

  polyShp = patch.Polygon(xydat, facecolor='C0', alpha=0.7-patwid, edgecolor=None)
  axs01.add_patch(polyShp)

axs01.plot(xval,yval,color='C0',linewidth=2)
axs01.set_ylabel('Cumulative Infections - Projected',fontsize=16)
axs01.set_xlim(2020, 2026)
axs01.set_ylim(   0, YMAX)

ticlocx = [2020, 2021, 2022, 2023, 2024, 2025]
ticlabx = ['','','','','','']
axs01.set_xticks(ticks=ticlocx)
axs01.set_xticklabels(ticlabx)

ticlocy = [0, 50e3, 100e3, 150e3, 200e3, 250e3, 300e3, YMAX]
ticlaby = ['0','','100k','','200k','','300k','']
axs01.set_yticks(ticks=ticlocy)
axs01.set_yticklabels(ticlaby)

for k1 in ticlocx:
  axs01.text(k1+0.5,-0.04*YMAX,str(k1),fontsize=11,ha='center')

plt.tight_layout()
plt.savefig('fig_forecast_01.png')
plt.close()

#*******************************************************************************
