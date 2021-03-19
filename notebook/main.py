#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import imp
import os, glob,sys,pprint
import numpy as np
import matplotlib.pyplot as plt

import synphot
from astropy.utils.data import get_pkg_data_filename
from astropy.io import ascii,fits
import pandas as pd


from pfsspecsim import pfsetc
from pfsspecsim import pfsspec


# Please specify a place where the codes are located 
pfs_calibstars_dir = '/Users/ishigakimiho/PFS/Github/pfs_calibstars'
sys.path.append(pfs_calibstars_dir)
import pfs_calibstars as cs



# ## Create databases of simulated spectra 

# In[14]:


libID = "HDS"
cs.produce_database(libID)


# ## Produce many simulated spectra 

# In[18]:


libID = "HDS"
nreal = 100
output_rootdir = "../outputs"

ObjIDs = ["00092", "00081"]

band = "sdss_g"

for ObjID in ObjIDs:
    
    for mag in [17., 18., 19., 20., 21., 22., 23.]: 
        
        exptime, nexp = cs.get_exptime("sdss_g", mag, "science")
        nexps = [nexp]

        cs.simulate_many_spectra(ObjID, libID, nexps, band, mag, nreal,                                  output_rootdir, setting = "Optimistic",                                  write_h5 = True)


# ## Simulated spectra for the Mohammad's grid

# In[6]:


# Use inputs and outputs from etc.run()

sim = pfsspec.Pfsspec()
sim.set_param('outDir', outpath)
sim.set_param('etcFile', outsncname)
sim.set_param('MAG_FILE', magfilename)
sim.set_param('EXP_NUM',nexp)
sim.set_param('asciiTable', file[:-5] + ".sim")
sim.set_param('nrealize',1)
sim.make_sim_spec()

w, f = np.loadtxt(outpath + file[:-5] + ".sim.dat", usecols = (0, 1), unpack = True)
plt.plot(w, f)
plt.xlim(380., 510.)
plt.ylim(0.0, np.max(f[w < 500.]))
plt.show()




# In[ ]:




