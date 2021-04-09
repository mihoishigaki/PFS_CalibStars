import numpy as np

import sys

import h5py


def write_to_hdf5(filelist, datadict, outsimnameh5):


    
    f = h5py.File(outsimnameh5, "w")

    for file in filelist:
        id = (file.split(".")[-2])

        wv, fl, err = np.loadtxt(file, usecols = (0, 1, 2), unpack = True)

        grp = f.create_group(id)
        grp.create_dataset("wave", data = wv, dtype = 'f')
        grp.create_dataset("flux", data = fl, dtype = 'f')
        grp.create_dataset("error", data = err, dtype = 'f')

        for key in datadict.keys():
            grp.attrs[key] = datadict[key]
            
        
    f.close()

    return()



def read_h5(h5file):


    with h5py.File(h5file, 'r') as f:

        sn = f['0'].attrs['snr']
        
        waves = ()
        fluxs = ()
        errors = ()

        nobj = len(f.keys())
        
        for key in f.keys():

            wave = f[key]['wave'][:]
            flux = f[key]['flux'][:]
            error = f[key]['error'][:]

            waves = np.append(waves, wave)
            fluxs = np.append(fluxs, flux)
            errors = np.append(errors, error)
            


    
    waves = waves.reshape(nobj, len(wave))
    fluxs = fluxs.reshape(nobj, len(wave))
    errors = errors.reshape(nobj, len(wave))
    
    return(sn, waves, fluxs, errors)

