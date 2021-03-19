import numpy as np

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



