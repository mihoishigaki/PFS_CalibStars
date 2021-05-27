import numpy as np

import sys, os

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



def create_zip_tar():

	import glob, subprocess
	
	files = glob.glob("../../pfs_calibstars_data/outputs/Synspec/*/MR/MP*/*pfsObject/", recursive = True)
	
	for file in files:
		tarfilename = file[:-1] + ".tar"
		command = "tar -cvf " + tarfilename + " " + file
		proc = subprocess.run(command, shell = True)
		command = "gzip " + tarfilename
		proc = subprocess.run(command, shell = True)
		print(file, tarfilename)
	#proc = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE)
	#output = proc.communicate()
	#dirlist = output.split('\n')
	#print(output[0].split('\n'))
	return()


#create_zip_tar()




