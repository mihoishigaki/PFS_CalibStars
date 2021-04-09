import numpy as np


import matplotlib.pyplot as plt
import scipy.optimize as so
import sys, logging, glob

import types


import pfs_calibstars.utility as ut
import pfs_calibstars.simulate as sm

import pandas as pd

from astropy.stats import sigma_clip


def calc_chi(params, model_func, x, y, yerr):


    model = model_func(x, params)
    chi = (y - model) / yerr
    
    return(chi)


def solve_leastsq(x, y, yerr, params0, model_func):


    output = so.leastsq(
        calc_chi,
        params0,
        args = (model_func, x, y, yerr),
        full_output = True)

    params, covar, info, mesg, ier = output

    
    if covar is None:
        params = (np.zeros(len(params0))).fill(np.nan)
        params_errs = np.zeros_like(params0)
        chi2 = np.nan

    else:
        params_errs = np.sqrt(covar.diagonal())
        chi2 = np.sum(np.power(calc_chi(params, model_func, x, y, yerr), 2))
        

        
    dof = len(x) - 1 - len(params0)

    
    

    return([params, params_errs, chi2, dof])



def corrfunc(x, params):

    # A simple gaussian function is assumed

    a, mu, sig, bg = params

    y = a * np.exp(-0.5 * ((x - mu)/sig)**2.) + bg

    return(y)
    



def cross_correlate_template(wobs, fobs, ferrobs, wsyn, fsyn, \
                             plot = False, normalize = False,
                             lower_velocity_limit = -400., \
                             upper_velocity_limit = 400., \
                             velocity_step = 0.2):

    last_reported_progress = -1


    # speed of light in km/s
    c = 299792.458


    velocity = np.arange(lower_velocity_limit, \
                         upper_velocity_limit + velocity_step, \
                         velocity_step)
    shifts = np.int32(velocity / velocity_step)
    num_shifts = len(shifts)



    # Initialize a cross-correlation function
    ccf = np.zeros(num_shifts)
    ccf_err = np.zeros(num_shifts)

    if normalize == True:
        depth = np.abs(np.max(fsyn) - fsyn)
        depth_obs = np.abs(np.max(fobs) - fobs)
    else:
        depth = fsyn - 1.
        depth_obs = fobs - 1.


    if plot == True:
        plt.plot(wsyn, depth)
        plt.plot(wobs, depth_obs)
        plt.xlim(840., 890.)
        plt.ylim(-0.5, 0.5)
        plt.show()
        
    for i, vel in enumerate(velocity):
        
        factor = np.sqrt((1. - vel/c)/(1. + vel/c))
        shifted_template = np.interp(wobs, wsyn/factor, depth, \
                                     left = 0.0, right = 0.0)

        ccf[i] = np.correlate(depth_obs, shifted_template)[0]
        ccf_err[i] = np.correlate(ferrobs, shifted_template)[0]
        current_work_progress = ((i*1.0)/num_shifts) * 100

        #if cs.report_progress(current_work_progress, last_reported_progress):
            #last_reported_progress = current_work_progress
            #logging.info("%.2f%%" % current_work_progress)



    max_ccf = np.max(ccf)
    ccf = ccf/max_ccf # Normalize
    ccf_err = ccf_err/max_ccf # Propagate errors

    #plt.plot(velocity, ccf)
    #plt.xlabel("RV [km/s]")
    #plt.plot([100., 100.], [0., 2.0], linestyle = ":", color = "gray")
    #plt.xlim(-100, 300.)
    #plt.ylim(0.10, 1.10)


    
    # Interpolate synthetic spectra at observed wavelengths
    #i_func = interpolate.interp1d(wsyn, fsyn, kind = "cubic")
    #fsyn_i = i_func(w)

    #i_func = interpolate.interp1d(wobs, fobs, kind = "cubic")
    #fobs_i = i_func(w)
    
    
    # Fit the correlation with a gaussian function
    
    
    ## Initial guess of the parameters

    a0 = -1.0
    mu0 = velocity[np.argmax(ccf)]
    sig0 = 10.
    bg0 = 1.0
    params0 = np.array([a0, mu0, sig0, bg0])
    params, params_errs, chi2, dof = \
        solve_leastsq(velocity, ccf, ccf_err, params0, corrfunc)


    print(params, params_errs, chi2)
    
    # get results


    if params is None:
        mu = np.nan
        mu_err = np.nan

    else:

        mu = params[1]
        sig = np.abs(params[2])

        mu_err = params_errs[1]
        sig_err = np.abs(params_errs[2])

        fwhm = 2.35 * sig
        fwhme = 2.35 * sig_err

        label = "Mu = " + str("%4.2f(+/-%4.2f) [km/s]" % (mu, mu_err)) + \
            " Sigma = " + str("%4.2f(+/-%4.2f)" % (sig, sig_err) )
        print(label)
        label = "Mu = " + str(r"%4.1f $\pm$ %3.1f [km/s]" % (mu, mu_err)) 

        #plt.text(-95, 0.95, label)
        #plt.show()

        #model_y = corrfunc(lags, params)
        #plt.errorbar(lags, corr, yerr = corr_err, fmt = "bo")
        #plt.plot(lags, model_y, 'r-')
        #plt.xlim(mu - 5.*sig, mu + 5.*sig)
        #plt.show()


        # Calculate w_shift as a function of lags:
        
        #w_shift = ()
        #for k, ll in enumerate(lags):

        #    ww = w[k] - w[lags == 0]
        
        #    w_shift = np.append(w_shift, ww)

        #plt.plot(lags, w_shift)
        #plt.show()
        #i_func = interpolate.interp1d(velocity, w_shift, kind = "cubic")
        #w_shift_i = i_func(mu)
        #w_shift_i_n = i_func(mu - sig)
        #w_shift_i_p = i_func(mu + sig)
        #print(w_shift_i, w_shift_i_n, w_shift_i_p)
        #w_shift_sig = 0.5 * (np.abs(w_shift_i_n - w_shift_i) \
            #                     + np.abs(w_shift_i_p- w_shift_i))
        
    
        
        #wshift = wobs[lags == 0] - wobs[lags == lag]
        #z = w_shift_i / w[lags == lag]
        #v = z * c
        #verr = (w_shift_sig / w[lags == lag]) * c
        #print("Velocity: %.1f +- %.1f km/s"%(v, verr))
        
    return(mu, mu_err)





def measure_rv_mags(ObjIDs, objtypes, libID, data_dir, band, mags, outdir, \
                    rv0, cont_file, wcover_file):


    for k, ObjID in enumerate(ObjIDs):


        templates = glob.glob("../../pfs_calibstars_data/database/" + libID \
            + "/" + ObjID + "/Synspec/*_Vrad0.0.txt")


        #templatepath = "../../pfs_calibstars_data/database/" + libID \
        #    + "/" + ObjID \
        #                         + "/Simulated/MR/MP0.00_FA0.000_MTangle060_TP1.00/"
        
        #templates = glob.glob(templatepath + "*_" + ObjID + ".0.dat")
        
        if len(templates) == 0:
            print("template: ", templates, " not found")
            sys.exit()

        elif len(templates) > 1:
            print("more than one templates!")
            sys.exit()


        rverrs = ()
        sns = ()
            
        for mag in mags:

        
            h5files = glob.glob(data_dir + "/" + libID + "/" + \
                                ObjID + "/*/*/*" + band + \
                                "_%.1f_*.h5"%(mag))

            if len(h5files) == 0:
                print("h5 file not found!")
                sys.exit()
                
            elif len(h5files) >1:
                print("more than 1 h5 files.")
                sys.exit()

            sn, rv, rvsig = measure_rv(h5files[0], templates[0], rv0, \
                                   cont_file, wcover_file)

            rverrs = np.append(rverrs, rvsig)
            sns = np.append(sns, sn)


        data = {'gmag': mags, 'rvsig': rverrs, 'SN': sns}
        df = pd.DataFrame(data)
            
        outfile = outdir + "/mag_rv_sn_" + objtypes[k] + ".csv" 
        df.to_csv(outfile, index = False)

        
    return()


def measure_rv(h5file, template_spec, rv0, cont_file, wcover_file, plot = False):

    
    

    wsyn0, fsyn0 = np.loadtxt(template_spec, usecols = (0, 1), unpack = True)


    filt = filter_wl(wsyn0, cont_file, wcover_file, ['b', 'r', 'n'], \
                     cont_only = False)

    
    #filt = (wsyn0 > wmin - 10.) & (wsyn0 < wmax + 10.) 

    wsyn = wsyn0[filt]


    filt_cont = filter_wl(wsyn0, cont_file, wcover_file, ['b', 'r', 'n'], \
                          cont_only = True)
    wsyn_cont = wsyn0[filt_cont]
    fsyn_cont = fsyn0[filt_cont]
    
    
    fsyn, fsyn_e = sm.normalize_spec\
        (wsyn_cont, fsyn_cont, wsyn, fsyn0[filt], \
         np.zeros(len(fsyn0[filt])), medfilt_width = 20, \
         maxfilt_width = 20, order = 7)
    
    
    sn, waves, fluxes, errors = ut.read_h5(h5file)

    rvs = ()
    rverrs = ()
    for i in range(0, len(waves)):

        # Normalize each arm separately
        wave0 = ()
        flux = ()
        error = ()
        
        for arm in ['b', 'r', 'n']:
            
            filt = filter_wl(waves[i], cont_file, wcover_file, \
                                  arm, cont_only = False)
            wv = (waves[i])[filt]
            fl = (fluxes[i])[filt]
            err = (errors[i])[filt] 
             
            filt_cont = filter_wl(waves[i], cont_file, wcover_file, \
                             arm, cont_only = True)
            wl_cont = (waves[i])[filt_cont]
            fl_cont = (fluxes[i])[filt_cont]
            

            fl_norm, err_norm = sm.normalize_spec\
                (wl_cont, fl_cont, wv, fl, err, medfilt_width = 20, \
                 maxfilt_width = 20, order = 4)
            
            wave0 = np.append(wave0, wv)
            flux= np.append(flux, fl_norm)
            error = np.append(error, err_norm)
            


        # Apply a wavelength shift 
        wave = rvshift(wave0, -1.*rv0)

        
        #flux_n, error_n = sm.normalize_spec(wave, flux, wave, flux, error)
        
        rv, rverr = cross_correlate_template\
            (wave, flux, error, \
                                 wsyn, fsyn, plot = plot)

        rvs = np.append(rvs, rv)
        rverrs = np.append(rverrs, rverr)
        

    rvs = (np.array(rvs))[~np.isnan(np.array(rvs))]
    rvs_filtered = sigma_clip(rvs, sigma = 2, maxiters = 5)
    print("rvs_filtered", rvs_filtered)
    rv_mean = np.mean(rvs_filtered)
    rv_sig = np.std(rvs_filtered)

    return(sn, rv_mean, rv_sig)



def rvshift(w0, rv):

    # Speed of light in km/s
    c = 299792.4580
    w = w0 * np.sqrt((1.-(rv/c))/(1.+(rv/c)))

    return(w)



def filter_wl(w0, cont_file, wcover_file, arms, cont_only):

    # e.g., arms = ['b', 'r', 'n']
    # If cont_only == False, continuum region is not taken into account. 


    if np.isscalar(arms):
        arms = [arms]

    
    df = pd.read_csv(wcover_file)


    
    for i, arm in enumerate(arms):
    
        # Read wavelength coverage:

        wmin = (df['wmin'][df['id'] == arm].values)[0]
        wmax = (df['wmax'][df['id'] == arm].values)[0]


        if i == 0:
            filt_wc = (w0 > wmin) & (w0 < wmax)

        else:
            filt_wc = filt_wc | \
                    ((w0 > wmin) & (w0 < wmax)) 
            


    if cont_only == False:

        return(filt_wc)


    
        
    df = pd.read_csv(cont_file)

    
    for i in range(0, len(df.index)):

        wmin = df['wmin'][i]
        wmax = df['wmax'][i]


        if i == 0:
            filt_cont = (w0 > wmin) & (w0 < wmax)
        else:
        
            filt_cont = filt_cont | \
                ((w0 > wmin) & (w0 < wmax))

            
    filt = filt_wc & filt_cont 
        
    return(filt)
    
