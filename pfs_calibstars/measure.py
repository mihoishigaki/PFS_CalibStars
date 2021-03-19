import numpy as np


import matplotlib.pyplot as plt
import scipy.optimize as so
import sys, logging

#pfs_calibstars_dir = '/Users/ishigakimiho/PFS/Github/pfs_calibstars'

#sys.path.append(pfs_calibstars_dir)
#from pfs_calibstars import calibstars as cs






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

    print(params)
    print(covar)
    if len(covar)==1:
        params_errs = np.zeros_like(params)
    else:
        params_errs = np.sqrt(covar.diagonal())
        
    dof = len(x) - 1 - len(params0)
    chi2 = np.sum(np.power(calc_chi(params, model_func, x, y, yerr), 2))

    return([params, params_errs, chi2, dof])



def corrfunc(x, params):

    # A simple gaussian function is assumed

    a, mu, sig, bg = params

    y = a * np.exp(-0.5 * ((x - mu)/sig)**2.) + bg

    return(y)
    



def cross_correlate_template(wobs, fobs, ferrobs, wsyn, fsyn, \
                             normalize = True, lower_velocity_limit = -400., \
                             upper_velocity_limit = 400., \
                             velocity_step = 0.5):

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
    
    plt.plot(wsyn, depth)
    plt.plot(wobs, depth_obs)
    plt.show()

    
    for i, vel in enumerate(velocity):
        factor = np.sqrt((1. - vel/c)/(1. + vel/c))
        shifted_template = np.interp(wobs, wsyn/factor, depth, \
                                     left = 0.0, right = 0.0)

        ccf[i] = np.correlate(depth_obs, shifted_template)[0]
        ccf_err[i] = np.correlate(ferrobs, shifted_template)[0]
        current_work_progress = ((i*1.0)/num_shifts) * 100

        if cs.report_progress(current_work_progress, last_reported_progress):
            last_reported_progress = current_work_progress
            logging.info("%.2f%%" % current_work_progress)



    max_ccf = np.max(ccf)
    ccf = ccf/max_ccf # Normalize
    ccf_err = ccf_err/max_ccf # Propagate errors

    plt.plot(velocity, ccf)
    plt.xlabel("RV [km/s]")
    plt.plot([100., 100.], [0., 2.0], linestyle = ":", color = "gray")
    plt.xlim(-100, 300.)
    plt.ylim(0.10, 1.10)


    
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

    # get results 
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
    plt.text(-95, 0.95, label)
    plt.show()

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
    
    return(velocity, ccf, ccf_err)

