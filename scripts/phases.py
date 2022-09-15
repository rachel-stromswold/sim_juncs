import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.fft import irfft, rfft, rfftfreq

import matplotlib
import matplotlib.pyplot as plt

SKIP = 200
LOG_CUTOFF = -4

#returns the angle theta transformed to be in the range (-pi,pi]
def fix_angle(theta):
    if theta > np.pi:
        theta -= 2*np.pi*np.floor(theta/(2*np.pi))
        if theta > np.pi:
            theta -= 2*np.pi
    elif theta <= -np.pi:
        theta += 2*np.pi*(np.floor(-theta/(2*np.pi)) + 1)
        if theta > np.pi:
            theta -= 2*np.pi
    return theta

    def __init__(self, t_pts):
        '''A class that uses least squares fitting to a Gaussian pulse envelope to estimate phases
        t_pts: a numpy array of time points to be used for the fit'''
        self.t_pts = t_pts
        self.a_pts = np.zeros(t_pts.shape[0])

#return a Gaussian envelope over the specified time points with an amplitude amp, a mean t_0 and a standard deviation sqrt(width_sq/2)
def gauss_env(x, t_pts):
    return x[0]*np.exp(-(t_pts - x[1])**2/x[2])

#return a Gaussian pulse series with parameters specified by x (amplitude, mean, 2*std_dev**2, phase)
def gauss_series(x, t_pts):
    time_difs = t_pts - x[1]
    return x[0]*np.exp(-time_difs**2/x[2])*np.cos(x[3]*time_difs + x[4])

def opt_envelope(t_pts, a_pts):
    '''Performs a least squares fit to find a rough estimate of the shape of the pulse envelope. This is used by the opt_pulse to optimize the pulse shape inside the envelope
    t_pts: a numpy array of times at which each point occurred
    a_pts: the values corresponding to the t_pts
    returns: a tuple containing (np.array([times of maxima, maxima]), global time maximum, global maximum, the spacing between adjacent peaks)
    '''
    #find the extrema
    ext_ts = []
    ext_vs = []
    l_max_t = 0
    l_max_v = 0
    #typ_spacing is used later for estimating omega
    typ_spacing = 0.0
    for i in range(1, len(a_pts)-1):
        #extrema
        if (a_pts[i] > a_pts[i-1] and a_pts[i] > a_pts[i+1]) or (a_pts[i] < a_pts[i-1] and a_pts[i] < a_pts[i+1]):
            ext_ts.append(t_pts[i])
            ext_vs.append(abs(a_pts[i]))
            #check if this is a global maxima
            if ext_vs[-1] > l_max_v:
                l_max_t = t_pts[i]
                l_max_v = ext_vs[-1]
                #set the typical spacing by the distance to the nearest earlier peak
                if len(ext_ts) > 1:
                    typ_spacing = ext_ts[-1] - ext_ts[-2]
    #make extra sure there aren't any divisions by zero
    if typ_spacing == 0.0:
        #this means there isn't enough data to perform an estimate of the spacing, just use omega_0
        if len(ext_ts) < 2:
            typ_spacing = np.pi/omega_0
        else:
            typ_spacing = ext_ts[1] - ext_ts[1]
    #use the full width at half max as a crude estimate of standard deviation
    fwhm = 0
    for i, v in enumerate(ext_vs):
        if fwhm == 0 and v > l_max_v/2:
            fwhm = l_max_t-ext_ts[i]
        elif fwhm != 0 and v < l_max_v/2:
            fwhm = (fwhm + ext_ts[i]-l_max_t)/2
            break
    #now take the logarithm to transform fitting a gaussian to fitting a parabola. We normalize so that the amplitude should be roughly 1.0. This is easily accounted for at the end
    ext_ts = np.array(ext_ts)
    ext_vs = np.log(np.array(ext_vs)/l_max_v)
    #remove everything that is below a constant cutoff threshold LOG_CUTOFF
    i_min = -1
    for i, v in enumerate(ext_vs):
        if i_min < 0:
            if v > LOG_CUTOFF:
                i_min = i
        elif v < LOG_CUTOFF:
            #make sure we include at least this data point
            ext_ts = ext_ts[i_min:i]
            ext_vs = ext_vs[i_min:i]
            break
    #least squares of the logs of the envelope points
    def ff(x):
        return np.sum( (np.log(abs(x[0])) - (ext_ts - x[1])**2/x[2] - ext_vs)**2 )
    #gradient of the least squares of the envelope points
    def jac_env(x):
        diff_ser = (np.log(abs(x[0])) - (ext_ts - x[1])**2/x[2] - ext_vs)
        ret = np.zeros(3)
        ret[0] = 2*np.sum( diff_ser/abs(x[0]) )
        ret[1] = 4*np.sum( diff_ser*(ext_ts - x[1])/x[2] )
        ret[2] = 2*np.sum( diff_ser*((ext_ts - x[1])/x[2])**2 )
        return ret

    #note that we fit to a form e**(-(t-t_0)^2/w) so using fwhm as an estimate for sigma we have w=2*sigma**2
    x0 = np.array([1.0, l_max_t, 2*fwhm*fwhm])
    res = opt.minimize(ff, x0, jac=jac_env)

    #we normalized all of our points by the global maxima, so we need to correct that
    res.x[0] = np.abs(res.x[0])*l_max_v
    #now create crude estimates of omega and phi by looking at the spacing between peaks and the phase offset between the data peak and the estimate for t_0
    est_omega = np.pi/typ_spacing
    est_phi = fix_angle(est_omega*(l_max_t - res.x[1]))

    return res, est_omega, est_phi

def opt_pulse_env(t_pts, a_pts):
    '''Set values for the a_pts for each time point
    a_pts: the values of the electric field at each corresponding time point in t_pts. Note that len(a_pts) must be the same as len(t_pts)'''
    if a_pts.shape != t_pts.shape:
        raise ValueError("t_pts and a_pts must have the same shape")

    #square error of the pulse envelope
    def fp(x):
        return np.sum( (gauss_series(x, t_pts) - a_pts)**2 )
    def jac_fp(x):
        p_i = gauss_series(x, t_pts) - a_pts
        exp_cos = np.exp(-(t_pts-x[1])/x[2])*np.cos(x[3]*(t_pts-x[1])+x[4])
        exp_sin = np.exp(-(t_pts-x[1])/x[2])*np.sin(x[3]*(t_pts-x[1])+x[4])       
        t_dif_by_w = (t_pts-x[1])/x[2]
        ret = np.zeros(5)
        ret[0] = 2*np.sum(p_i*exp_cos)
        ret[1] = 2*x[0]*np.sum(p_i*(2*t_dif_by_w*exp_cos + exp_sin))
        ret[2] = 2*x[0]*np.sum(p_i*exp_cos*t_dif_by_w**2)
        ret[3] = -2*x[0]*x[2]*np.sum(p_i*exp_sin*t_dif_by_w)
        ret[4] = -2*x[0]*np.sum(p_i*exp_sin)
        return ret

    #fit the pulse envelope, this gives us a good estimate of three points
    env_res, est_omega, est_phi = opt_envelope(t_pts[::SKIP], a_pts[::SKIP])
    x0 = np.array([env_res.x[0], env_res.x[1], env_res.x[2], est_omega, est_phi])
    #super hacky way to account for pi phase shifts
    fp_old = fp(x0)
    x0[4] += np.pi
    if fp(x0) > fp_old:
        x0[4] = est_phi
    #now we actually perform the minimization
    res = opt.minimize(fp, x0, jac=jac_fp)
    res.x[4] = fix_angle(res.x[4])
    return res, env_res
