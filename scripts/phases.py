import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.fft import irfft, rfft, rfftfreq

import matplotlib
import matplotlib.pyplot as plt

SKIP = 200
LOG_CUTOFF = -6
N_STDS = 2

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

def cut_devs(ext_ts, ext_vs, l_max_t, l_max_v, keep_n=-1):
    '''estimate the standard deviation of a distribution by taking the full width at half maximum and then return this value along with copies of t_pts and a_pts truncated to only include points within keep_n deviations from the center.
    ext_ts: x axis values of the distribution
    ext_vs: pdf
    keep_n: the number of standard deviations (estimated by fwhms) to keep, use a negative value to not perform any truncation
    returns: a tuple of three values (fwhm, ext_tsated, ext_vsated)
    '''
    fwhm = 0
    l_max_i = 0
    for i, v in enumerate(ext_vs):
        if fwhm == 0 and v > l_max_v/2:
            fwhm = l_max_t-ext_ts[i]
        elif fwhm != 0 and v < l_max_v/2:
            fwhm = (fwhm + ext_ts[i]-l_max_t)/2
            break
        if v == l_max_v:
            l_max_i = i
    #check if a truncation must be performed
    if keep_n < 0:
        return fwhm, ext_ts, ext_vs
    #do it
    i_min = -1
    for i, t in enumerate(ext_ts):
        if i_min < 0:
            if t > l_max_t - fwhm*keep_n:
                i_min = i
        elif t > l_max_t + fwhm*keep_n:
            return fwhm, ext_ts[i_min:i], ext_vs[i_min:i]
    #make sure that something is returned
    if i_min < 0:
        return fwhm, ext_ts, ext_vs
    else:
        return fwhm, ext_ts[i_min:], ext_vs[i_min:]

#return a Gaussian envelope over the specified time points with an amplitude amp, a mean t_0 and a standard deviation sqrt(width_sq/2)
def gauss_env(x, t_pts):
    return x[0]*np.exp(-(t_pts - x[1])**2/x[2])

#return a Gaussian pulse series with parameters specified by x (amplitude, mean, 2*std_dev**2, phase)
def gauss_series(x, t_pts):
    time_difs = t_pts - x[1]
    return x[0]*np.exp(-time_difs**2/x[2])*np.cos(x[3]*time_difs + x[4])

def opt_envelope(t_pts, a_pts, keep_n=-1):
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
    l_max_i = 0
    #typ_spacing is used later for estimating omega
    typ_spacing = 0.0
    for i in range(1, len(a_pts)-1):
        #extrema
        if (a_pts[i] > a_pts[i-1] and a_pts[i] > a_pts[i+1]) or (a_pts[i] < a_pts[i-1] and a_pts[i] < a_pts[i+1]):
            ext_ts.append(t_pts[i])
            ext_vs.append(abs(a_pts[i]))
            #check if this is a global maxima
            if ext_vs[-1] > l_max_v:
                l_max_i = i
                l_max_t = t_pts[i]
                l_max_v = ext_vs[-1]
                #set the typical spacing by the distance to the nearest earlier peak
                if len(ext_ts) > 1:
                    typ_spacing = ext_ts[-1] - ext_ts[-2]
    #remove everything that is below a constant cutoff threshold LOG_CUTOFF
    i_min = -1
    cut = l_max_v*np.exp(LOG_CUTOFF)
    for i, v in enumerate(ext_vs):
        if i_min < 0:
            if v > cut:
                i_min = i
        elif v < cut and i >= l_max_i:
            #make sure we include at least this data point
            ext_ts = ext_ts[i_min:i]
            ext_vs = ext_vs[i_min:i]
            break
    #make extra sure there aren't any divisions by zero
    if typ_spacing == 0.0:
        #this means there isn't enough data to perform an estimate of the spacing, just use omega_0
        if len(ext_ts) < 2:
            typ_spacing = np.pi/omega_0
        else:
            typ_spacing = ext_ts[1] - ext_ts[1]
    #use the full width at half max as a crude estimate of standard deviation
    fwhm, ext_ts, ext_vs = cut_devs(ext_ts, ext_vs, l_max_t, l_max_v, keep_n=keep_n)
    #now take the logarithm to transform fitting a gaussian to fitting a parabola. We normalize so that the amplitude should be roughly 1.0. This is easily accounted for at the end
    ext_ts = np.array(ext_ts)
    ext_vs = np.log(np.array(ext_vs)/l_max_v)
    #least squares of the logs of the envelope points
    def ff(x):
        return np.sum( (np.log(abs(x[0])) - (ext_ts - x[1])**2/x[2] - ext_vs)**2 )
    #gradient of the least squares of the envelope points
    def jac_env(x):
        diff_ser = (np.log(abs(x[0])) - (ext_ts - x[1])**2/x[2] - ext_vs)
        t_dif_by_w = (ext_ts - x[1])
        ret = np.zeros(3)
        ret[0] = 2*np.sum( diff_ser/abs(x[0]) )
        ret[1] = 4*np.sum( diff_ser*t_dif_by_w )/x[2]
        ret[2] = 2*np.sum( diff_ser*(t_dif_by_w**2) )/(x[2]**2)
        return ret

    #note that we fit to a form e**(-(t-t_0)^2/w) so using fwhm as an estimate for sigma we have w=2*sigma**2
    x0 = np.array([1.0, l_max_t, 2*fwhm*fwhm])
    res = opt.minimize(ff, x0, jac=jac_env)

    #we normalized all of our points by the global maxima, so we need to correct that
    res.x[0] = np.abs(res.x[0])*l_max_v
    #now create crude estimates of omega and phi by looking at the spacing between peaks and the phase offset between the data peak and the estimate for t_0
    est_omega = np.pi/typ_spacing
    est_phi = fix_angle(est_omega*(res.x[1] - l_max_t))

    return res, est_omega, est_phi

def opt_pulse_env(t_pts, a_pts, keep_n=-1):
    '''Set values for the a_pts for each time point
    a_pts: the values of the electric field at each corresponding time point in t_pts. Note that len(a_pts) must be the same as len(t_pts)'''
    if a_pts.shape != t_pts.shape:
        raise ValueError("t_pts and a_pts must have the same shape")

    #fit the pulse envelope, this gives us a good estimate of three points
    env_res, est_omega, est_phi = opt_envelope(t_pts[::SKIP], a_pts[::SKIP], keep_n=keep_n)
    #truncate so that we only look at N_STDS deviations from the center of the pulse
    if keep_n > 0:
        dt = (t_pts[-1]-t_pts[0])/t_pts.shape[0]
        i_0 = int(env_res.x[1]/dt)
        i_shift = int(keep_n*np.sqrt(env_res.x[2]/2) / dt)
        #make sure we don't go past the ends of the array
        if i_0 < i_shift:
            i_shift = int(dt*env_res.x[1])
        if i_0+i_shift > t_pts.shape[0]:
            i_shift = int(t_pts.shape[0] - dt*env_res.x[1])
        t_pts = t_pts[i_0-i_shift:i_0+i_shift]
        a_pts = a_pts[i_0-i_shift:i_0+i_shift]

    #it is helpful to try to get the amplitude close to unity so that the gradient in other directions is not suppressed
    if env_res.x[0] != 0.0:
        a_pts = a_pts / env_res.x[0]

    #square error of the pulse envelope
    def fp(x):
        return np.sum( (gauss_series(x, t_pts) - a_pts)**2 )
    def jac_fp(x):
        p_i = gauss_series(x, t_pts) - a_pts
        exp_cos = np.exp(-(t_pts-x[1])**2/x[2])*np.cos(x[3]*(t_pts-x[1])+x[4])
        exp_sin = np.exp(-(t_pts-x[1])**2/x[2])*np.sin(x[3]*(t_pts-x[1])+x[4])       
        t_dif_by_w = (t_pts-x[1])
        ret = np.zeros(5)
        ret[0] = 2*np.sum(p_i*exp_cos)
        ret[1] = 2*x[0]*np.sum(p_i*(2*t_dif_by_w*exp_cos/x[2] + exp_sin))
        ret[2] = 2*x[0]*np.sum(p_i*exp_cos*t_dif_by_w**2)/(x[2]**2)
        ret[3] = -2*x[0]*np.sum(p_i*exp_sin*t_dif_by_w)
        ret[4] = -2*x[0]*np.sum(p_i*exp_sin)
        return ret

    x0 = np.array([1.0, env_res.x[1], env_res.x[2], est_omega, est_phi])
    #super hacky way to account for pi phase shifts
    fp_old = fp(x0)
    x0[4] += np.pi
    if fp(x0) > fp_old:
        x0[4] = est_phi
    else:
        est_phi += np.pi

    #now we actually perform the minimization
    res = opt.minimize(fp, x0, jac=jac_fp)
    res.x[0] = res.x[0]*env_res.x[0]
    res.x[4] = fix_angle(res.x[4])

    return res, env_res, est_omega, est_phi

def get_params(res):
    '''Convert a scipy optimize result returned from opt_pulse_env into a human readable set of parameters
    '''
    amp = res.x[0]
    t_0 = res.x[1]
    sig = np.sqrt(res.x[2]/2)
    omega = res.x[2]
    cep = res.x[4]
    if amp < 0:
        amp *= -1
        if cep > 0:
            cep -= np.pi
        else:
            cep += np.pi
    cep = fix_angle(cep)
    return amp, t_0, sig, omega, cep
