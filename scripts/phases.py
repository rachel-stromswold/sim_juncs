import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.fft import irfft, rfft, rfftfreq

import matplotlib
import matplotlib.pyplot as plt

SKIP = 200
AMP_CUTOFF = 0.0025
N_STDS = 2
WIDTH_SCALE = 1000

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

def cut_devs(ext_ts, ext_vs, l_max_t, l_max_v, keep_n=-1, both_tails=False):
    '''estimate the standard deviation of a distribution by taking the full width at half maximum and then return this value along with copies of t_pts and a_pts truncated to only include points within keep_n deviations from the center.
    ext_ts: x axis values of the distribution
    ext_vs: pdf
    keep_n: the number of standard deviations (estimated by fwhms) to keep, use a negative value to not perform any truncation
    returns: a tuple of three values (fwhm, ext_tsated, ext_vsated)
    '''
    fwhm = 0
    l_max_i = 0
    if both_tails:
        last_i = 0
        for i, v in enumerate(ext_vs):
            #ensure that the fwhm is set before the maxima
            if v == l_max_v and fwhm == 0 and i > 0:
                fwhm = l_max_t-ext_ts[i-1]
            if fwhm == 0 and v > l_max_v/2:
                fwhm = l_max_t-ext_ts[i]
            elif fwhm != 0 and v > l_max_v/2:
                last_i = i
        fwhm = (fwhm + ext_ts[last_i]-l_max_t)/2
    else:
        for i, v in enumerate(ext_vs):
            if fwhm != 0:
                fwhm = np.abs(fwhm)
                break
            #ensure that the fwhm is set before the maxima
            if v > l_max_v/2:
                fwhm = l_max_t-ext_ts[i]
            elif v == l_max_v and i > 0:
                fwhm = l_max_t-ext_ts[i-1]
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
    return np.abs(x[0])*np.exp(-0.5*((t_pts - x[1])/x[2])**2)

#return a Gaussian pulse series with parameters specified by x (amplitude, mean, 2*std_dev**2, phase)
def gauss_series(x, t_pts):
    time_difs = t_pts - x[1]
    return np.abs(x[0])*np.exp(-0.5*(time_difs/x[2])**2)*np.cos(x[3]*time_difs - x[4])

#return a Gaussian pulse series with parameters specified by x (amplitude, mean, 2*std_dev**2, phase)
def double_gauss_series(x, t_pts):
    time_difs_0 = t_pts - x[1]
    time_difs_1 = t_pts - x[6]
    return (np.abs(x[0])*np.exp(-0.5*(time_difs_0/x[2])**2) + np.abs(x[5])*np.exp(-0.5*(time_difs_1/x[7])**2))*np.cos(x[3]*time_difs_0 - x[4])

def find_local_maxima(t_pts, a_pts, cutoff=AMP_CUTOFF):
    '''Given t_pts and a_pts return a list off all local extrema and their corresponding times. Also return the time for the global maxima, the value of the global maxima and the average spacing between maxima
    t_pts: time points to search
    a_pts: values to search
    cutoff: if this is set to a value greater than zero, then all points with an amplitude less than cutoff*global_max_a will be removed
    '''
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
    #remove everything that is below a user specified cutoff threshold
    if cutoff > 0:
        cutoff *= l_max_v
        i_min = -1
        for i, v in enumerate(ext_vs):
            if i_min < 0:
                if v > cutoff:
                    i_min = i
            elif v < cutoff and i >= l_max_i:
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

    return ext_ts, ext_vs, l_max_t, l_max_v, typ_spacing

def opt_envelope(t_pts, a_pts, keep_n=-1, fig_name=''):
    '''Performs a least squares fit to find a rough estimate of the shape of the pulse envelope. This is used by the opt_pulse to optimize the pulse shape inside the envelope
    t_pts: a numpy array of times at which each point occurred
    a_pts: the values corresponding to the t_pts
    keep_n: if set to a number greater than 0, then only the first keep_n standard deviations are used for inferring the envelope
    fig_name: the filename to save the best fit curve to
    returns: a tuple containing (np.array([times of maxima, maxima]), global time maximum, global maximum, the spacing between adjacent peaks)
    '''
    #find the extrema
    ext_ts, ext_vs, l_max_t, l_max_v, typ_spacing = find_local_maxima(t_pts, a_pts)
    #use the full width at half max as a crude estimate of standard deviation
    fwhm, ext_ts, ext_vs = cut_devs(ext_ts, ext_vs, l_max_t, l_max_v, keep_n=keep_n)

    #now take the logarithm to transform fitting a gaussian to fitting a parabola. We normalize so that the amplitude should be roughly 1.0. This is easily accounted for at the end
    ext_ts = np.array(ext_ts)
    ext_vs = np.log(np.array(ext_vs)/l_max_v)
    #least squares of the logs of the envelope points
    def ff(x):
        return np.sum( (np.log(abs(x[0])) - 0.5*((ext_ts - x[1])/x[2])**2 - ext_vs)**2 )
    #gradient of the least squares of the envelope points
    def jac_env(x):
        diff_ser = (np.log(abs(x[0])) - 0.5*((ext_ts - x[1])/x[2])**2 - ext_vs)
        t_dif_by_w = (ext_ts - x[1])
        ret = np.zeros(3)
        ret[0] = 2*np.sum( diff_ser/abs(x[0]) )
        ret[1] = 2*np.sum( diff_ser*t_dif_by_w )/(x[2]**2)
        ret[2] = 2*np.sum( diff_ser*(t_dif_by_w**2) )/(x[2]**3)
        return ret

    #note that we fit to a form e**(-(t-t_0)^2/w) so using fwhm as an estimate for sigma we have w=2*sigma**2
    x0 = np.array([1.0, l_max_t, fwhm])
    res = opt.minimize(ff, x0, jac=jac_env)

    if fig_name != '':
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.scatter(ext_ts, ext_vs, color='black')
        int_ts = np.linspace(ext_ts[0], ext_ts[-1])
        ax.plot(int_ts, np.log(abs(res.x[0])) - 0.5*((int_ts - res.x[1])/res.x[2])**2)
        fig.savefig(fig_name+".png", dpi=300)
        plt.close(fig)

    #we normalized all of our points by the global maxima, so we need to correct that
    res.x[0] = np.abs(res.x[0])*l_max_v
    #now create crude estimates of omega and phi by looking at the spacing between peaks and the phase offset between the data peak and the estimate for t_0
    est_omega = np.pi/typ_spacing
    est_phi = fix_angle(est_omega*(l_max_t - res.x[1]))

    return res, est_omega, est_phi

def opt_double_envelope(t_pts, a_pts, fig_name=''):
    '''Optimize an envelope that is a superposition of two gaussians'''
    #find the extrema
    ext_ts, ext_vs, l_max_t, l_max_v, typ_spacing = find_local_maxima(t_pts, a_pts)
    fwhm,_,_ = cut_devs(ext_ts, ext_vs, l_max_t, l_max_v, both_tails=True)

    #normalize the amplitude
    ext_ts = np.array(ext_ts)
    ext_vs = np.array(ext_vs)/l_max_v

    local_maxes = []
    for i in range(1, len(ext_vs)-1):
        if ext_vs[i] > 0.3 and ext_vs[i] > ext_vs[i-1] and ext_vs[i] > ext_vs[i+1]:
            local_maxes.append(i)

    #the square error
    def doub_gauss_ser(x, t_pts):
        return x[0]*np.exp(-0.5*((t_pts - x[1])/x[2])**2) + x[3]*np.exp(-0.5*((t_pts - x[4])/x[5])**2)
    def ff(x):
        return np.sum( (x[0]*np.exp(-0.5*((ext_ts - x[1])/x[2])**2) + x[3]*np.exp(-0.5*((ext_ts - x[4])/x[5])**2) - ext_vs)**2 )
    def jac_env(x):
        t_dif_0 = (ext_ts - x[1])
        t_dif_1 = (ext_ts - x[4])
        exp_0 = np.exp(-0.5*(t_dif_0/x[2])**2)
        exp_1 = np.exp(-0.5*(t_dif_0/x[5])**2)
        diff_ser = x[0]*exp_0 + x[3]*exp_1 - ext_vs
        ret = np.zeros(6)
        ret[0] = 2*np.sum(diff_ser*exp_0)
        ret[1] = 2*x[0]*np.sum(diff_ser*t_dif_0*exp_0)/(x[2]**2)
        ret[2] = 2*x[0]*np.sum(diff_ser*exp_0*(t_dif_0)**2)/(x[2]**3)
        ret[3] = 2*np.sum(diff_ser*exp_1)
        ret[4] = 2*x[3]*np.sum(diff_ser*t_dif_1*exp_1)/(x[5]**2)
        ret[5] = 2*x[3]*np.sum(diff_ser*exp_1*(t_dif_1)**2)/(x[5]**3)
        return ret

    #create an initial guess by supposing that the widths are each one half
    x0 = np.array([1.0, l_max_t, fwhm/2, 1.0, l_max_t+fwhm, fwhm/2])
    #try a different guess if we have multiple local maxima
    if len(local_maxes) > 1:
        est_fwhm = 0.5*(ext_ts[local_maxes[-1]] - ext_ts[local_maxes[0]] - fwhm)
        if est_fwhm > 0:
            x0 = np.array([ext_vs[local_maxes[0]], ext_ts[local_maxes[0]], est_fwhm, ext_vs[local_maxes[-1]], ext_ts[local_maxes[-1]], est_fwhm])
    #perform the optimization and make figures
    res = opt.minimize(ff, x0)
    if fig_name != '':
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.scatter(ext_ts, ext_vs, color='black')
        t_range = ext_ts[-1] - ext_ts[0]
        t_cent = (ext_ts[-1] + ext_ts[0])/2
        int_ts = np.linspace(t_cent-t_range, t_cent+t_range, num=200)
        ax.plot(int_ts, doub_gauss_ser(res.x, int_ts))
        ax.plot(int_ts, doub_gauss_ser(x0, int_ts), color='red')
        fig.savefig(fig_name+".png", dpi=300)
        plt.close(fig)

    res.x[0] = np.abs(res.x[0])*l_max_v
    res.x[3] = np.abs(res.x[3])*l_max_v
    #now create crude estimates of omega and phi by looking at the spacing between peaks and the phase offset between the data peak and the estimate for t_0
    est_omega = np.pi/typ_spacing
    est_phi = fix_angle(est_omega*(ext_ts[local_maxes[0]] - res.x[1]))
    return res, est_omega, est_phi

def opt_pulse(t_pts_0, a_pts_0, env_x, est_omega, est_phi, keep_n=-1):
    '''Set values for the a_pts for each time point
    a_pts: the values of the electric field at each corresponding time point in t_pts. Note that len(a_pts) must be the same as len(t_pts)
    env_x: an array representing the value found from an envelope fit (i.e if you call res=opt_envelope(...) this would be res.x
    est_omega: a guess of the value of omega based on the envelope fitting
    est_phi: a guess of the value of the CEP based on the envelope fitting
    returns: the fit result along with the lower and upper times used to perform the fitting'''

    #truncate so that we only look at N_STDS deviations from the center of the pulse
    if keep_n > 0:
        dt = (t_pts_0[-1]-t_pts_0[0])/t_pts_0.shape[0]
        i_0 = int(env_x[1]/dt)
        i_shift = int(keep_n*env_x[2] / dt)
        #make sure we don't go past the ends of the array
        if i_0 < i_shift:
            i_shift = i_0
        if i_0+i_shift > t_pts_0.shape[0]:
            i_shift = t_pts_0.shape[0] - i_0 - 1
        t_pts = t_pts_0[i_0-i_shift:i_0+i_shift]
        a_pts = a_pts_0[i_0-i_shift:i_0+i_shift]
        if len(t_pts) == 0:
            t_pts = t_pts_0
            a_pts = a_pts_0
    else:
        t_pts = t_pts_0
        a_pts = a_pts_0

    #it is helpful to try to get the amplitude close to unity so that the gradient in other directions is not suppressed
    if env_x[0] != 0.0:
        a_pts = a_pts / env_x[0]

    #square error of the pulse envelope
    def fp(x):
        if x[2] == 0:
            x[2] = 0.1
        return np.sum( (gauss_series(x, t_pts) - a_pts)**2 )
    def jac_fp(x):
        if x[2] == 0:
            x[2] = 0.1
        p_i = gauss_series(x, t_pts) - a_pts
        exp_cos = np.exp(-0.5*((t_pts-x[1])/x[2])**2)*np.cos(x[3]*(t_pts-x[1])-x[4])
        exp_sin = np.exp(-0.5*((t_pts-x[1])/x[2])**2)*np.sin(x[3]*(t_pts-x[1])-x[4])       
        t_dif_by_w = (t_pts-x[1])
        ret = np.zeros(5)
        ret[0] = 2*np.sum(p_i*exp_cos)
        ret[1] = 2*x[0]*np.sum(p_i*(t_dif_by_w*exp_cos/(x[2]**2) + exp_sin))
        ret[2] = 2*x[0]*np.sum(p_i*exp_cos*t_dif_by_w**2)/(x[2]**3)
        ret[3] = -2*x[0]*np.sum(p_i*exp_sin*t_dif_by_w)
        ret[4] = 2*x[0]*np.sum(p_i*exp_sin)
        return ret

    x0 = np.array([1.0, env_x[1], env_x[2], est_omega, est_phi])
    #super hacky way to account for pi phase shifts
    fp_old = fp(x0)
    x0[4] += np.pi
    if fp(x0) > fp_old:
        x0[4] = est_phi
    else:
        est_phi += np.pi
    x0[4] = fix_angle(x0[4])

    #now we actually perform the minimization
    res = opt.minimize(fp, x0, jac=jac_fp)

    #renormalize thingies
    res.x[0] = res.x[0]*env_x[0]
    res.x[4] = fix_angle(res.x[4])

    return res, t_pts[0], t_pts[-1]

def opt_double_pulse(t_pts, a_pts, env_x, est_omega, est_phi, keep_n=-1):
    '''Set values for the a_pts for each time point
    a_pts: the values of the electric field at each corresponding time point in t_pts. Note that len(a_pts) must be the same as len(t_pts)
    env_x: an array representing the value found from an envelope fit (i.e if you call res=opt_envelope(...) this would be res.x
    est_omega: a guess of the value of omega based on the envelope fitting
    est_phi: a guess of the value of the CEP based on the envelope fitting
    returns: the fit result along with the lower and upper times used to perform the fitting'''
    #it is helpful to try to get the amplitude close to unity so that the gradient in other directions is not suppressed
    if env_x[0] != 0.0:
        a_pts = a_pts / env_x[0]

    #square error of the pulse envelope
    def fp(x):
        return np.sum( (double_gauss_series(x, t_pts) - a_pts)**2 )
    def jac_fp(x):
        p_i = double_gauss_series(x, t_pts) - a_pts
        exp_cos_0 = x[0]*np.exp(-0.5*((t_pts-x[1])/x[2])**2)*np.cos(x[3]*(t_pts-x[1])-x[4])
        exp_sin_0 = x[0]*np.exp(-0.5*((t_pts-x[1])/x[2])**2)*np.sin(x[3]*(t_pts-x[1])-x[4])
        t_dif_0 = (t_pts-x[1])
        exp_cos_1 = x[5]*np.exp(-0.5*((t_pts-x[6])/x[7])**2)*np.cos(x[3]*(t_pts-x[1])-x[4])
        exp_sin_1 = x[5]*np.exp(-0.5*((t_pts-x[6])/x[7])**2)*np.sin(x[3]*(t_pts-x[1])-x[4])
        t_dif_1 = (t_pts-x[6])
        ret = np.zeros(8)
        ret[0] = 2*np.sum(p_i*exp_cos_0)/x[0]
        ret[1] = 2*np.sum(p_i*(t_dif_0*exp_cos_0/(x[2]**2) + x[6]*(exp_sin_0+exp_sin_1)))
        ret[2] = 2*np.sum(p_i*exp_cos_0*t_dif_0**2)/(x[2]**3)
        ret[3] =-2*np.sum(p_i*t_dif_1*(exp_sin_0+exp_sin_1))
        ret[4] = 2*np.sum(p_i*(exp_sin_0+exp_sin_1))
        ret[5] = 2*np.sum(p_i*exp_cos_1)/x[5]
        ret[6] = 2*np.sum(p_i*(t_dif_1*exp_cos_1/(x[7]**2)))
        ret[7] = 2*np.sum(p_i*exp_cos_1*t_dif_1**2)/(x[7]**3)
        return ret

    x0 = np.array([1.0, env_x[1], env_x[2], est_omega, est_phi, env_x[3]/env_x[0], env_x[4], env_x[5]])
    #super hacky way to account for pi phase shifts
    fp_old = fp(x0)
    x0[4] += np.pi
    if fp(x0) > fp_old:
        x0[4] = est_phi
    else:
        est_phi += np.pi
    x0[4] = fix_angle(x0[4])

    #now we actually perform the minimization
    res = opt.minimize(fp, x0)

    #renormalize thingies
    res.x[0] = np.abs(res.x[0])*env_x[0]
    res.x[5] = np.abs(res.x[5])*env_x[0]
    #res.x[4] = fix_angle(res.x[4])

    return res, t_pts[0], t_pts[-1]

def opt_pulse_env(t_pts_0, a_pts_0, a_sigmas=0.1, keep_n=-1, fig_name=''):
    if a_pts_0.shape != t_pts_0.shape:
        raise ValueError("t_pts and a_pts must have the same shape")
    env_fig_name = ''
    if fig_name:
        env_fig_name = fig_name+"_env"

    t_pts_skip = t_pts_0[::SKIP]
    a_pts_skip = a_pts_0[::SKIP]
    #set defaults
    res = None
    low_t = t_pts_0[0]
    hi_t = t_pts_0[-1]
    used_double = False
    #fit the pulse envelope, this gives us a good estimate of three points
    env_res, est_omega, est_phi = opt_envelope(t_pts_skip, a_pts_skip, keep_n=keep_n, fig_name=env_fig_name)
    #try fitting two gaussian pulses if that was bad
    if env_res.fun > 0.5:
        new_env_res, est_omega, est_phi = opt_double_envelope(t_pts_skip, a_pts_skip, fig_name=env_fig_name+"_2")
        #compute the bayes factor for the single and double pulses. We take both models to have identical priors. Note that under assumption of Gaussian errors the error cancels
        orig_er_sq = np.sum( (env_res.x[0]*np.exp(-0.5*((t_pts_skip - env_res.x[1])/env_res.x[2])**2) - a_pts_skip)**2 )
        bayes = np.exp(-0.5*new_env_res.fun/(a_sigmas**2))/np.exp(-0.5*orig_er_sq/(a_sigmas**2))
        print("Found bad envelope fit, testing bayes factor=%f" % bayes)
        #substantial evidence as defined by Jeffreys
        if bayes > 3.16:
            env_res = new_env_res
            res, low_t, hi_t = opt_double_pulse(t_pts_0, a_pts_0, env_res.x, est_omega, est_phi, keep_n=keep_n)
            used_double = True
        '''#now we try to suppress the second envelope by dividing by fitting with the first envelope
        res, low_t, hi_t = opt_pulse(t_pts_0, a_pts_0, np.resize(env_res.x, 3), est_omega, est_phi, keep_n=keep_n)'''

    if res is None:
        res, low_t, hi_t = opt_pulse(t_pts_0, a_pts_0, env_res.x, est_omega, est_phi, keep_n=keep_n)

    #plot figures if requested
    if fig_name != '':
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.plot(t_pts_0, a_pts_0, color='black')
        ax.plot(t_pts_0, gauss_env(env_res.x, t_pts_0), linestyle=':', color='gray')
        if used_double:
            x0 = np.array([env_res.x[0], env_res.x[1], env_res.x[2], est_omega, est_phi, env_res.x[3], env_res.x[4], env_res.x[5]])
            ax.plot(t_pts_0, double_gauss_series(x0, t_pts_0), color='blue')
            ax.plot(t_pts_0, double_gauss_series(res.x, t_pts_0), color='red')
        else:
            ax.plot(t_pts_0, gauss_series(res.x, t_pts_0), color='red')
        ax.vlines([low_t, hi_t], -1.2, 1.2)
        ax.set_ylim((-1.2, 1.2))
        ax.annotate(r"$er^2={:.1f}$".format(res.fun), (0.01, 0.84), xycoords='axes fraction')
        fig.savefig(fig_name+".png", dpi=300)
        plt.close(fig)

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
