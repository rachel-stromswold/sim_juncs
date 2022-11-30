import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.fft import irfft, rfft, rfftfreq

import matplotlib
import matplotlib.pyplot as plt

EPSILON = 0.125
H_EPSILON = EPSILON/2

SKIP = 200
AMP_CUTOFF = 0.0025
N_STDS = 2
WIDTH_SCALE = 1000

DOUB_SWAP_MAT = np.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]])

#for a local maxima with amplitude A_l/A_g >= LOC_MAX_CUT, the location of A_l will be considered for double envelope optimization. Here A_l is the amplitude of the local maxima and A_g is the amplitude of the global maxima
LOC_MAX_CUT = 0.2
#The number of standard deviations away from the local maxima to consider when performing least squares fits
DEF_KEEP_N = 2.5

'''def jac_num(func, x):
    #Numerically estimate the jacobian of func at point f
    dim = len(x)
    ret = np.zeros(dim)
    for i in range(dim):
        x[i] += H_EPSILON
        ret[i] = func(x)
        x[i] -= EPSILON
        ret[i] = (ret[i] - func(x))/EPSILON
        x[i] += H_EPSILON
    return ret'''

def write_reses(fname, res_list):
    with open(fname, "w") as res_f:
        for res in res_list:
            for reskey in res.keys():
                res_f.write(reskey)
                val = res.get(reskey)
                #formating matrices is a special case
                if isinstance(val, np.ndarray) and len(val.shape) == 2:
                    res_f.write("\n")
                    for row in val:
                        res_f.write("\t[")
                        for v in row:
                            if v >= 0:
                                res_f.write(" ")
                            res_f.write(" {:.3E} ".format(v))
                        res_f.write("]\n")
                else:
                    res_f.write(": {}\n".format(val))
            res_f.write("\n==========\n")

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

class EnvelopeFitter:
    def __init__(self, t_pts, a_pts, cutoff=AMP_CUTOFF):
        '''Given t_pts and a_pts return a list off all local extrema and their corresponding times. Also return the time for the global maxima, the value of the global maxima and the average spacing between maxima
        t_pts: time points to search
        a_pts: values to search
        cutoff: if this is set to a value greater than zero, then all points with an amplitude less than cutoff*global_max_a will be removed
        '''
        self.ext_ts = []
        self.ext_vs = []
        self.max_t = 0
        self.max_v = 0
        self.l_max_i = 0
        #typ_spacing is used later for estimating omega
        self.typ_spacing = 0.0
        for i in range(1, len(a_pts)-1):
            #extrema
            if (a_pts[i] > a_pts[i-1] and a_pts[i] > a_pts[i+1]) or (a_pts[i] < a_pts[i-1] and a_pts[i] < a_pts[i+1]):
                self.ext_ts.append(t_pts[i])
                self.ext_vs.append(abs(a_pts[i]))
                #check if this is a global maxima
                if self.ext_vs[-1] > self.max_v:
                    self.l_max_i = i
                    self.max_t = t_pts[i]
                    self.max_v = self.ext_vs[-1]
                    #set the typical spacing by the distance to the nearest earlier peak
                    if len(self.ext_ts) > 1:
                        self.typ_spacing = self.ext_ts[-1] - self.ext_ts[-2]
        #remove everything that is below a user specified cutoff threshold
        if cutoff > 0:
            cutoff *= self.max_v
            i_min = -1
            for i, v in enumerate(self.ext_vs):
                if i_min < 0:
                    if v > cutoff:
                        i_min = i
                elif v < cutoff and i >= self.l_max_i:
                    #make sure we include at least this data point
                    self.ext_ts = self.ext_ts[i_min:i]
                    self.ext_vs = self.ext_vs[i_min:i]
                    break
        #make extra sure there aren't any divisions by zero
        if self.typ_spacing == 0.0:
            #this means there isn't enough data to perform an estimate of the spacing, just use omega_0
            if len(self.ext_ts) < 2:
                self.typ_spacing = np.pi/omega_0
            else:
                self.typ_spacing = self.ext_ts[1] - self.ext_ts[0]

    def cut_devs(self, keep_n=DEF_KEEP_N, both_tails=False):
        '''estimate the standard deviation of a distribution by taking the full width at half maximum and then return this value along with copies of t_pts and a_pts truncated to only include points within keep_n deviations from the center.
        ext_ts: x axis values of the distribution
        ext_vs: pdf
        keep_n: the number of standard deviations (estimated by fwhms) to keep, use a negative value to not perform any truncation
        returns: a tuple of three values (fwhm, ext_tsated, ext_vsated)
        '''
        fwhm = 0
        if both_tails:
            last_i = 0
            for i, v in enumerate(self.ext_vs):
                #ensure that the fwhm is set before the maxima
                if v == self.max_v and fwhm == 0 and i > 0:
                    fwhm = self.max_t-self.ext_ts[i-1]
                if fwhm == 0 and v > self.max_v/2:
                    #use this point if there weren't any previously. Otherwise lerp
                    if i == 0:
                        fwhm = self.max_t-self.ext_ts[i]
                    else:
                        t_mid = self.ext_ts[i-1] + (0.5*self.max_v-self.ext_vs[i-1])*(self.ext_ts[i]-self.ext_ts[i-1])/(self.ext_vs[i]-self.ext_vs[i-1])
                        fwhm = self.max_t-t_mid
                elif fwhm != 0 and v > self.max_v/2:
                    last_i = i
            #use the last point if there weren't any after. Otherwise lerp
            if last_i == len(self.ext_vs)-1:
                fwhm = (fwhm + self.ext_ts[last_i] - self.max_t)/2
            else:
                t_mid = self.ext_ts[last_i] + (self.ext_vs[last_i]-0.5*self.max_v)*(self.ext_ts[last_i+1]-self.ext_ts[last_i])/(self.ext_vs[last_i]-self.ext_vs[last_i+1])
                fwhm = (fwhm + t_mid - self.max_t)/2
        else:
            for i, v in enumerate(self.ext_vs):
                if fwhm != 0:
                    fwhm = np.abs(fwhm)
                    break
                #ensure that the fwhm is set before the maxima
                if v > self.max_v/2:
                    fwhm = self.max_t-self.ext_ts[i]
                elif v == self.max_v and i > 0:
                    fwhm = self.max_t-self.ext_ts[i-1]
        #check if a truncation must be performed
        if keep_n < 0:
            return fwhm, self.ext_ts, self.ext_vs
        #do it
        i_min = -1
        for i, t in enumerate(self.ext_ts):
            if i_min < 0:
                if t > self.max_t - fwhm*keep_n:
                    i_min = i
            elif t > self.max_t + fwhm*keep_n:
                return fwhm, self.ext_ts[i_min:i], self.ext_vs[i_min:i]
        #make sure that something is returned
        if i_min < 0:
            return fwhm, self.ext_ts, self.ext_vs
        else:
            return fwhm, self.ext_ts[i_min:], self.ext_vs[i_min:]

    def opt_envelope(self, keep_n=DEF_KEEP_N, fig_name=''):
        '''Performs a least squares fit to find a rough estimate of the shape of the pulse envelope. This is used by the opt_pulse to optimize the pulse shape inside the envelope
        t_pts: a numpy array of times at which each point occurblue
        a_pts: the values corresponding to the t_pts
        keep_n: if set to a number greater than 0, then only the first keep_n standard deviations are used for inferring the envelope
        fig_name: the filename to save the best fit curve to
        returns: a tuple containing (np.array([times of maxima, maxima]), global time maximum, global maximum, the spacing between adjacent peaks)
        '''
        #use the full width at half max as a crude estimate of standard deviation
        fwhm, trunc_ts, trunc_vs = self.cut_devs(keep_n=keep_n)

        #now take the logarithm to transform fitting a gaussian to fitting a parabola. We normalize so that the amplitude should be roughly 1.0. This is easily accounted for at the end
        trunc_ts = np.array(trunc_ts)
        trunc_vs = np.log(np.array(trunc_vs)/self.max_v)
        #least squares of the logs of the envelope points
        def ff(x):
            return np.sum( (np.log(abs(x[0])) - 0.5*((trunc_ts - x[1])/x[2])**2 - trunc_vs)**2 )
        #gradient of the least squares of the envelope points
        def jac_env(x):
            diff_ser = (np.log(abs(x[0])) - 0.5*((trunc_ts - x[1])/x[2])**2 - trunc_vs)
            t_dif_by_w = (trunc_ts - x[1])
            ret = np.zeros(3)
            ret[0] = 2*np.sum( diff_ser/abs(x[0]) )
            ret[1] = 2*np.sum( diff_ser*t_dif_by_w )/(x[2]**2)
            ret[2] = 2*np.sum( diff_ser*(t_dif_by_w**2) )/(x[2]**3)
            return ret

        #note that we fit to a form e**(-(t-t_0)^2/w) so using fwhm as an estimate for sigma we have w=2*sigma**2
        x0 = np.array([1.0, self.max_t, fwhm])
        res = opt.minimize(ff, x0, jac=jac_env)

        if fig_name != '':
            fig = plt.figure()
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            ax.scatter(trunc_ts, trunc_vs, color='black')
            int_ts = np.linspace(trunc_ts[0], trunc_ts[-1])
            ax.plot(int_ts, np.log(abs(res.x[0])) - 0.5*((int_ts - res.x[1])/res.x[2])**2)
            fig.savefig(fig_name+".png", dpi=300)
            plt.close(fig)

        #we normalized all of our points by the global maxima, so we need to correct that
        res.x[0] = np.abs(res.x[0])*self.max_v
        #now create crude estimates of omega and phi by looking at the spacing between peaks and the phase offset between the data peak and the estimate for t_0
        est_omega = np.pi/self.typ_spacing
        est_phi = fix_angle(est_omega*(self.max_t - res.x[1]))

        return res, est_omega, est_phi

    def opt_double_envelope(self, fig_name=''):
        '''Optimize an envelope that is a superposition of two gaussians'''
        fwhm,_,_ = self.cut_devs(both_tails=True)
        print("fwhm = %f" % fwhm)

        #normalize the amplitude
        l_ext_ts = np.array(self.ext_ts)
        l_ext_vs = np.array(self.ext_vs)/self.max_v

        local_maxes = []
        for i in range(1, len(l_ext_vs)-1):
            if l_ext_vs[i] > LOC_MAX_CUT and l_ext_vs[i] > l_ext_vs[i-1] and l_ext_vs[i] > l_ext_vs[i+1]:
                local_maxes.append(i)

        #the square error
        def doub_gauss_ser(x, t_pts):
            return x[0]*np.exp(-0.5*((t_pts - x[1])/x[2])**2) + x[3]*np.exp(-0.5*((t_pts - x[4])/x[5])**2)
        def ff(x):
            return np.sum( (np.abs(x[0])*np.exp(-0.5*((l_ext_ts - x[1])/x[2])**2) + np.abs(x[3])*np.exp(-0.5*((l_ext_ts - x[4])/x[5])**2) - l_ext_vs)**2 )
        def jac_env(x):
            t_dif_0 = (l_ext_ts - x[1])
            t_dif_1 = (l_ext_ts - x[4])
            exp_0 = np.abs(x[0])*np.exp(-0.5*(t_dif_0/x[2])**2)
            exp_1 = np.abs(x[3])*np.exp(-0.5*(t_dif_1/x[5])**2)
            diff_ser = np.abs(x[0])*exp_0 + np.abs(x[3])*exp_1 - l_ext_vs
            ret = np.zeros(6)
            ret[0] = 2*np.sum(diff_ser*exp_0)/np.abs(x[0])
            ret[1] = 2*np.sum(diff_ser*t_dif_0*exp_0)/(x[2]**2)
            ret[2] = 2*np.sum(diff_ser*exp_0*(t_dif_0)**2)/(x[2]**3)
            ret[3] = 2*np.sum(diff_ser*exp_1)/np.abs(x[3])
            ret[4] = 2*np.sum(diff_ser*t_dif_1*exp_1)/(x[5]**2)
            ret[5] = 2*np.sum(diff_ser*exp_1*(t_dif_1)**2)/(x[5]**3)
            return ret

        #create an initial guess by supposing that the widths are each one half
        x0 = np.array([1.0, self.max_t, fwhm/2, 1.0, self.max_t+fwhm, fwhm/2])
        #try a different guess if we have multiple local maxima
        if len(local_maxes) > 1:
            est_fwhm = 0.5*(l_ext_ts[local_maxes[-1]] - l_ext_ts[local_maxes[0]] - fwhm)
            if est_fwhm > 0:
                x0 = np.array([l_ext_vs[local_maxes[0]], l_ext_ts[local_maxes[0]], est_fwhm, l_ext_vs[local_maxes[-1]], l_ext_ts[local_maxes[-1]], est_fwhm])
        #perform the optimization and make figures
        res = opt.minimize(ff, x0, jac=jac_env)
        if fig_name != '':
            fig = plt.figure()
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            ax.scatter(l_ext_ts, l_ext_vs, color='black')
            t_range = l_ext_ts[-1] - l_ext_ts[0]
            t_cent = (l_ext_ts[-1] + l_ext_ts[0])/2
            int_ts = np.linspace(t_cent-t_range, t_cent+t_range, num=200)
            ax.plot(int_ts, doub_gauss_ser(x0, int_ts), color='red')
            ax.plot(int_ts, doub_gauss_ser(res.x, int_ts), color='blue')
            fig.savefig(fig_name+".png", dpi=300)
            plt.close(fig)

        res.x[0] = np.abs(res.x[0])*self.max_v
        res.x[3] = np.abs(res.x[3])*self.max_v
        #now create crude estimates of omega and phi by looking at the spacing between peaks and the phase offset between the data peak and the estimate for t_0
        est_omega = np.pi/self.typ_spacing
        est_phi = fix_angle(est_omega*(l_ext_ts[local_maxes[0]] - res.x[1]))
        return res, est_omega, est_phi

    def sq_err_fit_single(self, x):
        return np.sum( (abs(x[0])*np.exp(-0.5*((self.ext_ts - x[1])/x[2])**2) - self.ext_vs)**2 )

    def sq_err_fit_double(self, x):
        return np.sum( (np.abs(x[0])*np.exp(-0.5*((self.ext_ts - x[1])/x[2])**2) + np.abs(x[3])*np.exp(-0.5*((self.ext_ts - x[4])/x[5])**2) - self.ext_vs)**2 )

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
    max_t = 0
    max_v = 0
    l_max_i = 0
    #typ_spacing is used later for estimating omega
    typ_spacing = 0.0
    for i in range(1, len(a_pts)-1):
        #extrema
        if (a_pts[i] > a_pts[i-1] and a_pts[i] > a_pts[i+1]) or (a_pts[i] < a_pts[i-1] and a_pts[i] < a_pts[i+1]):
            ext_ts.append(t_pts[i])
            ext_vs.append(abs(a_pts[i]))
            #check if this is a global maxima
            if ext_vs[-1] > max_v:
                l_max_i = i
                max_t = t_pts[i]
                max_v = ext_vs[-1]
                #set the typical spacing by the distance to the nearest earlier peak
                if len(ext_ts) > 1:
                    typ_spacing = ext_ts[-1] - ext_ts[-2]
    #remove everything that is below a user specified cutoff threshold
    if cutoff > 0:
        cutoff *= max_v
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

    return ext_ts, ext_vs, max_t, max_v, typ_spacing

def opt_pulse(t_pts_0, a_pts_0, env_x, est_omega, est_phi, keep_n=DEF_KEEP_N):
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

def fix_double_pulse_order(res):
    '''Adjust the double pulse fit so that the first envelope comes first in time'''
    #ensure that the pulses are ordeblue the way we expect
    if res.x.shape[0] == 8 and res.x[1] > res.x[6]:
        res.x = np.dot(DOUB_SWAP_MAT, res.x)
        res.jac = np.dot(DOUB_SWAP_MAT, res.jac)
        res.hess_inv = np.dot(DOUB_SWAP_MAT, res.hess_inv)
        #the phase adjustment is slightly more complicated
        res.x[4] += res.x[3]*(res.x[1]-res.x[6])
    return res

def opt_double_pulse(t_pts, a_pts, env_x, est_omega, est_phi, keep_n=DEF_KEEP_N):
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
        t_dif_0 = (t_pts-x[1])
        t_dif_1 = (t_pts-x[6])
        exp_cos_0 = np.abs(x[0])*np.exp(-0.5*(t_dif_0/x[2])**2)*np.cos(x[3]*t_dif_0-x[4])
        exp_sin_0 = np.abs(x[0])*np.exp(-0.5*(t_dif_0/x[2])**2)*np.sin(x[3]*t_dif_0-x[4])
        exp_cos_1 = np.abs(x[5])*np.exp(-0.5*(t_dif_1/x[7])**2)*np.cos(x[3]*t_dif_0-x[4])
        exp_sin_1 = np.abs(x[5])*np.exp(-0.5*(t_dif_1/x[7])**2)*np.sin(x[3]*t_dif_0-x[4])
        ret = np.zeros(8)
        ret[0] = 2*np.sum(p_i*exp_cos_0)/x[0]
        ret[1] = 2*np.sum(p_i*(t_dif_0*exp_cos_0/(x[2]**2) + x[3]*(exp_sin_0+exp_sin_1)))
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
    res = opt.minimize(fp, x0, jac=jac_fp)

    #renormalize thingies
    res.x[0] = np.abs(res.x[0])*env_x[0]
    res.x[5] = np.abs(res.x[5])*env_x[0]
    #res.x[4] = fix_angle(res.x[4])

    return fix_double_pulse_order(res), t_pts[0], t_pts[-1]

def opt_pulse_env(t_pts_0, a_pts_0, a_sigmas_sq=0.1, keep_n=DEF_KEEP_N, fig_name=''):
    if a_pts_0.shape != t_pts_0.shape:
        raise ValueError("t_pts and a_pts must have the same shape")
    env_fig_name = ''
    env_fig_name_2 = ''
    if fig_name != '':
        env_fig_name = fig_name+"_env"
        env_fig_name_2 = fig_name+"_env_2"

    t_pts_skip = t_pts_0[::SKIP]
    a_pts_skip = a_pts_0[::SKIP]

    #the maximum of the derivative on a_pts(t_pts) and use this to get a conservative estimate of the variance on each of the sample points. Note that for an RV Y with P(Y|X)=N(mX+x0, s) we have V(Y) = s^2 + m^2 V(X). We take X to follow a continuous uniform distribution between (t_pts[i], t_pts[i+SKIP])
    max_diff = np.max(np.diff(a_pts_skip))
    skip_sigma_sq = a_sigmas_sq + (max_diff**2)/12
    #find the extrema
    env_fit = EnvelopeFitter(t_pts_skip, a_pts_skip)
    #set defaults
    env_res = None
    res = None
    low_t = t_pts_0[0]
    hi_t = t_pts_0[-1]
    err = 1.0
    used_double = False
    #try fitting to both possible envelopes
    env_res_0, est_omega_0, est_phi_0 = env_fit.opt_envelope(keep_n=keep_n, fig_name=env_fig_name)
    env_res_1, est_omega_1, est_phi_1 = env_fit.opt_double_envelope(fig_name=env_fig_name_2) 
    err_0 = env_fit.sq_err_fit_single(env_res_0.x)
    err_1 = env_fit.sq_err_fit_double(env_res_1.x)
    #compute the bayes factor for the single and double pulses. We take both models to have identical priors. Note that under assumption of Gaussian errors the error cancels
    ln_bayes =  (0.5/skip_sigma_sq)*(err_0 - err_1)
    print("\t\tenvelope sigma={}, bayes factor={}".format(skip_sigma_sq, np.exp(ln_bayes)))
    #substantial evidence as defined by Jeffreys
    if np.isnan(ln_bayes) or ln_bayes > 1.15:
        env_res = env_res_1
        res, low_t, hi_t = opt_double_pulse(t_pts_skip, a_pts_skip, env_res_1.x, est_omega_1, est_phi_1, keep_n=keep_n)
        err = err_1
        used_double = True
    else:
        env_res = env_res_0
        res, low_t, hi_t = opt_pulse(t_pts_skip, a_pts_skip, env_res_0.x, est_omega_0, est_phi_0, keep_n=keep_n)
        err = err_0
        used_double = False

    #plot figures if requested
    if fig_name != '':
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.plot(t_pts_0, a_pts_0, color='black')
        ax.plot(t_pts_0, gauss_env(env_res.x, t_pts_0), linestyle=':', color='red')
        ax.plot(t_pts_0, gauss_env(res.x, t_pts_0), linestyle=':', color='blue')
        if used_double:
            x0 = np.array([env_res.x[0], env_res.x[1], env_res.x[2], est_omega_1, est_phi_1, env_res.x[3], env_res.x[4], env_res.x[5]])
            ax.plot(t_pts_0, double_gauss_series(x0, t_pts_0), color='red')
            ax.plot(t_pts_0, double_gauss_series(res.x, t_pts_0), color='blue')
        else:
            x0 = np.array([env_res.x[0], env_res.x[1], env_res.x[2], est_omega_0, est_phi_0])
            ax.plot(t_pts_0, gauss_series(x0, t_pts_0), color='red')
            ax.plot(t_pts_0, gauss_series(res.x, t_pts_0), color='blue')
        ax.vlines([low_t, hi_t], -1.2, 1.2)
        ax.set_ylim((-1.2, 1.2))
        ax.annotate(r"$er^2={:.1f}$".format(res.fun), (0.01, 0.84), xycoords='axes fraction')
        fig.savefig(fig_name+".png", dpi=300)
        plt.close(fig)

    if fig_name != '':
        res_name = fig_name+"_res.txt"
        write_reses(res_name, [res,env_res])

    return res, env_res

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
