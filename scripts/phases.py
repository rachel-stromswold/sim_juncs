import numpy as np
import argparse
import h5py
import scipy.optimize as opt
#from scipy.fft import ifft, fft, fftfreq
import scipy.fft as fft
import scipy.signal as ssig
import time
import pickle
import os.path

import matplotlib
import matplotlib.pyplot as plt

EPSILON = 0.125
H_EPSILON = EPSILON/2
AMP_CUTOFF = 0.025
PEAK_AMP_THRESH = 0.1
NOISE_THRESH = 0.1
OUTLIER_FACTOR = 20
N_STDS = 2
WIDTH_SCALE = 1000
verbose = 4

#The highest waveguide mode (n) that is considered in fits sin((2n-1)pi/L)
HIGHEST_MODE=3

DOUB_SWAP_MAT = np.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]])

#for a local maxima with amplitude A_l/A_g >= LOC_MAX_CUT, the location of A_l will be considered for double envelope optimization. Here A_l is the amplitude of the local maxima and A_g is the amplitude of the global maxima
LOC_MAX_CUT = 1e-6
#The number of standard deviations away from the local maxima to consider when performing least squares fits
DEF_KEEP_N = 2.5

MAX_LOWPASS_EVALS = 4

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

#checks whether the point at index j is a local maxima of the function narr
def is_max(narr, j):
    if j >= narr.shape[0]:
        return False
    if (j == 0):
        return (narr[j] > narr[j+1])
    if (j == narr.shape[0]-1):
        return (narr[j] > narr[j-1])
    return (narr[j] > narr[j-1] and narr[j] > narr[j+1])

def deriv_num_1(t_pts, f_pts, i):
    '''Numerically estimate the first derivative of f_pts (with corresponding points in time t_pts) evaluated at the ith time index'''
    if len(f_pts) == 0:
        return 0
    elif i == 0:
        return (f_pts[i+1]-f_pts[i])/(t_pts[i+1]-t_pts[i])
    elif i == len(f_pts)-1:
        return (f_pts[i]-f_pts[i-1])/(t_pts[i]-t_pts[i-1])
    #we don't know that times are equally spaced so we consider both adjacent points seperately
    return 0.5*( (f_pts[i+1] - f_pts[i])/(t_pts[i+1]-t_pts[i]) + (f_pts[i] - f_pts[i-1])/(t_pts[i]-t_pts[i-1]) )

def deriv_num_2(t_pts, f_pts, i):
    '''Numerically estimate the second derivative of f_pts (with corresponding points in time t_pts) evaluated at the ith time index'''
    if len(f_pts) == 0:
        return 0
    elif i == 0:
        return 2*(f_pts[i+1]-f_pts[i])/((t_pts[i+1]-t_pts[i])**2)
    elif i == len(f_pts)-1:
        return 2*(f_pts[i-1]-f_pts[i])/((t_pts[i]-t_pts[i-1])**2)
    pside = (f_pts[i+1] - f_pts[i])/(t_pts[i+1]-t_pts[i])
    mside = (f_pts[i] - f_pts[i-1])/(t_pts[i]-t_pts[i-1])
    tdif = t_pts[i+1]-t_pts[i-1]
    ret =(pside-mside)/tdif
    return ret

def est_width(t_pts, f_pts, i, default=1):
    '''Treating i as t_pts[i] as the center of a Gaussian distribution, estimate the width based on f_pts'''
    der = deriv_num_2(t_pts, f_pts, i)
    #scan right and then left to find a point with a negative derivative
    m = 1
    while der >= 0:
        #if we move past the end of the array return 1
        if m > i or i+m >= len(t_pts):
            return default
        der = min(deriv_num_2(t_pts,f_pts,i-m), deriv_num_2(t_pts,f_pts,i+m))
        m += 1
    return np.sqrt(-1/der)

class EnvelopeFitter:
    def __init__(self, t_pts, a_pts, cutoff=AMP_CUTOFF, outlier=OUTLIER_FACTOR):
        '''Given t_pts and a_pts return a list off all local extrema and their corresponding times. Also return the time for the global maxima, the value of the global maxima and the average spacing between maxima
        t_pts: time points to search
        a_pts: values to search
        cutoff: if this is set to a value greater than zero, then all points with an amplitude less than cutoff*global_max_a will be removed
        '''
        #make sure we're dealing with real amplitudes instead of complex numbers when finding maxima
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

        #convert to numpy arrays to make life easier
        self.ext_ts = np.array(self.ext_ts)
        self.ext_vs = np.array(self.ext_vs)

    def find_fwhm(self, both_tails=False):
        fwhm = 0
        if both_tails:
            last_i = 0
            for i, v in enumerate(self.ext_vs):
                #ensure that the fwhm is set before the maxima
                if v == self.max_v and fwhm == 0:
                    if i > 0:
                        fwhm = self.max_t-self.ext_ts[i-1]
                    else:
                        fwhm = self.ext_ts[i+1]-self.max_t
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
                elif v == self.max_v:
                    if i > 0:
                        fwhm = self.max_t-self.ext_ts[i-1]
                    else:
                        fwhm = self.ext_ts[i+1]-self.max_t
        return fwhm

    def cut_devs(self, keep_n=DEF_KEEP_N, both_tails=False):
        '''estimate the standard deviation of a distribution by taking the full width at half maximum and then return this value along with copies of t_pts and a_pts truncated to only include points within keep_n deviations from the center.
        ext_ts: x axis values of the distribution
        ext_vs: pdf
        keep_n: the number of standard deviations (estimated by fwhms) to keep, use a negative value to not perform any truncation
        returns: a tuple of three values (fwhm, ext_tsated, ext_vsated)
        '''
        fwhm = self.find_fwhm(both_tails=both_tails)
        if fwhm == 0:
            fwhm = np.sqrt(-1/deriv_num_2(self.ext_ts, self.ext_vs, np.argmax(self.ext_vs)))
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

    def get_single_guess(self, keep_n):
        #use the full width at half max as a crude estimate of standard deviation
        fwhm, trunc_ts, trunc_vs = self.cut_devs(keep_n=keep_n)
        #now take the logarithm to transform fitting a gaussian to fitting a parabola. We normalize so Open Choicethat the amplitude should be roughly 1.0. This is easily accounted for at the end
        trunc_ts = np.array(trunc_ts)
        trunc_vs = np.log(np.array(trunc_vs)/self.max_v)
        if verbose > 2:
            print("\tsingle_guess: peak={}, width={}".format(self.max_t, fwhm))
        #note that we fit to a form e**(-(t-t_0)^2/w) so using fwhm as an estimate for sigma we have w=2*sigma**2
        return np.array([1.0, self.max_t, fwhm]), trunc_ts, trunc_vs

    def opt_envelope(self, keep_n=DEF_KEEP_N, fig_name=''):
        '''Performs a least squares fit to find a rough estimate of the shape of the pulse envelope. This is used by the opt_pulse to optimize the pulse shape inside the envelope
        t_pts: a numpy array of times at which each point occurblue
        a_pts: the values corresponding to the t_pts
        keep_n: if set to a number greater than 0, then only the first keep_n standard deviations are used for inferring the envelope
        fig_name: the filename to save the best fit curve to
        returns: a tuple containing (np.array([times of maxima, maxima]), global time maximum, global maximum, the spacing between adjacent peaks)
        '''
        #use the full width at half max as a crude estimate of standard deviation
        x0, trunc_ts, trunc_vs = self.get_single_guess(keep_n)

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
        if verbose > 2:
            print("\tsingle optimized: {}", res.x)
        return res, est_omega, est_phi

    def find_response(self, kern_x):
        '''Find the response of the convolution of the time series with the kernel specified by kern_x.
        returns: a tuple with the response (f*kern_x)(t) and the local extrema of the resulting convolution
        '''
        n_t_pts = self.ext_ts.shape[0]
        response = np.zeros(n_t_pts)
        #since time samples aren't regularly spaced, convolve in the time domain. n is small right? so O(n^2) isn't thaat bad :p
        for i in range(n_t_pts):
            for j in range(1, n_t_pts):
                response[i] += 0.5*( self.ext_vs[i]*kern_x(self.ext_ts[i] - self.ext_ts[j]) + self.ext_vs[i]*kern_x(self.ext_ts[i] - self.ext_ts[j-1]) ) \
                    *(self.ext_ts[j]-self.ext_ts[j-1])
        peaks = ssig.argrelextrema(response, np.greater)[0]
        return response, peaks


    def get_double_guess(self):
        '''Get a starting guess for the parameters of the double pulse shape
        '''
        fwhm = np.abs( self.find_fwhm(both_tails=True) )
        #ensure that a valid value for the fwhm is used
        if np.isnan(fwhm) or fwhm == 0.0:
            fwhm = 1.0
        if verbose > 1:
            print("\tfwhm = %f" % fwhm)
                
        l_ext_ts = np.array(self.ext_ts)
        l_ext_vs = np.array(self.ext_vs)/self.max_v
        #response, peaks = self.find_response(lambda t: np.exp( -t**2/(2*fwhm**2) ))
        response, peaks = self.find_response(lambda t: np.exp( -t**2/(2*fwhm**2) ) + 2/(1+np.exp(-t/fwhm)))
        #find the two highest peaks in the response function
        local_maxes = [0, 1]
        if response[local_maxes[1]] > response[local_maxes[0]]:
            tmp = local_maxes[0]
            local_maxes[0] = local_maxes[1]
            local_maxes[1] = tmp
        for i in peaks:
            #only proceed if this is a maxima of l_ext_vs
            if response[i] > response[local_maxes[0]]:
                local_maxes[1] = local_maxes[0]
                local_maxes[0] = i
            elif response[i] > response[local_maxes[1]]:
                local_maxes[1] = i
        #now estimate the standard deviation from the second derivative
        local_sigs = [est_width(l_ext_ts, l_ext_vs, local_maxes[0]), est_width(l_ext_ts, l_ext_vs, local_maxes[1])]
        '''local_sigs = [1, 1]
        for i in range(2):
            der = deriv_num_2(l_ext_ts, l_ext_vs, local_maxes[i])
            local_sigs[i] = np.sqrt(-1/der) if der < 0 else fwhm #avoid nans by'''
        if verbose > 2:
            print("\tdouble guess: peaks={}, widths={}".format([l_ext_ts[i] for i in local_maxes], local_sigs))
        #create an initial guess by supposing that the widths are each one half
        x0 = np.array([1.0, self.max_t, fwhm/2, 1.0, self.max_t+fwhm, fwhm/2])
        if local_maxes[0] != local_maxes[1]:
            x0 = np.array([l_ext_vs[local_maxes[0]], l_ext_ts[local_maxes[0]], local_sigs[0], l_ext_vs[local_maxes[1]], l_ext_ts[local_maxes[1]], local_sigs[1]])
        return x0, l_ext_ts, l_ext_vs

    def opt_double_envelope(self, fig_name=''):
        '''Optimize an envelope that is a superposition of two gaussians'''
        x0, l_ext_ts, l_ext_vs = self.get_double_guess()
        #the square error
        def doub_gauss_ser(x, t_pts):
            return x[0]*np.exp(-0.5*((t_pts - x[1])/x[2])**2) + x[3]*np.exp(-0.5*((t_pts - x[4])/x[5])**2)
        def ff(x):
            return np.sum( (np.abs(x[0])*np.exp(-0.5*((l_ext_ts - x[1])/x[2])**2) + np.abs(x[3])*np.exp(-0.5*((l_ext_ts - x[4])/x[5])**2) - l_ext_vs)**2 )
        def jac_env(x):
            #negative amplitudes are invalid, fix them to be positive
            x[0] = np.abs(x[0])
            x[3] = np.abs(x[3])
            t_dif_0 = (l_ext_ts - x[1])
            t_dif_1 = (l_ext_ts - x[4])
            exp_0 = np.exp(-0.5*(t_dif_0/x[2])**2)
            exp_1 = np.exp(-0.5*(t_dif_1/x[5])**2)
            diff_ser = x[0]*exp_0 + x[3]*exp_1 - l_ext_vs
            ret = np.zeros(6)
            ret[0] = 2*np.sum(diff_ser*exp_0)
            ret[1] = 2*x[0]*np.sum(diff_ser*t_dif_0*exp_0)/(x[2]**2)
            ret[2] = 2*x[0]*np.sum(diff_ser*exp_0*(t_dif_0)**2)/(x[2]**3)
            ret[3] = 2*np.sum(diff_ser*exp_1)
            ret[4] = 2*x[3]*np.sum(diff_ser*t_dif_1*exp_1)/(x[5]**2)
            ret[5] = 2*x[3]*np.sum(diff_ser*exp_1*(t_dif_1)**2)/(x[5]**3)
            return ret
        #perform the optimization and make figures
        res = opt.minimize(ff, x0)
        if fig_name != '':
            fig = plt.figure()
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            ax.scatter(l_ext_ts, l_ext_vs, color='black')
            t_range = l_ext_ts[-1] - l_ext_ts[0]
            t_cent = (l_ext_ts[-1] + l_ext_ts[0])/2
            int_ts = np.linspace(t_cent-t_range, t_cent+t_range, num=200)
            ax.plot(int_ts, np.abs(doub_gauss_ser(x0, int_ts)), color='red')
            ax.plot(int_ts, np.abs(doub_gauss_ser(res.x, int_ts)), color='blue')
            #plot the response
            fwhm = self.find_fwhm(both_tails=True)
            conv,_ = self.find_response(lambda t: np.exp( -t**2/(2*fwhm**2) ))
            ax.plot(l_ext_ts, conv/self.max_v)
            fig.savefig(fig_name+".png", dpi=300)
            plt.close(fig)

        res.x[0] = np.abs(res.x[0])*self.max_v
        res.x[3] = np.abs(res.x[3])*self.max_v
        #now create crude estimates of omega and phi by looking at the spacing between peaks and the phase offset between the data peak and the estimate for t_0
        est_omega = np.pi/self.typ_spacing
        est_phi = fix_angle(est_omega*(self.max_t - res.x[1]))
        if verbose > 3:
            print("\tdouble input: {}", x0)
            print("\tdouble optimized: {}", res.x)
        return res, est_omega, est_phi

    def sq_err_fit_single(self, x):
        return np.sum( (abs(x[0])*np.exp(-0.5*((self.ext_ts - x[1])/x[2])**2) - self.ext_vs)**2 )

    def sq_err_fit_double(self, x):
        return np.sum( (np.abs(x[0])*np.exp(-0.5*((self.ext_ts - x[1])/x[2])**2) + np.abs(x[3])*np.exp(-0.5*((self.ext_ts - x[4])/x[5])**2) - self.ext_vs)**2 )

def peaks_to_pulse(t_pts, peaks):
    res = np.zeros(len(t_pts))
    for p in peaks:
        res += p[2]*np.exp( -((t_pts-p[0])/p[1])**2 )
    return res

def find_peaks(t_pts, vs, ns, width_steps=5):
    '''Helper function for est_env_new
    t_pts: the sampled points in time, this should have the same dimension as vs and ns
    vs: the absolute value of the actual pulse. This WILL break if there are any negative values.
    ns: the envelope of the pulse obtained through some signal processing
    '''
    dt = t_pts[1] - t_pts[0]
    #now that we have an envolope, find its maxima
    pulse_peaks = ssig.argrelmax(vs)[0]
    env_peaks = ssig.argrelmax(ns, order=4)[0]
    max_peak = np.max(ns.take(env_peaks))
    #return a single default peak if no extrema were found to prevent out of bounds errors
    if len(pulse_peaks) == 0 or len(env_peaks) == 0:
        return np.array([[0, 1, 1]])
    
    #an array of peaks. The columns are time, width, amplitude, and the residual error. The error is cumulative i.e. The pulse incorporates all peaks before the current one and computes the square error.
    peaks_arr = np.zeros((len(env_peaks),4))
    i = 0
    for v in env_peaks:
        now = t_pts[v]
        #find the nearest maxima in the actual pulse. The envelopes aren't great at getting amplitudes, but they are useful for widths and centers.
        p_ind = 0
        closest_sep = t_pts[-1] - t_pts[0]
        for p in pulse_peaks:
            sep = abs(t_pts[p] - now)
            if sep < closest_sep:
                p_ind = p
                closest_sep = sep
            else:
                break
        if vs[p_ind] > PEAK_AMP_THRESH*max_peak:
            #make sure that we're inside the array but we take at least one step
            stp = max(min(min(width_steps, v), len(t_pts)-v-1), 1) 
            #do a quadratic fit to get a guess of the Gaussian width
            samp_ts = t_pts[v-stp:v+stp] - t_pts[v]
            samp_ls = np.log(max_peak)-np.log(ns[v-stp:v+stp])
            reg,cov = opt.curve_fit(lambda x,a,b: a*x**2 + b, samp_ts, samp_ls)
            peaks_arr[i, 1] = 1/np.sqrt(reg[0])
            #only add this to the array if the value isn't nan
            if not np.isnan(peaks_arr[i, 1]):
                peaks_arr[i, 0] = now
                peaks_arr[i, 2] = vs[p_ind]
                peaks_arr[i, 3] = np.sum((peaks_to_pulse(t_pts[pulse_peaks], peaks_arr[0:i+1]) - vs[pulse_peaks])**2)
                i += 1
    #sort peaks by intensity
    peaks_arr = peaks_arr[(-peaks_arr)[:, 2].argsort()]
    peaks_arr.resize((i,4))
    return peaks_arr

def est_env_new(t_pts, v_pts, fig_name='', lowpass_inc=1.0, max_n_peaks=4, in_freq=1.3, scan_length=5):
    '''Based on a time series sampled at t_pts with values v_pts, try to find a Gaussian pulse envelope.'''
    #take the fourier transform
    freqs = fft.fftfreq(len(v_pts), d=t_pts[1]-t_pts[0])
    vf0 = 2*fft.fft(v_pts)
    vf = np.copy(vf0)
    avg_len = 0.0
    while True:
        phi_est = np.angle(vf[0])
        mags = np.abs(vf*np.heaviside(freqs, 1))
        #Shift the pulse so that the peak magnitude is at zero frequency. Taking the inverse fourier transform of this shifted magnitude gives pulse amplitude. The frequency at which the magnitude is maximized is the same as the unmodulated sine wave.
        f0_ind = np.argmax(mags)
        #comp = lambda x1,x2: True if abs(x1) > abs(x2) else False
        v_abs = np.abs(v_pts)
        t0_ind = np.argmax(v_abs)
        #take the inverse fourier transform and shift to the peak of the pulse
        env_pts = np.roll( np.abs(fft.ifft(np.roll(mags, -f0_ind))), t0_ind )
        peaks_arr = find_peaks(t_pts, v_abs, env_pts)
        #if there are many peaks that is a sign that this is a noisy signal. We apply a lowpass filter in this case. Default to 2*peak_magnitude to ensure that a lowpass filter is applied in the case where the array is empty.
        noisiness = np.max(mags[4*f0_ind-scan_length:4*f0_ind+scan_length]) if 4*f0_ind-scan_length < len(mags) else 2*mags[f0_ind]
        if noisiness < NOISE_THRESH*mags[f0_ind] or avg_len >= MAX_LOWPASS_EVALS*lowpass_inc:
            break
        else:
            avg_len += lowpass_inc
            vf = vf0*np.sinc(freqs*avg_len)
            v_pts = np.real(fft.ifft(vf))
            if verbose > 0:
                print("\tNoisy data, applying low pass filter strength={}".format(avg_len))

    #make plots if requested
    if fig_name != '':
        if fig_name != '':
            fig, ax = plt.subplots(2)
            #plot time domain stuff
            ax[0].set_title("Time domain")
            ax[0].plot(t_pts, v_pts, color='black')
            for i in range(len(peaks_arr)):
                ax[0].plot(t_pts, peaks_to_pulse(t_pts, peaks_arr[i:i+1]), color='gray')
            ax[0].plot(t_pts, peaks_to_pulse(t_pts, peaks_arr), color='green', linestyle=':')
            #plot frequency domain stuff
            ax[1].set_title("Frequency domain")
            ax[1].plot(freqs, vf, color='black')
            ax[1].plot(freqs, mags, color='gray', linestyle=':')
            fig.savefig(fig_name+".png", dpi=300)
            plt.close(fig)
 
    #We want to return the actual time series that was used for analysis. So we check if any lowpass filters were applied.
    return v_pts, peaks_arr, 2*np.pi*freqs[f0_ind], phi_est

def est_env_new_wip(t_pts, v_pts, fig_name='', lowpass_strength=1.0, max_n_peaks=4):
    '''Based on a time series sampled at t_pts with values v_pts, try to find a Gaussian pulse envelope.'''
    freqs = fft.fftfreq(len(v_pts), d=t_pts[1]-t_pts[0])
    vf = 2*fft.fft(v_pts)

    #if there is a large number of identified peaks, that is a sign that the data is very noisy. We repeatedly apply a lowpass filter until we are below the maximum allowed number of peaks.
    while True:
        phi_est = np.angle(vf[0])
        mags = np.abs(vf*np.heaviside(freqs, 1))
        v_abs = np.abs(v_pts)
        f0_ind = np.argmax(v_abs)
        t0_ind = np.argmax(mags)
        env_pts = np.roll(np.abs( fft.ifft(np.roll(mags, -f0_ind)) ), t0_ind)
        #find the peaks and check if we need to go again
        peaks_arr = find_peaks(t_pts, v_abs, env_pts)
        if peaks_arr.shape[0] <= max_n_peaks:
            break
        else:
            vf = vf*np.sinc(freqs*lowpass_strength)
            v_pts = np.real(fft.ifft(vf))
            if verbose > 2:
                print("\tApplying lowpass filter")

    #make plots if requested
    if fig_name != '':
        if fig_name != '':
            fig, ax = plt.subplots(2)
            #plot time domain stuff
            ax[0].set_title("Time domain")
            ax[0].plot(t_pts, v_pts, color='black')
            for i in range(len(peaks_arr)):
                ax[0].plot(t_pts, peaks_to_pulse(t_pts, peaks_arr[i:i+1]), color='gray')
            ax[0].plot(t_pts, peaks_to_pulse(t_pts, peaks_arr), color='green', linestyle=':')
            #plot frequency domain stuff
            ax[1].set_title("Frequency domain")
            ax[1].plot(freqs, vf, color='black')
            ax[1].plot(freqs, mags, color='gray', linestyle=':')
            fig.savefig(fig_name+".png", dpi=300)
            plt.close(fig)
    return env_pts, peaks_arr, 2*np.pi*freqs[f0_ind], phi_est

#all of this is for debuging, first pick a cluster then get the double guess
def doub_gauss_ser(x, t_pts):
    return x[0]*np.exp(-0.5*((t_pts - x[1])/x[2])**2) + x[3]*np.exp(-0.5*((t_pts - x[4])/x[5])**2)

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

def opt_pulse(t_pts, a_pts, env_x, est_omega, est_phi, keep_n=DEF_KEEP_N):
    '''Set values for the a_pts for each time point
    a_pts: the values of the electric field at each corresponding time point in t_pts. Note that len(a_pts) must be the same as len(t_pts)
    env_x: an array representing the value found from an envelope fit (i.e if you call res=opt_envelope(...) this would be res.x
    est_omega: a guess of the value of omega based on the envelope fitting
    est_phi: a guess of the value of the CEP based on the envelope fitting
    returns: the fit result along with the lower and upper times used to perform the fitting'''

    #truncate so that we only look at N_STDS deviations from the center of the pulse
    if keep_n > 0:
        dt = (t_pts[-1]-t_pts[0])/t_pts.shape[0]
        i_0 = int(env_x[1]/dt)
        i_shift = int(keep_n*env_x[2] / dt)
        #make sure we don't go past the ends of the array
        if i_0 < i_shift:
            i_shift = i_0
        if i_0+i_shift > t_pts.shape[0]:
            i_shift = min(t_pts.shape[0] - i_0 - 1, i_0)
        t_pts = t_pts[i_0-i_shift:i_0+i_shift]
        a_pts = a_pts[i_0-i_shift:i_0+i_shift]
        if len(t_pts) == 0:
            t_pts = t_pts
            a_pts = a_pts
    else:
        t_pts = t_pts
        a_pts = a_pts

    #it is helpful to try to get the amplitude close to unity so that the gradient in other directions is not suppressed
    if env_x[0] != 0.0:
        a_pts = a_pts / env_x[0]

    #square error of the pulse envelope
    def fp(x):
        if x[2] == 0:
            x[2] = 0.1
        return np.sum( (gauss_series(x, t_pts) - a_pts)**2 )
    def jac_fp(x):
        p_i = gauss_series(x, t_pts) - a_pts
        exp_cos = np.exp(-0.5*((t_pts-x[1])/x[2])**2)*np.cos(x[3]*(t_pts-x[1])-x[4])
        exp_sin = np.exp(-0.5*((t_pts-x[1])/x[2])**2)*np.sin(x[3]*(t_pts-x[1])-x[4])       
        t_dif_by_w = (t_pts-x[1])
        ret = np.zeros(5)
        ret[0] = 2*np.sum(p_i*exp_cos)
        ret[1] = 2*x[0]*np.sum(p_i*(t_dif_by_w*exp_cos/(x[2]**2) + x[3]*exp_sin))
        ret[2] = 2*x[0]*np.sum(p_i*exp_cos*t_dif_by_w**2)/(x[2]**3)
        ret[3] = -2*x[0]*np.sum(p_i*exp_sin*t_dif_by_w)
        ret[4] = 2*x[0]*np.sum(p_i*exp_sin)
        return np.real(ret)

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
    #ensure that the pulses are ordeblue the way we expect, fix the order such that the biggest pulse comes first
    if res.x.shape[0] == 8 and res.x[0] < res.x[5]:
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
        #negative amplitudes are invalid, fix them to be positive
        x[0] = np.abs(x[0])
        x[5] = np.abs(x[5])
        #the difference between the fit/data series
        p_i = double_gauss_series(x, t_pts) - a_pts
        t_dif_0 = (t_pts-x[1])
        t_dif_1 = (t_pts-x[6])
        exp_cos_0 = x[0]*np.exp(-0.5*(t_dif_0/x[2])**2)*np.cos(x[3]*t_dif_0-x[4])
        exp_sin_0 = x[0]*np.exp(-0.5*(t_dif_0/x[2])**2)*np.sin(x[3]*t_dif_0-x[4])
        exp_cos_1 = x[5]*np.exp(-0.5*(t_dif_1/x[7])**2)*np.cos(x[3]*t_dif_0-x[4])
        exp_sin_1 = x[5]*np.exp(-0.5*(t_dif_1/x[7])**2)*np.sin(x[3]*t_dif_0-x[4])
        ret = np.zeros(8)
        ret[0] = 2*np.sum(p_i*exp_cos_0)/x[0]
        ret[1] = 2*np.sum(p_i*(t_dif_0*exp_cos_0/(x[2]**2) + x[3]*(exp_sin_0+exp_sin_1)))
        ret[2] = 2*np.sum(p_i*exp_cos_0*t_dif_0**2)/(x[2]**3)
        ret[3] =-2*np.sum(p_i*t_dif_0*(exp_sin_0+exp_sin_1))
        ret[4] = 2*np.sum(p_i*(exp_sin_0+exp_sin_1))
        ret[5] = 2*np.sum(p_i*exp_cos_1)/x[5]
        ret[6] = 2*np.sum(p_i*(t_dif_1*exp_cos_1/(x[7]**2)))
        ret[7] = 2*np.sum(p_i*exp_cos_1*t_dif_1**2)/(x[7]**3)
        return np.real(ret)

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

    #return fix_double_pulse_order(res), t_pts[0], t_pts[-1]
    return res, t_pts[0], t_pts[-1]

def opt_pulse_full_old(t_pts, a_pts, omega_fft=-1, a_sigmas_sq=0.1, keep_n=DEF_KEEP_N, fig_name=''):
    if a_pts.shape != t_pts.shape:
        raise ValueError("t_pts and a_pts must have the same shape")
    if t_pts.shape[0] == 0:
        raise ValueError("empty time series supplied!")
    env_fig_name = ''
    env_fig_name_2 = ''
    if fig_name != '':
        env_fig_name = fig_name+"_env"
        env_fig_name_2 = fig_name+"_env_2"

    #the maximum of the derivative on a_pts(t_pts) and use this to get a conservative estimate of the variance on each of the sample points. Note that for an RV Y with P(Y|X)=N(mX+x0, s) we have V(Y) = s^2 + m^2 V(X). We take X to follow a continuous uniform distribution between (t_pts[i], t_pts[i+SKIP])
    max_diff = np.max(np.diff(np.abs(a_pts)))
    skip_sigma_sq = a_sigmas_sq + (max_diff**2)/12
    #find the extrema
    env_fit = EnvelopeFitter(t_pts, a_pts)
    #set defaults
    env_res = None
    res = None
    low_t = t_pts[0]
    hi_t = t_pts[-1]
    err = 1.0
    used_double = False
    #try fitting to both possible envelopes
    env_res_0, est_omega_0, est_phi_0 = env_fit.opt_envelope(keep_n=keep_n, fig_name=env_fig_name)
    env_res_1, est_omega_1, est_phi_1 = env_fit.opt_double_envelope(fig_name=env_fig_name_2) 
    err_0 = env_fit.sq_err_fit_single(env_res_0.x)
    err_1 = env_fit.sq_err_fit_double(env_res_1.x)
    #compute the bayes factor for the single and double pulses. We take both models to have identical priors. Note that under assumption of Gaussian errors the error cancels
    ln_bayes =  (0.5/skip_sigma_sq)*(err_0 - err_1)
    if verbose > 1:
        print("\tenvelope sigma={}, bayes factor={}".format(skip_sigma_sq, np.exp(ln_bayes)))
    #substantial evidence as defined by Jeffreys
    if np.isnan(ln_bayes) or ln_bayes > 1.15:
        env_res = env_res_1
        if omega_fft <= 0:
            omega_fft = est_omega_1
        res, low_t, hi_t = opt_double_pulse(t_pts, a_pts, env_res_1.x, omega_fft, est_phi_1, keep_n=keep_n)
        err = err_1
        used_double = True
    else:
        env_res = env_res_0
        if omega_fft <= 0:
            omega_fft = est_omega_0
        res, low_t, hi_t = opt_pulse(t_pts, a_pts, env_res_0.x, omega_fft, est_phi_0, keep_n=keep_n)
        err = err_0
        used_double = False
    if verbose > 0:
        print("\tenvelope type={}".format("double" if used_double else "single"))

    #plot figures if requested
    if fig_name != '':
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.plot(t_pts, a_pts, color='black')
        ax.plot(t_pts, gauss_env(env_res.x, t_pts), linestyle=':', color='red')
        ax.plot(t_pts, gauss_env(res.x, t_pts), linestyle=':', color='blue')
        if used_double:
            x0 = np.array([env_res.x[0], env_res.x[1], env_res.x[2], est_omega_1, est_phi_1, env_res.x[3], env_res.x[4], env_res.x[5]])
            ax.plot(t_pts, double_gauss_series(x0, t_pts), color='red')
            ax.plot(t_pts, double_gauss_series(res.x, t_pts), color='blue')
        else:
            x0 = np.array([env_res.x[0], env_res.x[1], env_res.x[2], est_omega_0, est_phi_0])
            ax.plot(t_pts, gauss_series(x0, t_pts), color='red')
            ax.plot(t_pts, gauss_series(res.x, t_pts), color='blue')
        '''ax.vlines([low_t, hi_t], -1.2, 1.2)
        ax.set_ylim((-1.2, 1.2))'''
        ax.annotate(r"$er^2={:.1f}$".format(res.fun), (0.01, 0.84), xycoords='axes fraction')
        fig.savefig(fig_name+".png", dpi=300)
        plt.close(fig)

    if fig_name != '':
        res_name = fig_name+"_res.txt"
        write_reses(res_name, [res,env_res])

    return res, env_res

def opt_pulse_full(t_pts, a_pts, err_sq, lowpass_inc=1.0, keep_n=DEF_KEEP_N, fig_name='', in_freq=1.3):
    if a_pts.shape != t_pts.shape:
        raise ValueError("t_pts and a_pts must have the same shape")
    if t_pts.shape[0] == 0:
        raise ValueError("empty time series supplied!")
    env_fig_name = ''
    if fig_name != '':
        env_fig_name = fig_name+"_env"

    #do the signal processing to find the envelope and decompose it into a series of peaks
    a_pts, peak_arr, f0, phi = est_env_new(t_pts, a_pts, fig_name=env_fig_name, lowpass_inc=lowpass_inc, in_freq=in_freq)

    if verbose > 5:
        print(peak_arr)

    used_double = False
    res = None
    low_t = t_pts[0]
    hi_t = t_pts[-1]

    #substantial evidence as defined by Jeffreys
    ln_bayes = (0.5/err_sq)*(peak_arr[0,3] - peak_arr[1,3]) if peak_arr.shape[0] > 1 else 0.0
    if verbose > 1:
        print("\tenvelope sigma={}, bayes factor={}".format(err_sq, np.exp(ln_bayes)))
    if np.isnan(ln_bayes) or ln_bayes > 1.15:
        guess_x = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1], peak_arr[1,2], peak_arr[1,0], peak_arr[1,1]])
        res, low_t, hi_t = opt_double_pulse(t_pts, a_pts, guess_x, f0, phi, keep_n=keep_n)
        used_double = True
    else:
        guess_x = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1]])
        res, low_t, hi_t = opt_pulse(t_pts, a_pts, guess_x, f0, phi, keep_n=keep_n)

    if verbose > 0:
        print("\tenvelope type={}".format("double" if used_double else "single"))

    #plot figures if requested
    if fig_name != '':
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.plot(t_pts, a_pts, color='black')
        ax.plot(t_pts, gauss_env(guess_x, t_pts), linestyle=':', color='red')
        ax.plot(t_pts, gauss_env(res.x, t_pts), linestyle=':', color='blue')
        if used_double:
            x0 = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1], f0, phi, peak_arr[1,2], peak_arr[1,0], peak_arr[1,1]])
            ax.plot(t_pts, double_gauss_series(x0, t_pts), color='red')
            ax.plot(t_pts, double_gauss_series(res.x, t_pts), color='blue')
        else:
            x0 = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1], f0, phi])
            ax.plot(t_pts, gauss_series(x0, t_pts), color='red')
            ax.plot(t_pts, gauss_series(res.x, t_pts), color='blue')
        '''ax.vlines([low_t, hi_t], -1.2, 1.2)
        ax.set_ylim((-1.2, 1.2))'''
        ax.annotate(r"$er^2={:.1f}$".format(res.fun), (0.01, 0.84), xycoords='axes fraction')
        fig.savefig(fig_name+".png", dpi=300)
        plt.close(fig)

    if fig_name != '':
        res_name = fig_name+"_res.txt"
        write_reses(res_name, [res])

    return res, None

def get_params(res):
    '''Convert a scipy optimize result returned from opt_pulse_full into a human readable set of parameters
    '''
    amp = res.x[0]
    t_0 = res.x[1]
    sig = np.sqrt(res.x[2]/2)
    omega = res.x[3]
    cep = res.x[4]
    if amp < 0:
        amp *= -1
        if cep > 0:
            cep -= np.pi
        else:
            cep += np.pi
    cep = fix_angle(cep)
    cepr = cep
    if len(res.x) > 6:
        cepr = cep + omega*(res.x[5] - res.x[1])
    return amp, t_0, sig, omega, cep, cepr

N_RES_PARAMS = 6
class cluster_res:
    def __init__(self, xs):
        self.n_pts = len(xs)
        self.xs = np.array(xs)
        self.res_arr = np.zeros((2*N_RES_PARAMS + 1, self.n_pts))

    def fix_xs(self, x_0, scale):
        '''transforms all x points from x to x' where x'=(x-x_0)/scale
        '''
        self.xs = (self.xs-x_0)/scale

    def get_amp(self):
        return self.res_arr[0]
    def get_amp_err(self):
        return self.res_arr[1]
    def get_t0(self):
        return self.res_arr[2]
    def get_t0_err(self):
        return self.res_arr[3]
    def get_sig(self):
        return self.res_arr[4]
    def get_sig_err(self):
        return self.res_arr[5]
    def get_omega(self):
        return self.res_arr[6]
    def get_omega_err(self):
        return self.res_arr[7]
    def get_phase(self):
        return self.res_arr[8]
    def get_phase_err(self):
        return self.res_arr[9]
    def get_phase_ref(self):
        return self.res_arr[10]
    def get_phase_ref_err(self):
        return self.res_arr[11]
    def get_err_sq(self):
        return self.res_arr[12]

    def set_point(self, jj, res, err_2):
        '''set the point at index jj to have the paramters from the scipy optimization result res
        '''
        amp, t0, sig, omega, cep, cepr = get_params(res)
        self.res_arr[0,jj] = amp
        self.res_arr[1,jj] = np.sqrt(res.hess_inv[0][0]/(err_2))/2
        self.res_arr[2,jj] = t0
        self.res_arr[3,jj] = np.sqrt(res.hess_inv[1][1]/(err_2))
        self.res_arr[4,jj] = sig
        self.res_arr[5,jj] = np.sqrt(res.hess_inv[2][2]/(err_2*8*sig))
        self.res_arr[6,jj] = omega
        self.res_arr[7,jj] = np.sqrt(res.hess_inv[3][3]/(err_2))
        self.res_arr[8,jj] = cep/np.pi
        self.res_arr[9,jj] = np.sqrt(res.hess_inv[4][4]/(err_2*np.pi))
        self.res_arr[10,jj] = cepr/np.pi
        self.res_arr[11,jj] = np.sqrt(res.hess_inv[4][4]/(err_2*np.pi))
        self.res_arr[12,jj] = res.fun

    def trim_to(self, trim_arr):
        '''trim the result so that it only contains points from indices specified by trim_arr
        '''
        self.xs = np.take(self.xs, trim_arr)
        self.res_arr = np.take(self.res_arr, trim_arr, axis=1)

    def mirror_pts(self, x_0=0, scale=0):
        '''returns: a copy of the cluster_res that has twice as many points mirrored about the x=0 plane.
        x_0: if scale is not zero, then fix_xs(x_0, scale) is applied to the resulting data
        scale: if scale is not zero, then fix_xs(x_0, scale) is applied to the resulting data 
        '''
        nr = cluster_res(self.xs)
        if scale != 0:
            nr.fix_xs(x_0, scale)

        #we want to use spatial symmetry to infer the points which have not been collected
        i_zer = nr.n_pts
        for ii, zz in enumerate(nr.xs):
            if zz >= 0:
                i_zer = ii
                break
        #keep track of the last point that satisfies z<=0
        i_cent = i_zer-1
        new_size = 2*i_zer
        #if zero is included then this is a special case where there is an odd number of points
        if i_zer < nr.n_pts and nr.xs[i_zer] == 0:
            i_cent = i_zer
            new_size += 1
        if verbose > 1:
            print("new_size={}, i_zer={}, i_cent={}".format(new_size, i_zer, i_cent))
        if new_size > nr.n_pts:
            new_xs = np.zeros(new_size)
            nr.res_arr = np.pad(self.res_arr, ((0,0),(0,new_size-nr.n_pts)))
            for ii in range(nr.n_pts - i_zer):
                #update errors
                nr.res_arr[2*N_RES_PARAMS, i_cent-ii] = (nr.res_arr[2*N_RES_PARAMS, i_cent-ii] + nr.res_arr[2*N_RES_PARAMS, i_cent+ii])/2
                #update everything else
                for k in range(N_RES_PARAMS):
                    vid = 2*k
                    eid = 2*k+1
                    nr.res_arr[vid, i_cent-ii] = (nr.res_arr[vid,i_cent-ii] + nr.res_arr[vid,i_zer+ii])/2
                    nr.res_arr[eid, i_cent-ii] = np.sqrt(nr.res_arr[eid,i_cent-ii]**2 + nr.res_arr[eid,i_zer+ii]**2)
            for ii in range(new_size - i_zer):
                #update errors
                nr.res_arr[2*N_RES_PARAMS, i_cent+ii] = nr.res_arr[2*N_RES_PARAMS, i_cent-ii]
                #update everything else
                new_xs[i_cent-ii] = nr.xs[i_cent-ii]
                new_xs[i_zer+ii] = -nr.xs[i_cent-ii]
                for k in range(N_RES_PARAMS):
                    vid = 2*k
                    eid = 2*k+1
                    nr.res_arr[vid,i_cent+ii] = nr.res_arr[vid, i_cent-ii]
                    nr.res_arr[eid,i_cent+ii] = nr.res_arr[eid, i_cent-ii]
        else:
            new_xs = np.array(nr.xs[:new_size])
            nr.res_arr = np.resize(self.res_arr, (2*N_RES_PARAMS, new_size))
        nr.n_pts = new_size
        nr.xs = new_xs
        return nr

class phase_finder:
    '''
    Initialize a phase finder reading from the  specified by fname and a juction width and gap specified by width and height respectively.
    fname: h5 file to read time samples from
    width: the width of the junction
    height: the thickness of the junction
    pass_alpha: The scale of the low pass filter applied to the time series. This corresponds to taking a time average of duration 2/pass_alpha on either side
    '''
    def __init__(self, fname, width, height, pass_alpha=1.0, slice_dir='x', prefix='.'):
        self.prefix = prefix
        self.slice_dir = slice_dir
        self.slice_name = 'z'
        if self.slice_dir == 'z':
            self.slice_name = 'x'
        #open the h5 file and identify all the clusters
        self.f = h5py.File(fname, "r")
        self.clust_names = []
        for key in self.f.keys():
            #make sure that the cluster has a valid name and it actually has points
            if 'cluster' in key and len(self.f[key]) > 1 and len(self.f[key]['locations']) > 0:
                self.clust_names.append(key)
        self.n_clusts = len(self.clust_names)
        #create a numpy array for time points, this will be shared across all points
        keylist = list(self.f[self.clust_names[0]].keys())
        self.n_t_pts = len(self.f[self.clust_names[0]][keylist[1]]['time'])
        #convert wavelength in um to frequency in 1/fs
        self.in_freq = .299792458 / self.f['info']['sources']['wavelen']
        #read information to figure out time units and step sizes
        t_min = self.f['info']['time_bounds'][0]
        t_max = self.f['info']['time_bounds'][1]
        self.dt = (t_max-t_min)/self.n_t_pts
        self.t_pts = np.linspace(t_min, t_max, num=self.n_t_pts)
        #these are used for applying a low pass filter to the time data, the low-pass is a sinc function applied in frequency space
        self.lowpass_inc = pass_alpha
        self.f_pts = fft.fftfreq(self.n_t_pts, d=self.dt)
        self.low_filter = np.sin(self.f_pts*np.pi*pass_alpha) / (self.f_pts*np.pi*pass_alpha)
        self.low_filter[0] = 1
        #this normalizes the square errors to have reasonable units
        self.sq_er_fact = self.dt/(t_max-t_min)

    def get_src_wavelen(self):
        return self.f['info']['sources']['wavelen']
    def get_src_width(self):
        return self.f['info']['sources']['width']
    def get_src_phase(self):
        return self.f['info']['sources']['phase']
    def get_src_start_time(self):
        return self.f['info']['sources']['start_time']
    def get_src_end_time(self):
        return self.f['info']['sources']['end_time']
    def get_src_amplitude(self):
        return self.f['info']['sources']['amplitude']

    def get_clust_span(self, clust):
        if isinstance(clust, int):
            clust = self.clust_names[0]
        return self.f[clust]['locations'][self.slice_dir][0], self.f[clust]['locations'][self.slice_dir][-1]
    def get_clust_location(self, clust):
        return self.f[clust]['locations'][self.slice_name][0]

    def get_point_times(self, clust, ind, low_pass=False):
        '''#fetch a list of points and their associated coordinates
        points = list(self.f[clust].keys())[1:]
        v_pts = np.array(self.f[clust][points[ind]]['time']['Re']) + 1j*np.array(self.f[clust][points[ind]]['time']['Im'])
        # (dt*|dE/dt|)^2 evaluated at t=t_max (one factor of two comes from the imaginary part, the other from the gaussian. We do something incredibly hacky to account for the spatial part. We assume that it has an error identical to the time error, and add the two in quadrature hence the other factor of two.
        #err_2 = 8*np.max(np.diff(v_pts)**2)
        err_2 = 0.02
        vf_pts = fft.fft(v_pts)
        if low_pass:
            vf_pts = fft.fft(v_pts)*self.low_filter
            v_pts = fft.ifft(vf_pts)
        peak_i = np.argmax(np.abs(vf_pts))
        peak_omega = 2*np.pi*self.f_pts[peak_i]
        return v_pts, np.abs(peak_omega), err_2'''
        points = list(self.f[clust].keys())[1:]
        #TODO: generate an estimate of err_2
        err_2 = 0.02
        return np.array(self.f[clust][points[ind]]['time']['Re']), err_2

    def read_cluster(self, clust, save_fit_figs=False):
        '''Read an h5 file and return three two dimensional numpy arrays of the form ([amplitudes, errors], [sigmas, errors], [phases, errors]
        '''
        #for performance metrics
        n_evals = 0
        t_start = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
        #fetch a list of points and their associated coordinates
        points = list(self.f[clust].keys())[1:]
        xs = np.array(self.f[clust]['locations'][self.slice_dir])
        n_x_pts = len(points)
        if xs.shape[0] != n_x_pts:
            raise ValueError("points and locations do not have the same size")

        ret = cluster_res(xs)
        good_js = []
        phase_cor = 0
        for j in range(len(points)):
            v_pts, err_2 = self.get_point_times(clust, j, low_pass=False)
            if verbose > 0:
                #print("clust={}, j={}\n\tfield error^2={}\n\tomega_fft={}".format(clust, j, err_2*self.n_t_pts, guess_omega))
                print("clust={}, j={}\n\tfield error^2={}\n".format(clust, j, err_2*self.n_t_pts))
            #before doing anything else, save a plot of just the time series
            fig_name = ""
            if save_fit_figs:
                plt.plot(self.t_pts, np.real(v_pts))
                plt.savefig("{}/fit_figs/t_series_{}_{}.pdf".format(self.prefix,clust,j))
                plt.clf()
                fig_name = "{}/fit_figs/fit_{}_{}".format(self.prefix,clust,j)
            good_fit = True
            #now actually try optimizing
            res,_ = opt_pulse_full(self.t_pts, np.real(v_pts), err_2, keep_n=2.5, fig_name=fig_name, lowpass_inc=self.lowpass_inc, in_freq=self.in_freq)
            try:
                #res,res_env = opt_pulse_full(self.t_pts, np.real(v_pts), a_sigmas_sq=err_2, keep_n=2.5, fig_name=fig_name, omega_fft=guess_omega)
                res,_ = opt_pulse_full(self.t_pts, np.real(v_pts), err_2, keep_n=2.5, fig_name=fig_name, lowpass_inc=self.lowpass_inc)
            except:
                if verbose > 0:
                    print("\tfitting to raw time series failed, applying low pass filter")
                if verbose > 0:
                    print("clust={}, j={}\n\tfield error^2={}\n\tomega_fft={}".format(clust, j, err_2*self.n_t_pts, guess_omega))
                v_pts, guess_omega, err_2 = self.get_point_times(clust, j, low_pass=True)
                #res,res_env = opt_pulse_full_old(self.t_pts, np.real(v_pts), a_sigmas_sq=err_2, keep_n=2.5, fig_name=fig_name, omega_fft=guess_omega)
                res,_ = opt_pulse_full(self.t_pts, np.real(v_pts), err_2, keep_n=2.5, fig_name=fig_name, omega_fft=guess_omega)
                '''try:
                    print(guess_omega)
                    res,res_env = opt_pulse_full(self.t_pts, np.real(v_pts), a_sigmas_sq=err_2, keep_n=2.5, fig_name=fig_name, omega_fft=guess_omega)
                except:
                    good_fit = False'''
            if good_fit:
                n_evals += 1
                if verbose > 0:
                    #print("\tsquare errors = {}, {}\n\tx={}\n\tdiag(H^-1)={}".format(res.fun, res_env.fun, res.x, np.diagonal(res.hess_inv)))
                    print("\tsquare errors = {}\n\tx={}\n\tdiag(H^-1)={}".format(res.fun, res.x, np.diagonal(res.hess_inv)))
                #only include this point if the fit was good
                if res.fun*self.sq_er_fact < 50.0:
                    good_js.append(j)
                    ret.set_point(j, res, err_2)
            elif verbose > 0:
                print("\tbad fit!".format(clust,j))
        ret.trim_to(good_js)
        if verbose > 1:
            print("\tfound {} good points".format(len(good_js), clust))
        #figure out the time required for calculations
        t_dif = time.clock_gettime_ns(time.CLOCK_MONOTONIC) - t_start
        if n_evals > 0:
            print("Completed optimizations in {:.5E} ns, average time per eval: {:.5E} ns".format(t_dif, t_dif/n_evals))
        else:
            print("No successful optimizations!")
        return ret

    def lookup_fits(self, clust_name, recompute=False, save_fit_figs=False):
        '''Load the fits from the pickle file located at fname or perform the fits if it doesn't exist. Return the results'''
        data_name = '{}/dat_{}'.format(self.prefix, clust_name)
        if recompute or not os.path.exists(data_name):
            #figure out data by phase fitting
            ret = self.read_cluster(clust_name, save_fit_figs=save_fit_figs)
            with open(data_name, 'wb') as fh:
                pickle.dump(ret, fh)
            return ret
        else:
            with open(data_name, 'rb') as fh:
                ret = pickle.load(fh)
            return ret

    def find_avg_vals(self, clust):
        ret = self.lookup_fits(clust)
        field_0 = np.sum(ret.get_amps()) / ret.n_pts
        sig_0 = np.sum(ret.get_sigs()) / ret.n_pts
        phase_0 = np.sum(ret.get_phases()) / ret.n_pts
        return field_0, sig_0, phase_0
