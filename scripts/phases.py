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
        amp = res.x[0]
        t0 = res.x[1]
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
    def __init__(self, fname, width, height, pass_alpha=1.0, slice_dir='x', prefix='.', keep_n=2.5, scan_length=5):
        self.prefix = prefix
        self.slice_dir = slice_dir
        self.slice_name = 'z'
        if self.slice_dir == 'z':
            self.slice_name = 'x'
        self.keep_n = keep_n
        self.scan_length = scan_length
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
        #fetch a list of points and their associated coordinates
        points = list(self.f[clust].keys())[1:]
        #TODO: generate an estimate of err_2
        #err_2 = 8*np.max(np.diff(v_pts)**2)
        err_2 = 0.02
        return np.array(self.f[clust][points[ind]]['time']['Re']), err_2

    def est_env_new(self, t_pts, v_pts, fig_name=''):
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
            v_abs = np.abs(v_pts)
            t0_ind = np.argmax(v_abs)
            #take the inverse fourier transform and shift to the peak of the pulse
            env_pts = np.roll( np.abs(fft.ifft(np.roll(mags, -f0_ind))), t0_ind )
            peaks_arr = find_peaks(t_pts, v_abs, env_pts)
            #if there are many peaks that is a sign that this is a noisy signal. We apply a lowpass filter in this case. Default to 2*peak_magnitude to ensure that a lowpass filter is applied in the case where the array is empty.
            start = 4*f0_ind-self.scan_length
            end = 4*f0_ind+self.scan_length
            noisiness = np.max(mags[start:end]) if start < len(mags) else 2*mags[f0_ind]
            if noisiness < NOISE_THRESH*mags[f0_ind] or avg_len >= MAX_LOWPASS_EVALS*self.lowpass_inc:
                break
            else:
                avg_len += self.lowpass_inc
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

    def opt_pulse(self, t_pts, a_pts, env_x, est_omega, est_phi):
        '''Set values for the a_pts for each time point
        a_pts: the values of the electric field at each corresponding time point in t_pts. Note that len(a_pts) must be the same as len(t_pts)
        env_x: an array representing the value found from an envelope fit (i.e if you call res=opt_envelope(...) this would be res.x
        est_omega: a guess of the value of omega based on the envelope fitting
        est_phi: a guess of the value of the CEP based on the envelope fitting
        returns: the fit result along with the lower and upper times used to perform the fitting'''

        #truncate so that we only look at N_STDS deviations from the center of the pulse
        if self.keep_n > 0:
            dt = (t_pts[-1]-t_pts[0])/t_pts.shape[0]
            i_0 = int(env_x[1]/dt)
            i_shift = int(self.keep_n*env_x[2] / dt)
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
    def opt_double_pulse(self, t_pts, a_pts, env_x, est_omega, est_phi):
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

    def opt_pulse_full(self, t_pts, a_pts, err_sq, fig_name=''):
        if a_pts.shape != t_pts.shape:
            raise ValueError("t_pts and a_pts must have the same shape")
        if t_pts.shape[0] == 0:
            raise ValueError("empty time series supplied!")
        env_fig_name = ''
        if fig_name != '':
            env_fig_name = fig_name+"_env"

        #do the signal processing to find the envelope and decompose it into a series of peaks
        a_pts, peak_arr, f0, phi = self.est_env_new(t_pts, a_pts, fig_name=env_fig_name)

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
            res, low_t, hi_t = self.opt_double_pulse(t_pts, a_pts, guess_x, f0, phi)
            used_double = True
        else:
            guess_x = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1]])
            res, low_t, hi_t = self.opt_pulse(t_pts, a_pts, guess_x, f0, phi)

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
            ax.annotate(r"$er^2={:.1f}$".format(res.fun), (0.01, 0.84), xycoords='axes fraction')
            fig.savefig(fig_name+".png", dpi=300)
            plt.close(fig)

        if fig_name != '':
            res_name = fig_name+"_res.txt"
            write_reses(res_name, [res])

        return res, None

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
                print("clust={}, j={}\n\tfield error^2={}\n".format(clust, j, err_2*self.n_t_pts))
            #before doing anything else, save a plot of just the time series
            fig_name = ""
            if save_fit_figs:
                plt.plot(self.t_pts, np.real(v_pts))
                plt.savefig("{}/fit_figs/t_series_{}_{}.pdf".format(self.prefix,clust,j))
                plt.clf()
                fig_name = "{}/fit_figs/fit_{}_{}".format(self.prefix,clust,j)
            #now actually try optimizing
            res,_ = self.opt_pulse_full(self.t_pts, np.real(v_pts), err_2, fig_name=fig_name)
            #save the point
            n_evals += 1
            if verbose > 0:
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
