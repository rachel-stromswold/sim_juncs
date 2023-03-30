import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.stats import linregress
#from scipy.fft import ifft, fft, fftfreq
import scipy.fft as fft
import scipy.signal as ssig
import time
import pickle
import os.path

import matplotlib
import matplotlib.pyplot as plt

import emcee
import corner
from multiprocessing import Pool

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
MAX_PARAM_EVALS=1

DOUB_SWAP_MAT = np.array([[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0]])

#lowpass filter strength applied to asymetric envelopes
DEF_LOW_SCALE = 1.0

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

def fix_angle_seq(angles):
    '''Given a sequence of angles from -pi to pi look for -pi -> pi rollovers and perform adjustments such that each angle is the closest equivalent angle to the previous.'''
    shift = 0
    n = len(angles)
    for i in range(1,n):
        dist_raw = abs(angles[i] - angles[i-1])
        if abs(angles[i]+2*np.pi - angles[i-1]) < dist_raw:
            for j in range(i,n):
                angles[j] += 2*np.pi
        elif abs(angles[i]-2*np.pi - angles[i-1]) < dist_raw:
            for j in range(i,n):
                angles[j] -= 2*np.pi
    return angles

def peaks_to_pulse(t_pts, peaks):
    res = np.zeros(len(t_pts))
    for p in peaks:
        res += p[2]*np.exp( -0.5*((t_pts-p[0])/p[1])**2 )
    return res

#all of this is for debuging, first pick a cluster then get the double guess
def double_gauss_env(x, t_pts):
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

class dummy_result:
    def __init__(self, amp, omega, phi):
        self.x = np.array([amp, 1, 1, omega, phi])
        #todo: give these useful values
        self.fun = 0.0
        self.hess_inv = np.zeros((5,5))
        self.jac = np.zeros(5)

class signal:
    def _guess_t0_ind(self):
        #guess t0 by looking at the magnitude peaks
        mag_peaks,_ = ssig.find_peaks(self.v_pts**2)
        avg_t = 0.0
        tot_mag = 0.0
        for i in mag_peaks:
            avg_t += i*self.v_abs[i]
            tot_mag += self.v_abs[i]
        if tot_mag <= 0.0:
            return 0
        return int( avg_t / tot_mag )

    def _get_freq_fwhm(self):
        #find the region in frequency that has magnitude greater than max_freq/2
        cut_amp = self.mags[self.f0_ind]/2
        f_min = self.f0_ind
        f_max = self.f0_ind
        j_bound = min(self.f0_ind, len(self.mags)-self.f0_ind-1)
        if verbose > 2 and len(self.mags)-self.f0_ind-1 < self.f0_ind:
            print("\tanomalous high frequency! peak f={}".format(self.freqs[self.f0_ind]))
        for j in range(j_bound):
            if self.mags[self.f0_ind - j] > cut_amp:
                f_min = self.f0_ind - j - 1
            if self.mags[self.f0_ind + j] > cut_amp:
                f_max = self.f0_ind + j + 1
        return [f_min, f_max]

    def _param_est(self, n_evals=0):
        freq_range = self._get_freq_fwhm()
        f_min = freq_range[0]
        f_max = freq_range[1]
        #using this frequency range, perform a linear regression to estimate phi and t0
        angles = fix_angle_seq(np.angle(self.vf[f_min:f_max]))
        res = linregress(self.freqs[f_min:f_max], angles)
        self._t0_corr = -res.slope/(2*np.pi)
        self._phi_corr = fix_angle(-res.intercept)
        #if there were a lot of oscillations within the fwhm, that is a sign that the signal may need to be shifted
        if np.abs(self._t0_corr) > 2/(f_max-f_min) and n_evals < MAX_PARAM_EVALS:
            if verbose > 0:
                print("\trecentering pulse!")
            self.t0_ind += int(self._t0_corr/self.dt)
            tmp_vf = 2*fft.rfft( np.roll(self.v_pts, -self.t0_ind) )
            self._apply_lowpass(tmp_vf, self._last_cent_f, self._low_stren)
            self.t0_ind -= int(self._t0_corr/self.dt)
            self._param_est(n_evals=n_evals+1)

    def _calc_sumdif(self, f0=-1):
        '''return the sum and difference series expanded about the peak in frequency space.'''
        if f0 < 0:
            #freq_range = self._get_freq_fwhm()
            #mid = (freq_range[0] + freq_range[1])//2
            #mid = self.f0_ind
            mid = self.f0_ind
        else:
            mid = min( len(self.freqs)-1, int(f0/(self.freqs[1]-self.freqs[0])) )
        drange = min(mid, len(self.freqs)//2 - mid-1)
        self.dif_ser = (self.vf[mid:mid+drange] - self.vf[mid:mid-drange:-1])
        self.sum_ser = (self.vf[mid:mid+drange] + self.vf[mid:mid-drange:-1])
        end_arr = self.vf[mid+drange:]
        self.dif_ser = np.concatenate((self.dif_ser, end_arr))
        self.sum_ser = np.concatenate((self.sum_ser, end_arr))

    def _apply_lowpass(self, vf0, cent_freq, strength):
        self.vf = vf0*np.sinc((self.freqs-cent_freq)*strength)
        self.mags = np.abs(self.vf)
        self.v_pts = np.roll( np.real(fft.irfft(self.vf)), self.t0_ind )
        self.v_pts = np.pad( self.v_pts, (0,len(self.t_pts)-len(self.v_pts)) )
        self.v_abs = np.abs(self.v_pts)
        self._last_cent_f = cent_freq

    def is_noisy(self, cutoff, scan_length=5, noise_thresh=0.2):
        '''Checks whether a signal is noisy by testing how strong high frequency components (arbitrarily chosen to be twice the central frequency) are compared to the central frequency.
        cutoff: The relative strength of the high frequency component to central frequency at which a cutoff is performed
        '''
        mag_peaks,_ = ssig.find_peaks(self.mags, prominence=noise_thresh*self.mags[self.f0_ind])
        if len(mag_peaks) < 2:
            return False
        for i in mag_peaks:
            if i > self.f0_ind*3:
                return True
        return False

    def __init__(self, t_pts, v_pts, scan_length=5, lowpass_inc=1.0, noise_thresh=0.1, max_f0=-1, lowpass_center=-1.0):
        '''Based on a time series sampled at t_pts with values v_pts, try to find a Gaussian pulse envelope.
        t_pts: a numpy array with the times of each sampled point. This must have the same dimensions as v_pts
        v_pts: a numpy array with the field at each sampled time
        lowpass_inc: in the case of a bad fit, a lowpass filter with increasingly high strength is applied. This specifies the increment at which the strength of the lowpass filter increases.
        scan_length: This is used when checking for noise
        noise_thresh: This is the ratio of peak frequency magnitude to high frequency components at which a signal is decided to be noisy.
        max_f0: This specifies the maximum allowable frequency at which the central frequency may occur. If the central frequency is found to be higher than this value, then a lowpass filter is applied.
        lowpass_center: if >0, then lowpass filters are centered around this frequency so as to avoid attenuating genuine signal
        '''
        self.t_pts = t_pts
        self.v_pts = v_pts
        #inferred values
        self.dt = t_pts[1]-t_pts[0]
        self.v_abs = np.abs(v_pts)
        #self.t0_ind = np.argmax(self.v_abs)
        self.t0_ind = self._guess_t0_ind()
        #take the fourier transform       
        self.freqs = fft.rfftfreq(len(v_pts), d=self.dt)
        max_f0_ind = len(self.freqs)//2
        if max_f0 > 0:
            max_f0_ind = int(max_f0 / (self.freqs[1]-self.freqs[0]))
        #Shift the pulse so that the peak magnitude is at zero frequency. Taking the inverse fourier transform of this shifted magnitude gives pulse amplitude. The frequency at which the magnitude is maximized is the same as the unmodulated sine wave.
        vf0 = 2*fft.rfft(np.roll(v_pts, -self.t0_ind))
        self.vf = np.copy(vf0) 
        self.mags = np.abs(self.vf)
        self._low_stren = 0.0
        self._last_cent_f = 0.0
        while True:
            self.phi = np.angle(self.vf[0])
            self.f0_ind = np.argmax(self.mags)
            #apply lowpass filters to noisy signals
            if self._low_stren >= MAX_LOWPASS_EVALS*lowpass_inc:
                break
            #there are two seperate cases we check for, one is a suspiciously high central frequency, and the other is strong signal well above the (otherwise normal) central frequency
            cent_freq = -1.0
            if self.f0_ind > max_f0_ind:
                cent_freq = 0.0
            elif self.is_noisy(noise_thresh, scan_length=scan_length):
                cent_freq = self.freqs[self.f0_ind]/2
            if cent_freq >= 0.0:
                self._low_stren += lowpass_inc
                if lowpass_center >= 0.0:
                    self._apply_lowpass(vf0, lowpass_center, self._low_stren)
                else:
                    self._apply_lowpass(vf0, cent_freq, self._low_stren)
                if verbose > 0:
                    print("\tNoisy data, applying low pass filter strength={}, f0={}".format(self._low_stren, cent_freq))
            else:
                break

        #perform post-processing
        self._param_est()
        if verbose > 1:
            print("\tt0_corr: {}\n\tphi_corr: {}".format(self._t0_corr, self._phi_corr))
        self.envelope = None
        self.envelope_asym = None
        self.asym_lowpass = DEF_LOW_SCALE
        self.asym_f0 = None
        self.peaks_arr = None
        self.f0 = self.freqs[self.f0_ind]
        self.t0_ind = np.argmax(self.v_abs)

    def get_envelope(self):
        if self.envelope is not None:
            return self.envelope
        #take the inverse fourier transform and shift to the peak of the pulse
        full_fft = np.pad(self.mags, (0, self.mags.shape[0]-2))
        self.envelope = np.abs( np.roll(fft.ifft(np.roll(full_fft, -self.f0_ind)), self.t0_ind) )
        #for odd numbers of time points, the inverse fourier transform will have one fewer points
        if len(self.t_pts) > len(self.envelope):
            self.envelope = np.pad( self.envelope, (0,len(self.t_pts)-len(self.envelope)) )
        return self.envelope

    def get_envelope_asym(self, t0=-1, f0=-1):
        '''Get the envelope of the packet under the assumption that said envelope is not purely even nor purely odd.'''
        if self.envelope_asym is not None and self.asym_f0 is not None and self.asym_f0 == f0:
            return self.envelope_asym
        if self.asym_f0 is None or f0 != self.asym_f0:
            self._calc_sumdif(f0=f0)
        if t0 < 0:
            t0 = self._t0_corr
        t_shift = self.t_pts[self.t0_ind]-self._t0_corr
        #use some math to figure out the real and imaginary parts of the envelope Fourier transform.
        plt_fs = self.freqs[0:len(self.dif_ser)]
        c_ser = np.cos(2*np.pi*plt_fs*t0)
        s_ser = np.sin(2*np.pi*plt_fs*t0)        
        rec_env = np.real(self.sum_ser)*c_ser - np.imag(self.dif_ser)*s_ser
        imc_env = np.real(self.sum_ser)*s_ser + np.imag(self.dif_ser)*c_ser
        res_env = np.imag(self.sum_ser)*c_ser + np.real(self.dif_ser)*s_ser
        ims_env = np.imag(self.sum_ser)*s_ser - np.real(self.dif_ser)*c_ser
        #compute the magnitude and the angle of the fourier transform of the envelope. if there are any nans just set them to zero
        mag_env = np.sqrt( rec_env**2 + res_env**2 + imc_env**2 + ims_env**2 )
        ang_env = np.arctan(imc_env/rec_env)
        ang_env[np.isnan(ang_env)] = 0
        shift_fact = np.exp( 1j*(ang_env-2*np.pi*plt_fs*t_shift) )
        env_fft = np.pad(mag_env*shift_fact, (0, self.vf.shape[0]-mag_env.shape[0]))
        return fft.irfft(env_fft)

    def sample_params(self, verr, countour_fname=''):
        '''Use MCMC sampling near the estimated parameters t0 shift, t0 correction, central frequency, and phase'''
        theta_0 = np.array([self._t0_corr, self.f0, self._phi_corr])
        t0_span = 0.5/self.f0
        print(t0_span)
        t0_range = [self._t0_corr-t0_span, self._t0_corr+t0_span]
        t_peaks = ssig.argrelmax(self.v_pts)[0]
        def log_prior(t_0, omega_0, phi):
            if t0_range[0] < t_0 < t0_range[1] and 0.0 < omega_0 < 2*np.pi*self.f0 and -np.pi < phi < np.pi:
                return 0.0
            return -np.inf
        def log_prob(theta, ts, vs, er):
            env = self.get_envelope_asym(t0=theta[0], f0=theta[1])
            #get parameters and the prior distribution
            t_0 = theta[0]
            omega_0 = 2*np.pi*theta[1]
            phi = theta[2]
            lp = log_prior(t_0, omega_0, phi)
            if not np.isfinite(lp):
                return -np.inf
            #actually compute the log likelihood
            return lp - (0.5/er**2)*np.sum( (env*np.cos(omega_0*(ts - t_0) + phi) - vs)**2 )
        pos = np.random.normal(theta_0, [t0_span/2, self.f0/4, np.pi/2], size=(32, 3))
        nwalkers, ndim = pos.shape
        with Pool() as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(self.t_pts, self.v_pts, verr))
            sampler.run_mcmc(pos, 2000, progress=True)
        flat_samples = sampler.get_chain(discard=200, thin=20, flat=True)
        if countour_fname != ''
            fig = corner.corner(flat_samples, labels=["t0", "f0", "phi"], truths=theta_0)
            fig.savefig(countour_fname)
        #find mean and one standard deviation from mcmc samples
        theta_f = np.zeros((ndim, 3))
        for i in range(ndim):
            percs = np.percentile(flat_samples[:,i], [16, 50, 84])
            q = np.diff(percs)
            theta_f[i, 0] = percs[1]
            theta_f[i, 1] = q[0]
            theta_f[i, 2] = q[1]
        return theta_f

    def get_skewness(self):
        #get the asymmetric envelope and normalize
        env = np.abs( self.get_envelope_asym() )
        norm = np.trapz(env)
        env = env / norm
        #Calculate mu and sigma. Use these to get the skewness
        mu = np.trapz(env*self.t_pts)
        sig = np.sqrt( np.trapz(env*(self.t_pts-mu)**2) )
        return np.trapz(env*(self.t_pts - mu)**3)/(sig**3)

    def get_peaks_new(self, peak_amp_thresh=0.1):
        if self.peaks_arr is not None:
            return self.peaks_arr
        this_env = self.get_envelope_asym()
        max_peak = np.max(this_env)
        pulse_peaks = ssig.argrelmax(self.v_abs)[0]
        env_peaks,env_props = ssig.find_peaks(this_env, prominence=max_peak*0.1)
        #return a single default peak if no extrema were found to prevent out of bounds errors
        if len(pulse_peaks) == 0 or len(env_peaks) == 0:
            return np.array([[0, 1, 1, 0]])
        #an array of peaks, each specified by time, width, and amplitude
        self.peaks_arr = np.zeros((len(env_peaks),4))
        n = 0
        for i, v in enumerate(env_peaks):
            now = self.t_pts[v]
            #Look only at points near this peak. Make sure that we're inside the array but we take at least one step. If there are any zeros in the envelope array then skip over this one before the program crashes because of NaNs
            stp = min(v-env_props['left_bases'][i], env_props['right_bases'][i]-v)
            if stp < 1 or 0 in this_env[v-stp:v+stp]:
                continue
            #do a quadratic fit to get a guess of the Gaussian width
            samp_ts = self.t_pts[v-stp:v+stp+1] - self.t_pts[v]
            samp_ls = np.log(max_peak) - np.log(np.abs(this_env[v-stp:v+stp+1]))
            reg,cov = opt.curve_fit(lambda x,a,b,c: a*x**2 + b*x + c, samp_ts, samp_ls)
            self.peaks_arr[n, 1] = self.dt*stp/np.sqrt(2)
            self.peaks_arr[n, 0] = now
            self.peaks_arr[n, 2] = this_env[v]
            n += 1
        #now sort peaks by their heights descending and calculate the square errors for including each
        self.peaks_arr.resize((n,4))
        self.peaks_arr = self.peaks_arr[np.argsort(-self.peaks_arr[:,2])]
        #self.peaks_arr = self.peaks_arr[np.argsort(self.peaks_arr[:,3])]
        for i in range(n):
            self.peaks_arr[i,3] = np.sum((peaks_to_pulse(self.t_pts[pulse_peaks], self.peaks_arr[0:i+1]) - self.v_abs[pulse_peaks])**2)

        return self.peaks_arr

    def _get_region(self, ns, i, width_steps):
        #make sure that we're inside the array but we take at least one step
        stp = max(min(min(width_steps, i), len(self.t_pts)-i-1), 1) 
        #do a quadratic fit to get a guess of the Gaussian width
        samp_ts = self.t_pts[i-stp:i+stp] - self.t_pts[i]
        samp_ls = ns[i-stp:i+stp]
        return samp_ts, samp_ls

    def get_peaks(self, peak_amp_thresh=0.1, width_steps=5):
        '''Helper function for est_env_new
        t_pts: the sampled points in time, this should have the same dimension as vs and ns
        vs: the absolute value of the actual pulse. This WILL break if there are any negative values.
        ns: the envelope of the pulse obtained through some signal processing
        '''
        if self.peaks_arr is not None:
            return self.peaks_arr
        ns = self.get_envelope()
        #now that we have an envolope, find its maxima
        pulse_peaks = ssig.argrelmax(self.v_abs)[0]
        max_peak = np.max(ns)
        env_peaks,env_props = ssig.find_peaks(ns, prominence=max_peak*0.05)
        #return a single default peak if no extrema were found to prevent out of bounds errors
        if len(pulse_peaks) == 0 or len(env_peaks) == 0:
            return np.array([[0, 1, 1]])
        
        #an array of peaks. The columns are time, width, amplitude, and the residual error. The error is cumulative i.e. The pulse incorporates all peaks before the current one and computes the square error.
        self.peaks_arr = np.zeros((len(env_peaks),4))
        n = 0
        for i, v in enumerate(env_peaks):
            if v >= len(self.t_pts):
                break
            now = self.t_pts[v]
            #find the nearest maxima in the actual pulse. The envelopes aren't great at getting amplitudes, but they are useful for widths and centers.
            p_ind = 0
            closest_sep = self.t_pts[-1] - self.t_pts[0]
            for p in pulse_peaks:
                sep = abs(self.t_pts[p] - now)
                if sep < closest_sep:
                    p_ind = p
                    closest_sep = sep
                else:
                    break
            if self.v_abs[p_ind] > peak_amp_thresh*max_peak:
                #make sure that we're inside the array but we take at least one step
                stp = min(v-env_props['left_bases'][i], env_props['right_bases'][i]-v)//4
                if stp <= 0:
                    continue
                samp_ts, samp_ls = self._get_region(ns, v, stp)
                reg,cov = opt.curve_fit(lambda x,a,b: a*x**2 + b, samp_ts, np.log(max_peak)-np.log(samp_ls))
                this_w = 1/np.sqrt(reg[0])
                #only add this to the array if the value isn't nan
                if not np.isnan(this_w):
                    self.peaks_arr[n, 0] = now
                    self.peaks_arr[n, 1] = this_w
                    self.peaks_arr[n, 2] = self.v_abs[p_ind]
                    #self.peaks_arr[n, 3] = np.sum((peaks_to_pulse(self.t_pts[pulse_peaks], self.peaks_arr[n:n+1]) - self.v_abs[pulse_peaks])**2)
                    n += 1

        #sort peaks by intensity
        self.peaks_arr.resize((n,4))
        self.peaks_arr = self.peaks_arr[(-self.peaks_arr)[:, 2].argsort()]
        #self.peaks_arr = self.peaks_arr[(self.peaks_arr)[:, 3].argsort()]
        for i in range(len(self.peaks_arr)):
            self.peaks_arr[i,3] = np.sum((peaks_to_pulse(self.t_pts[pulse_peaks], self.peaks_arr[0:i+1]) - self.v_abs[pulse_peaks])**2)
        return self.peaks_arr

    def opt_phase(self, env, peak_t):
        n_pts = min(min(env.shape[0], self.t_pts.shape[0]), self.v_pts.shape[0])
        nenv = env[:n_pts]
        ts = self.t_pts[:n_pts]
        vs = self.v_pts[:n_pts]
        x0 = np.array([1.0, 2*np.pi*self.f0, self.phi])
        def ff(x):
            return np.sum( (x[0]*nenv*np.cos(x[1]*(ts-peak_t)+x[2]) - vs)**2 )
        res = opt.minimize(ff, x0)
        return res

    def extract_phi_asym(self):
        '''returns: tuple with amplitude and phase'''
        env = self.get_envelope_asym()
        peak_i = np.argmax(env)
        peak_t = peak_i*self.dt
        res = self.opt_phase(env, peak_t)
        return env[peak_i]*res.x[0], res.x[2]

    def make_raw_plt(self, axs):
        #get the envelope and perform a fitting
        env = self.get_envelope_asym()
        peak_t = self.t_pts[self.t0_ind] - self._t0_corr
        axs[0].set_title("Frequency space representation of a simulated pulse")
        # Add the magnitude 
        axs[0].set_xlabel('frequency (1/fs)') 
        axs[0].set_ylabel('magnitude (arb. units)', color = 'black') 
        plt1 = axs[0].plot(self.freqs, np.abs(self.vf), color = 'black')
        plt1 = axs[0].plot(self.freqs[:len(self.sum_ser)], np.abs(self.sum_ser)/2, color = 'blue')
        plt1 = axs[0].plot(self.freqs[:len(self.dif_ser)], np.abs(self.dif_ser), color = 'orange')
        #plot the bounds used for fitting
        frange = self._get_freq_fwhm()
        axs[0].axvline(self.freqs[frange[0]], color='gray', linestyle=':')
        axs[0].axvline(self.freqs[frange[1]], color='gray', linestyle=':')
        axs[0].tick_params(axis ='y', labelcolor = 'black')
        # Adding Twin Axes
        ax2 = axs[0].twinx()
        ax2.set_ylabel('phase', color = 'green')
        ax2.set_ylim([-np.pi, np.pi])
        ax2.tick_params(axis ='y', labelcolor = 'green')
        ax2.plot(self.freqs, np.angle(self.vf), color='green')
        #save time domain
        axs[1].plot(self.t_pts, self.v_pts, color='black', label='simulated data')
        axs[1].plot(self.t_pts, env, color='teal', label='envelope')
        axs[1].plot(self.t_pts, self.get_envelope(), color='red', label='envelope')
        axs[1].axvline(peak_t, color='teal', linestyle=':')
        axs[1].plot(self.t_pts, env*np.cos(2*np.pi*self.f0*(self.t_pts-peak_t)+self._phi_corr), color='orange', label='extracted pulse')

    def make_fit_plt(self, ax):
        #extract gaussian peaks
        s_peaks = self.get_peaks()
        #make plots
        ax.set_title("Comparisons in time space")
        ax.set_xlabel("time (fs)")
        ax.set_ylabel("amplitude (arb. units)")
        ax.plot(self.t_pts, peaks_to_pulse(self.t_pts, s_peaks), color='darkturquoise', label='fitted peaks')
        ax.errorbar(s_peaks[:,0], s_peaks[:,2], xerr=s_peaks[:,1], color='darkturquoise', linestyle=None, marker='s')
        ax.plot(self.t_pts, self.v_pts, color='black', label='simulated data')
        ax.legend()

N_RES_PARAMS = 7
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
    def get_amp_ref(self):
        return self.res_arr[10]
    def get_amp_ref_err(self):
        return self.res_arr[11]
    def get_skew(self):
        return self.res_arr[12]
    def get_skew_err(self):
        return self.res_arr[13]
    def get_err_sq(self):
        return self.res_arr[2*N_RES_PARAMS]

    def set_point(self, jj, res, skew, err_2):
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
        self.res_arr[12,jj] = skew
        if len(res.x) > 5:
            self.res_arr[10,jj] = np.abs(res.x[5])
            self.res_arr[11,jj] = np.sqrt(res.hess_inv[5][5]/err_2)/2
        self.res_arr[2*N_RES_PARAMS,jj] = res.fun

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
    def __init__(self, fname, pass_alpha=1.0, slice_dir='x', prefix='.', keep_n=2.5, scan_length=5, skip_fits=False):
        self.prefix = prefix
        self.slice_dir = slice_dir
        self.slice_name = 'z'
        if self.slice_dir == 'z':
            self.slice_name = 'x'
        self.keep_n = keep_n
        self.scan_length = scan_length
        self.skip_fits = skip_fits
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
        #figure out the central frequency
        self.in_freq = self.f['info']['sources']['wavelen'][0] / .299792458
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
        init_fp = fp(x0)
        x0[4] += np.pi
        if fp(x0) > init_fp:
            x0[4] = est_phi
        else:
            est_phi += np.pi
            init_fp = fp(x0)
        x0[4] = fix_angle(x0[4])

        #now we actually perform the minimization
        res = opt.minimize(fp, x0, jac=jac_fp)
        #sometimes the signal optimizes to zero, try flipping around if this happens
        if res.fun > 10:
            x0[4] += np.pi
            x0[4] = fix_angle(x0[4])
            res2 = opt.minimize(fp, x0, jac=jac_fp)
            if res2.fun < res.fun:
                res = res2

        #renormalize thingies
        res.x[0] = res.x[0]*env_x[0]
        res.x[4] = fix_angle(res.x[4])
        return res, t_pts[0], t_pts[-1], init_fp
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
        init_fp = fp(x0)
        x0[4] += np.pi
        if fp(x0) > init_fp:
            x0[4] = est_phi
            init_fp = fp(x0)
        else:
            est_phi += np.pi
        x0[4] = fix_angle(x0[4])

        #now we actually perform the minimization
        res = opt.minimize(fp, x0)
        #sometimes the signal optimizes to zero, try flipping around if this happens
        if res.fun > 10:
            x0[4] += np.pi
            x0[4] = fix_angle(x0[4])
            res2 = opt.minimize(fp, x0, jac=jac_fp)
            if res2.fun < res.fun:
                res = res2

        #renormalize thingies
        res.x[0] = np.abs(res.x[0])*env_x[0]
        res.x[5] = np.abs(res.x[5])*env_x[0]
        return res, t_pts[0], t_pts[-1], init_fp

    def opt_pulse_full(self, t_pts, a_pts, err_sq, fig_name='', raw_ax=None, final_ax=None):
        if a_pts.shape != t_pts.shape:
            raise ValueError("t_pts and a_pts must have the same shape")
        if t_pts.shape[0] == 0:
            raise ValueError("empty time series supplied!")

        #do the signal processing to find the envelope and decompose it into a series of peaks
        psig = signal(t_pts, a_pts, lowpass_inc=2.0)
        w0 = 2*np.pi*psig.f0
        phi = psig.phi
        if self.skip_fits:
            print("skipping!")
            amp = np.max(np.abs(a_pts))
            res = dummy_result(amp, w0, psig._phi_corr)
            return res, psig.get_skewness()

        peak_arr = psig.get_peaks()

        if verbose > 2:
            print("\tinitial peaks: {}\n\tinitial omega_0: {}\n\tinitial phi: {}".format(peak_arr[0:2], w0, phi))

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
            res, low_t, hi_t, init_fun = self.opt_double_pulse(t_pts, a_pts, guess_x, w0, phi)
            used_double = True
        else:
            guess_x = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1]])
            res, low_t, hi_t, init_fun = self.opt_pulse(t_pts, a_pts, guess_x, w0, phi)

        if verbose > 0:
            print("\tenvelope type={}".format("double" if used_double else "single"))

        #plot figures if requested
        if fig_name != '':
            res_name = fig_name+"_res.txt"
            write_reses(res_name, [res, {"fun": init_fun, "x": guess_x}])
        #save plots of the frequency domain and envelope extraction
        if raw_ax is not None:
            psig.make_raw_plt(raw_ax)
        #save plots of the final fitted waveform
        if final_ax is not None:
            psig.make_fit_plt(final_ax)
            final_ax.set_ylim([-0.2, 0.2])
            if used_double:
                x0 = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1], 2*np.pi*psig.f0, psig.phi, peak_arr[1,2], peak_arr[1,0], peak_arr[1,1]])
                final_ax.plot(t_pts, double_gauss_env(x0, t_pts), color='red', linestyle=':')
                final_ax.plot(t_pts, double_gauss_env(res.x, t_pts), color='blue', linestyle=':')
                final_ax.plot(t_pts, double_gauss_series(x0, t_pts), color='red')
                final_ax.plot(t_pts, double_gauss_series(res.x, t_pts), color='blue')
            else:
                x0 = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1], 2*np.pi*psig.f0, psig.phi])
                final_ax.plot(t_pts, gauss_env(x0, t_pts), color='red', linestyle=':')
                final_ax.plot(t_pts, gauss_env(res.x, t_pts), color='blue', linestyle=':')
                final_ax.axvline(x=low_t, color='gray', linestyle=':')
                final_ax.axvline(x=hi_t, color='gray', linestyle=':')
                final_ax.plot(t_pts, gauss_series(x0, t_pts), color='red')
                final_ax.plot(t_pts, gauss_series(res.x, t_pts), color='blue')

        return res, psig.get_skewness()

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
            fin_ax = None
            raw_ax = None
            if save_fit_figs:
                #setup the plots
                fig_name = "{}/fit_figs/fit_{}_{}".format(self.prefix,clust,j)
                raw_fig, raw_ax = plt.subplots(2)
                fin_fig = plt.figure()
                fin_ax = fin_fig.add_axes([0.1, 0.1, 0.8, 0.8])
            #now actually try optimizing
            res,skew = self.opt_pulse_full(self.t_pts, np.real(v_pts), err_2, fig_name=fig_name, raw_ax=raw_ax, final_ax=fin_ax)
            #save the point
            n_evals += 1
            if verbose > 0:
                print("\tsquare errors = {}\n\tx={}\n\tdiag(H^-1)={}".format(res.fun, res.x, np.diagonal(res.hess_inv)))
            good_js.append(j)
            ret.set_point(j, res, skew, err_2)
            if verbose > 0:
                print("\tbad fit!".format(clust,j))
            #save fit figs to disk if requested
            if save_fit_figs:
                raw_fig.savefig("{}/fit_figs/fit_cluster_{}_{}_raw.svg".format(self.prefix, clust, j))
                fin_fig.savefig("{}/fit_figs/fit_cluster_{}_{}.svg".format(self.prefix, clust, j))
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
