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
        res += p[2]*np.exp( -0.5*((t_pts-p[0])/p[1])**2 )
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
def double_gauss_env(x, t_pts):
    return x[0]*np.exp(-0.5*((t_pts - x[1])/x[2])**2) + x[5]*np.exp(-0.5*((t_pts - x[6])/x[7])**2)

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

#the signal class extracts information about a pulse
class signal:
    def get_freq_fwhm(self):
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
        
    def _param_est(self, freq_range):
        f_min = freq_range[0]
        f_max = freq_range[1]
        #using this frequency range, perform a linear regression to estimate phi and t0
        angles = fix_angle_seq(np.angle(self.vf[f_min:f_max]))
        res = linregress(self.freqs[f_min:f_max], angles)
        t0 = -res.slope/(2*np.pi)
        phi = res.intercept
        return t0, phi
    
    def calc_sumdif(self):
        '''return the sum and difference series expanded about the peak in frequency space.'''
        drange = min(self.f0_ind, len(self.freqs)//2 - self.f0_ind-1)
        self.dif_ser = (self.vf[self.f0_ind:self.f0_ind+drange] - self.vf[self.f0_ind:self.f0_ind-drange:-1])
        self.sum_ser = (self.vf[self.f0_ind:self.f0_ind+drange] + self.vf[self.f0_ind:self.f0_ind-drange:-1])
        end_arr = self.vf[self.f0_ind+drange:]
        self.dif_ser = np.concatenate((self.dif_ser, end_arr))
        self.sum_ser = np.concatenate((self.sum_ser, end_arr))
    
    def apply_lowpass(self, vf0, strength):
        self.vf = vf0*np.sinc((self.freqs-self.freqs[self.f0_ind])*strength)
        self.v_pts = np.roll( np.real(fft.irfft(self.vf)), self.t0_ind )
        self.v_abs = np.abs(self.v_pts)
        
    def guess_center_ind(self):
        '''It is helpful to use signals that are shifted to attain an envelope maxima near t=0. Return a guess of the index where this occurs in self.v_pts'''
        return np.argmax(self.v_abs)
    
    def __init__(self, t_pts, v_pts, lowpass_inc=1.0, scan_length=5, noise_thresh=0.1):
        '''Based on a time series sampled at t_pts with values v_pts, try to find a Gaussian pulse envelope.'''
        self.t_pts = t_pts
        self.dt = t_pts[1]-t_pts[0]
        self.v_abs = np.abs(v_pts)
        #take the fourier transform
        self.freqs = fft.rfftfreq(len(v_pts), d=self.dt)
        self.v_pts = v_pts
        #we shift the signal to have a peak near t=0. The phase in frequency space oscillates at a rate proportional to t0, so this improves numerical stability.
        self.t0_ind = self.guess_center_ind()
        vf0 = fft.rfft(np.roll(v_pts, -self.t0_ind))
        self.vf = np.copy(vf0)
        avg_len = 0.0
        while True:
            self.mags = np.abs(self.vf*np.heaviside(self.freqs, 1))
            #Shift the pulse so that the peak magnitude is at zero frequency. Taking the inverse fourier transform of this shifted magnitude gives pulse amplitude. The frequency at which the magnitude is maximized is the same as the unmodulated sine wave.
            self.f0_ind = np.argmax(self.mags)
            #if there are many peaks that is a sign that this is a noisy signal. We apply a lowpass filter in this case. Default to 2*peak_magnitude to ensure that a lowpass filter is applied in the case where the array is empty.
            start = 4*self.f0_ind-scan_length
            end = 4*self.f0_ind+scan_length
            self.noisiness = np.max(self.mags[start:end]) if start < len(self.mags) else 2*self.mags[self.f0_ind]
            if self.noisiness < noise_thresh*self.mags[self.f0_ind] or avg_len >= MAX_LOWPASS_EVALS*lowpass_inc:
                break
            else:
                avg_len += lowpass_inc
                #shift the sinc function to the peak frequency to avoid erroneously killing off a real oscillation
                self.apply_lowpass(vf0, avg_len)
                if verbose > 0:
                    print("\tNoisy data, applying low pass filter strength={}".format(avg_len))
                                    
        freq_fwhm = self.get_freq_fwhm() 
        #self.envelope = np.roll( np.abs(fft.ifft(np.roll(self.mags, -self.f0_ind))), int(self._t0_corr/self.dt)+self.t0_ind )
        self.f0 = self.freqs[self.f0_ind]
        self._t0_corr, self._phi_corr = self._param_est(freq_fwhm)
        self.t0 = self._t0_corr + t_pts[self.t0_ind]
        self.phi = fix_angle(self._phi_corr + t_pts[self.t0_ind]*2*np.pi*self.f0)
        #now compute the sum and difference vf(omega_0 + delta)+vf(omega_0 - delta) and vf(omega_0 + delta)-vf(omega_0 - delta)
        self.calc_sumdif()
        self.peaks_arr = None
        self.envelope = None
        
    def get_envelope(self):
        '''Get the envelope of the packet under the assumption that said envelope is purely even.'''
        if self.envelope is not None:
            return self.envelope
        #It is much more convenient to work with a full fourier transform that includes negative frequency components. This is equivalent to multiplying Heaviside(f) with the actual transform and shifting so that the peak is at zero frequency.
        full_fft = np.pad(self.mags, (0, self.mags.shape[0]-2))
        self.envelope = np.roll(fft.ifft(np.roll(full_fft, -self.f0_ind)), self.t0_ind)
        return self.envelope
        
    def get_envelope_asym_old(self):
        '''This does the same thing as get_envelope_asym, but it uses a small angle approximation. You probably shouldn't ever use this except as a comparison.'''
        t0_ind = int(self.t0/self.dt)
        prd_ser = -np.conj(self.dif_ser)*self.sum_ser/np.cos(self.t0*self.freqs[0:drange])**2
        rat_ser = self.dif_ser/self.sum_ser
        re_efft = np.real( np.array([np.real(self.sum_ser[0]) if i == 0 else np.sqrt(np.abs(prd_ser[i]/rat_ser[i])) for i in range(len(prd_ser))]) )
        im_efft = np.real(np.sqrt(-prd_ser*rat_ser))
        env_fft = np.pad((re_efft - 1j*im_efft)/self.mags[np.argmax(self.mags)], (0, (t_pts.shape[0]+2)//2 - im_efft.shape[0]))
        full_env = np.roll(fft.irfft(env_fft), self.t0_ind)
        full_phi = self.phi + (t_pts[t0_ind] - t_pts[np.argmax(full_env)])/self.f0
        return full_env, full_phi
        
    def get_envelope_asym(self, lowpass_scale=0.5):
        '''Get the envelope of the packet under the assumption that said envelope is not purely even nor purely odd.'''
        #use some math to figure out the real and imaginary parts of the envelope Fourier transform.
        plt_fs = self.freqs[0:len(self.dif_ser)]
        re_env = np.abs(self.sum_ser)*np.cos(2*np.pi*plt_fs*self._t0_corr) - np.abs(self.dif_ser)*np.sin(2*np.pi*plt_fs*self._t0_corr)
        im_env = np.abs(self.sum_ser)*np.sin(2*np.pi*plt_fs*self._t0_corr) + np.abs(self.dif_ser)*np.cos(2*np.pi*plt_fs*self._t0_corr)
        env_fft = np.pad((re_env - 1j*im_env), (0, self.vf.shape[0]-re_env.shape[0]))
        #apply a lowpass filter first if requested
        if lowpass_scale > 0:
            env_fft = env_fft*np.sinc((self.freqs-self.freqs[self.f0_ind])*lowpass_scale)
        env_test = fft.irfft(env_fft)
        #return np.roll(env_test, self.t0_ind)
        
        #we start off with a reasonable guess for the index that we should shift the envelope by. Then we calculate the response between the envelope and actual data in points near this region. Continue until we reach half of the peak maximum. I found that this guess isn't great, so we calculate the convolution between the envelope and the peaks and shift in time correspondingly.
        guess_ind = int( -(np.angle(env_fft[2]) - np.angle(env_fft[0]))/((self.freqs[2] - self.freqs[0])*2*np.pi*self.dt) ) if len(env_fft) > 2 else 0
        #only the peaks in the pulse are relevant for fitting the envelope
        v_abs = np.abs(self.v_pts)
        pulse_peaks = ssig.argrelmax(v_abs)[0]
        def calc_response(ii):
            ret = 0.0
            for t_ind in pulse_peaks:
                #prevent out of bounds errors
                if np.abs(t_ind-ii) >= env_test.shape[0]:
                    break
                ret += v_abs[t_ind]*np.abs(env_test[t_ind-ii])
            return ret

        best_ind = guess_ind
        best_res = calc_response(guess_ind)
        for j in range(1, self.t_pts.shape[0]-guess_ind-1):
            p_res = calc_response(guess_ind+j)
            m_res = calc_response(guess_ind-j)
            if p_res > best_res:
                best_res = p_res
                best_ind = guess_ind+j
            if m_res > best_res:
                best_res = m_res
                best_ind = guess_ind-j
            if p_res < best_res/2 and m_res < best_res/2:
                break
        return np.roll(env_test, best_ind)

    def get_peaks(self, peak_amp_thresh=0.1, width_steps=5):
        '''Return an array containing peaks in the envelope
        '''
        if self.peaks_arr is not None:
            return self.peaks_arr
        v_abs = np.abs(self.v_pts)
        this_env = self.get_envelope_asym(lowpass_scale=1.0)
        #this_env = self.get_envelope()
        #now that we have an envelope, find its maxima
        pulse_peaks = ssig.argrelmax(v_abs)[0]
        env_peaks = ssig.argrelmax(this_env, order=4)[0]
        max_peak = np.max(this_env.take(env_peaks))
        #return a single default peak if no extrema were found to prevent out of bounds errors
        if len(pulse_peaks) == 0 or len(env_peaks) == 0:
            return np.array([[0, 1, 1, 0]])
        
        #an array of peaks, each specified by time, width, and amplitude
        self.peaks_arr = np.zeros((len(env_peaks),4))
        i = 0
        for v in env_peaks:
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

            if v_abs[p_ind] > peak_amp_thresh*max_peak:
                #make sure that we're inside the array but we take at least one step. If there are any zeros in the envelope array then skip over this one before the program crashes because of NaNs
                stp = max(min(min(width_steps, v), len(self.t_pts)-v-1), 1)
                if stp < 1 or 0 in this_env[v-stp:v+stp]:
                    continue
                #do a quadratic fit to get a guess of the Gaussian width
                samp_ts = self.t_pts[v-stp:v+stp+1] - self.t_pts[v]
                samp_ls = np.log(max_peak) - np.log(np.abs(this_env[v-stp:v+stp+1]))
                reg,cov = opt.curve_fit(lambda x,a,b,c: a*x**2 + b*x+c, samp_ts, samp_ls)
                #only add this to the array if the value isn't nan
                if reg[0] > 0:
                    self.peaks_arr[i, 1] = 1/np.sqrt(2*reg[0])
                    self.peaks_arr[i, 0] = now - reg[1]/(2*reg[0])
                    self.peaks_arr[i, 2] = v_abs[p_ind]*np.sign(this_env[v])
                    #self.peaks_arr[i, 2] = this_env[v]
                    self.peaks_arr[i, 3] = np.sum((peaks_to_pulse(self.t_pts[pulse_peaks], self.peaks_arr[0:i+1]) - v_abs[pulse_peaks])**2)
                    i += 1

        self.peaks_arr = self.peaks_arr[:i]
        return self.peaks_arr

    def opt_phase(self, env, peak_t):
        x0 = np.array([1.0, 2*np.pi*self.f0, self.phi])
        def ff(x):
            return np.sum( (x[0]*env*np.cos(x[1]*(self.t_pts-peak_t)+x[2]) - self.v_pts)**2 )
        res = opt.minimize(ff, x0)
        return res

    def make_raw_plt(self, axs):
        #get the envelope and perform a fitting
        env = self.get_envelope_asym()
        peak_t = np.argmax(env)*self.dt
        res = self.opt_phase(env, peak_t)
        axs.set_title("Frequency space representation of a simulated pulse")
        # Add the magnitude 
        axs[0].set_xlabel('frequency (1/fs)') 
        axs[0].set_ylabel('magnitude (arb. units)', color = 'black') 
        plt1 = axs[0].plot(self.freqs, np.abs(self.vf), color = 'black')
        plt1 = axs[0].plot(self.freqs[:len(self.sum_ser)], np.abs(self.sum_ser)/2, color = 'blue')
        plt1 = axs[0].plot(self.freqs[:len(self.dif_ser)], np.abs(self.dif_ser), color = 'orange')
        #plot the bounds used for fitting
        frange = self.get_freq_fwhm()
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
        axs[1].axvline(peak_t, color='teal', linestyle=':')
        axs[1].plot(self.t_pts, res.x[0]*env*np.cos(res.x[1]*(self.t_pts-peak_t)+res.x[2]), color='orange', label='extracted pulse')

    def make_fit_plt(self, ax):
        #extract gaussian peaks
        s_peaks = self.get_peaks()
        #make plots
        ax.set_title("Comparisons in time space")
        ax.set_xlabel("time (fs)")
        ax.set_ylabel("amplitude (arb. units)")
        ax.set_ylim([-1.0, 1.0])
        ax.plot(self.t_pts, peaks_to_pulse(self.t_pts, s_peaks), color='darkturquoise', label='fitted peaks')
        ax.plot(self.t_pts, self.v_pts, color='black', label='simulated data')
        ax.legend()

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
        '''if amp < 0:
            amp *= -1
            if cep > 0:
                cep -= np.pi
            else:
                cep += np.pi'''
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
    def __init__(self, fname, width, height, pass_alpha=1.0, slice_dir='x', prefix='.', scan_length=5):
        self.prefix = prefix
        self.slice_dir = slice_dir
        self.slice_name = 'z'
        if self.slice_dir == 'z':
            self.slice_name = 'x'
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
        res = opt.minimize(fp, x0)

        #renormalize thingies
        res.x[0] = res.x[0]*env_x[0]
        #res.x[4] = fix_angle(res.x[4])

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

        flipped_res = fix_double_pulse_order(res)
        if fp(flipped_res.x) < res.fun:
            return flipped_res, t_pts[0], t_pts[-1]
        else:
            return res, t_pts[0], t_pts[-1]

    def opt_pulse_full(self, t_pts, a_pts, err_sq, fig_name='', raw_ax=None, final_ax=None):
        '''Put all of the different fitting subroutines together!
        t_pts: the time samples
        a_pts: the electric field at each time sample, this must have the same size as t_pts
        err_sq: an estimate of the square error on each time point
        fig_name: the name of the results file to write fit information to
        raw_ax: if not None, then this is used as a pyplot axes and fit information about the envelope extraction is saved
        final_ax: if not None, then this is used as a pyplot axes and fit information about final time series is saved
        '''
        if a_pts.shape != t_pts.shape:
            raise ValueError("t_pts and a_pts must have the same shape")
        if t_pts.shape[0] == 0:
            raise ValueError("empty time series supplied!")
        #do the signal processing to find the envelope and decompose it into a series of peaks
        #a_pts, peak_arr, f0, phi = self.est_env_new(t_pts, a_pts, fig_name=env_fig_name)
        psig = signal(t_pts, a_pts)
        stp = 5
        peak_arr = psig.get_peaks(width_steps=20)

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
            res, low_t, hi_t = self.opt_double_pulse(t_pts, psig.v_pts, guess_x, 2*np.pi*psig.f0, psig.phi)
            used_double = True
        else:
            guess_x = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1]])
            res, low_t, hi_t = self.opt_pulse(t_pts, psig.v_pts, guess_x, 2*np.pi*psig.f0, psig.phi)

        if verbose > 0:
            print("\tenvelope type={}".format("double" if used_double else "single"))

        if fig_name != '':
            res_name = fig_name+"_res.txt"
            write_reses(res_name, [res])
        #save plots of the frequency domain and envelope extraction
        if raw_ax is not None:
            psig.make_fit_plt(raw_ax)
            psig.save_raw_plt(fig_name+"_env.png")
        #save plots of the final fitted waveform
        if final_ax is not None:
            psig.make_fit_plt(final_ax)
            if used_double:
                x0 = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1], 2*np.pi*psig.f0, psig.phi, peak_arr[1,2], peak_arr[1,0], peak_arr[1,1]])
                final_ax.plot(t_pts, double_gauss_env(x0, t_pts), color='red', linestyle=':')
                final_ax.plot(t_pts, double_gauss_env(res.x, t_pts), color='blue', linestyle=':')
                final_ax.plot(t_pts, double_gauss_series(x0, t_pts), color='red')
                final_ax.plot(t_pts, double_gauss_series(res.x, t_pts), color='blue')
            else:
                final_x0 = np.array([peak_arr[0,2], peak_arr[0,0], peak_arr[0,1], 2*np.pi*psig.f0, psig.phi])
                final_ax.plot(t_pts, gauss_env(x0, t_pts), color='red', linestyle=':')
                final_ax.plot(t_pts, gauss_env(res.x, t_pts), color='blue', linestyle=':')
                final_ax.plot(t_pts, gauss_series(x0, t_pts), color='red')
                final_ax.plot(t_pts, gauss_series(res.x, t_pts), color='blue')
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
            raw_ax = None
            fin_ax = None
            if save_fit_figs:
                plt.plot(self.t_pts, np.real(v_pts))
                plt.savefig("{}/fit_figs/t_series_{}_{}.pdf".format(self.prefix,clust,j))
                plt.clf()
                fig_name = "{}/fit_figs/fit_{}_{}".format(self.prefix,clust,j)
                #setup the plot
                raw_fig, raw_ax = plt.subplots(2)
                fin_fig = plt.figure()
                fin_ax = fin_fig.add_axes([0.1, 0.1, 0.8, 0.8])
            #now actually try optimizing
            res,_ = self.opt_pulse_full(self.t_pts, np.real(v_pts), err_2, fig_name=fig_name, raw_ax=raw_ax, final_ax=fin_ax)
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

    def lookup_fits(self, clust_name, recompute=False, save_fit_figs=False, save_fit_results=True):
        '''Load the fits from the pickle file located at fname or perform the fits if it doesn't exist. Return the results'''
        data_name = '{}/dat_{}'.format(self.prefix, clust_name)
        if recompute or not os.path.exists(data_name):
            #figure out data by phase fitting
            ret = self.read_cluster(clust_name, save_fit_figs=save_fit_figs)
            if save_fit_results:
                with open(data_name, 'wb') as fh:
                    pickle.dump(ret, fh)
            return ret
        else:
            with open(data_name, 'rb') as fh:
                ret = pickle.load(fh)
            return ret
