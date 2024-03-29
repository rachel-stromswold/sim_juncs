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

verbose = 5
#to avoid divisions by zero shift both the numerator and denominator by a small amount
DIV_S = 4000
#when fitting to a generalized Gaussian, expand the polynomial in the exponent up to this order.
DEF_DEG = 4

#lowpass filter strength applied to asymetric envelopes
DEF_LOW_SCALE = 1.0
#this is the parameter alpha used to specify the range on which linear regression is performed
CUT_ALPHA = 0.1
#we need to stop recursion for parameter estimation and lowpass evals. Note that a value of 1 indicates only one iteration is performed (no recursion)
MAX_PARAM_EVALS=1
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

def rtc_fft(vals, f0i):
    '''Convert a numpy array in the scipy.rfftfreq() format to one in the scipy.fftfreq() format rolled by -f0i'''
    return np.roll(np.pad(vals, (0,vals.shape[0]-2)), -f0i)

class extract_info:
    def __init__(self, amp, omega, phi, t0=1, w=1):
        self.x = np.array([amp, t0, w, omega, phi])
        #todo: give these useful values
        self.fun = 0.0
        self.hess_inv = np.zeros((5,5))
        self.jac = np.zeros(5)
    def keys(self):
        return ["x"]
    def get(self, reskey):
        if reskey == "x":
            return self.x
        return None

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

    def _get_freq_fwhm(self, f0i=-1, threshold=CUT_ALPHA):
        if f0i < 0:
            f0i = self.f0_ind
        #find the region in frequency that has magnitude greater than max_freq/2
        cut_amp = self.mags[f0i]*threshold
        f_min = f0i
        f_max = f0i
        j_bound = min(f0i, len(self.mags)-f0i-1)
        if verbose > 2 and len(self.mags)-f0i-1 < f0i:
            print("\tanomalous high frequency! peak f={}".format(self.freqs[f0i]))
        for j in range(j_bound):
            if self.mags[f0i - j] > cut_amp:
                f_min = f0i - j - 1
            if self.mags[f0i + j] > cut_amp:
                f_max = f0i + j + 1
        return f_min, f_max

    def _param_est(self, n_evals=0, threshold=CUT_ALPHA, axs=None):
        f_min,f_max = self._get_freq_fwhm(threshold=threshold)
        angles = fix_angle_seq(np.angle(self.vf[f_min:f_max]))
        #using this frequency range, perform a linear or cubic fit to estimate phi and t0
        res = linregress(self.freqs[f_min:f_max], angles)
        self._t0_corr = -res.slope/(2*np.pi)
        self._phi_corr = fix_angle(res.intercept)
        print(self._phi_corr)
        #self._phi_corr = np.sum(angles)/len(angles)
        if verbose > 2:
            print("\tparameter estimation: arg(a(t)) ~ {}f + {}".format(-self._t0_corr, self._phi_corr))
        if axs is not None and n_evals == MAX_PARAM_EVALS-1:
            dat_series = [fix_angle(ang)/np.pi for ang in angles]
            fit_series_0 = (res.slope*self.freqs[f_min:f_max] + res.intercept)/np.pi
            fit_series_1 = (res.slope*self.freqs[f_min:f_max] + res.intercept)/np.pi - 2
            axs.scatter(2*np.pi*self.freqs[f_min:f_max], dat_series, color='green', label="simulation")
            axs.plot(2*np.pi*self.freqs[f_min:f_max], fit_series_0, color='green', linestyle='--', label="linear fit")
            axs.plot(2*np.pi*self.freqs[f_min:f_max], fit_series_1, color='green', linestyle='--', label="linear fit")
            axs.annotate('$t_0 = ${:.2f} fs, $\\varphi = ${:.2f}\n$R^2 = ${:.2f}'.format(self._t0_corr, self._phi_corr, res.rvalue**2), xy=(0,0), xytext=(0.2, 0.80), xycoords='figure fraction')
            axs.set_ylim(-1,1)
        #if there were a lot of oscillations within the fwhm, that is a sign that the signal may need to be shifted
        if np.abs(self._t0_corr) > 2/(f_max-f_min) and n_evals < MAX_PARAM_EVALS:
            if verbose > 0:
                print("\trecentering pulse!")
            self.t0_ind += int(self._t0_corr/self.dt)
            tmp_vf = 2*fft.rfft( np.roll(self.v_pts, -self.t0_ind) )
            self._apply_lowpass(tmp_vf, self._last_cent_f, self._low_stren)
            self._param_est(n_evals=n_evals+1, threshold=threshold, axs=axs)
            self.t0 = self.t_pts[self.t0_ind]

    def _fitted_v(self, x, fs, deg):
        #get the generalized Gaussian parameterized by the vector x up to degree deg
        rexp = -x[1]*np.ones(len(fs))
        for j in range(deg//2):
            rexp = rexp - x[j+2]*(fs-x[0])**(2*j+2)
        return np.exp(rexp)
    def _fitted_e(self, x, fs, deg):
        #get the generalized Gaussian multiplied by omega (to get the a(omega as opposed to b(omega)) parameterized by the vector x up to degree deg
        return 2*np.pi*fs*self._fitted_v(x, fs, deg)
    def _opt_gauss_envelope(self, deg):
        '''Helper function for _fit_sym_amps and _est_f0 which optimizes the frequency space magnitudes to the form 2*pi*f*exp(sum_k a_2k(f-f0)**2k)
        deg: the degree of the polynomial in the exponent to use when fitting
        returns: a scipy optimize result and a frequency scaling factor that was used to perform the fit'''
        #setup an initial guess using the peak in frequency space
        max_ind = self.f0_ind//4
        df = self.freqs[1] - self.freqs[0]
        #things blow up really fast if we don't rescale the terms in the exponent to order unity
        fscale = 1/self.freqs[-1]
        fs = self.freqs*fscale
        x0 = np.zeros(deg//2 + 2)
        x0[0] = fs[max_ind]
        x0[1] = -np.log(self.mags[max_ind])
        f_min,f_max = self._get_freq_fwhm(f0i=max_ind, threshold=0.5)
        x0[2] = 0.5/(fs[f_max]-fs[f_min])**2
        #gradient of cost function
        def jac_fp(x, fe):
            #2*residual*omega*e^(-sum ...) appears in all the Jacobian terms
            rexp = -x[1]*np.ones(len(fs))
            for j in range(deg//2):
                rexp = rexp - x[j+2]*(fs-x[0])**(2*j+2)
            #res_w = 4*np.pi*self.freqs*(2*np.pi*self.freqs*fe(x, fs, deg)-self.mags)*np.exp(rexp)
            res_w = 4*np.pi*fs*(fe(x, fs, deg)-self.mags)*np.exp(rexp)
            ret = np.zeros(deg//2+2)
            ret[1] = -np.sum(res_w)*df
            for j in range(deg//2):
                ret[0] += (2*j+2)*x[j+2]*np.sum(res_w*(fs-x[0])**(2*j+1))*df
                ret[j+2] = -np.sum(res_w*(fs-x[0])**(2*j+2))*df
            return ret
        #res = opt.minimize(lambda x, fe: np.sum( (fe(x, fs, deg)-self.mags)**2 )*df, x0, method='Nelder-Mead', args=self._fitted_e)
        res = opt.minimize(lambda x, fe: np.sum( (fe(x, fs, deg)-self.mags)**2 )*df, x0, jac=jac_fp, args=self._fitted_e)
        if verbose > 4:
            print("optimized generalized Gaussian x0={}", x0)
            print(jac_fp(x0, self._fitted_e))
            print(res)
        return res, fscale
    def _fit_sym_amps(self, deg):
        '''Fit the frequency space magnitudes to the form 2*pi*f*exp(sum_k a_2k(f-f0)**2k)
        deg: the degree of the polynomial in the exponent to use when fitting'''
        res,sc = self._opt_gauss_envelope(deg)
        return self._fitted_e(res.x, self.freqs*sc, deg), self._fitted_v(res.x, self.freqs*sc, deg)*sc, res.x[0]/sc

    def _est_f0(self, threshold=CUT_ALPHA, mode='odd'):
        '''Estimate the value for f0 by finding the point where the third derivative of the angles vanishes
        returns: a tuple with the frequency and its closest index'''
        f0i = np.argmax(self.mags)
        f0 = self.freqs[f0i]
        if mode == 'odd':
            #find the point to be used as f0. This is the point where the second derivative of arg(E) crosses 0.
            f_min,f_max = self._get_freq_fwhm(f0i=f0i, threshold=threshold)
            angles = fix_angle_seq(np.angle(self.vf[f_min:f_max]))
            f0 = self.f0
            df = self.freqs[1]-self.freqs[0]
            prev_thd_der = (angles[2] - 2*angles[1] + angles[0])/df**2
            best_mag = 0
            for i in range(2, len(angles)-2):
                thd_der = (angles[i+2] - 3*angles[i+1] + 3*angles[i] - angles[i-1])/df**2
                #check for zero crossings
                if (prev_thd_der <= 0 and thd_der >= 0) or (prev_thd_der >= 0 and thd_der <= 0):
                    mn,mx = min(prev_thd_der, thd_der), max(prev_thd_der, thd_der)
                    #take the higher value of |E| as a tiebreaker
                    if self.mags[f_min+i] >= best_mag:
                        f0 = self.freqs[f_min+i]+df/2 if (mx == mn) else self.freqs[f_min+i] - mn*df/(mx-mn)
                        best_mag = self.mags[f_min+i]
            f0i = int( np.rint(f0 / df) )
        elif mode == 'avg':
            #we exclude high frequency noisy data when taking the average
            f0 = np.trapz(self.mags[:2*self.f0_ind]*self.freqs[:2*self.f0_ind])/np.trapz(self.mags[:2*self.f0_ind])
            f0i = int( np.rint(f0/(self.freqs[1]-self.freqs[0])) )
        elif mode == 'avg_v':
            #we exclude high frequency noisy data when taking the average
            vfv = self.mags[:2*self.f0_ind]/(2*np.pi*self.freqs+np.exp(-DIV_S*self.freqs/self.f0))
            f0 = np.trapz(vfv*self.freqs[:2*self.f0_ind])/np.trapz(vfv)
            f0i = int( np.rint(f0/(self.freqs[1]-self.freqs[0])) )
        elif mode == 'fit':
            res,sc = self._opt_gauss_envelope(DEF_DEG)
            f0 = res.x[0]/sc
            f0i = int( np.rint(f0/(self.freqs[1]-self.freqs[0])) )
        if f0i >= len(self.mags):
            f0i = self.mags-1
        print(f0i)
        return f0, f0i

    def get_ang_func(self, f0, fit_axs=None):
        f_min,f_max = self._get_freq_fwhm()
        angles = fix_angle_seq(np.angle(self.vf[f_min:f_max]))
        '''Fit arg(E(f+f0)) to an odd polynomial'''
        #fit to a polynomial with only odd powers + constant
        popt, pcov = opt.curve_fit(lambda f,p0,p1,p3,p5,p7: p7*(f-f0)**7 + p5*(f-f0)**5 + p3*(f-f0)**3 + p1*(f-f0) + p0, self.freqs[f_min:f_max], angles)
        #save a plot of the fit to the axes if supplied
        if fit_axs is not None:
            #setup the axes
            fit_axs.set_xlabel("$f$")
            fit_axs.set_ylabel("|$\tilde{E}(f)|")
            ax2 = fit_axs.twinx()
            ax2.set_ylabel("arg[$\tilde{E}(f)$]", color = 'green')
            ax2.tick_params(axis ='y', labelcolor = 'green')
            #actually plot
            fs = np.linspace(self.freqs[f_min], self.freqs[f_max])
            ax2.scatter(self.freqs[f_min:f_max], angles, color='green')
            ax2.plot(fs, popt[4]*(fs-f0)**7 + popt[3]*(fs-f0)**5 + popt[2]*(fs-f0)**3 + popt[1]*(fs-f0) + popt[0], color='green', linestyle='--')
            fit_axs.plot(self.freqs[f_min:f_max], self.mags[f_min:f_max])
            fit_axs.plot(self.freqs[f_min:f_max], self.mags[f_min:f_max]/(self.freqs[f_min:f_max]+f0))
            fit_axs.axvline(x=f0, color='gray')
        return lambda x: popt[4]*x**7 + popt[3]*x**5 + popt[2]*x**3 + popt[1]*x

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

    def __init__(self, t_pts, v_pts, scan_length=5, lowpass_inc=1.0, noise_thresh=0.1, lowpass_center=-1.0, w0_type='avg', param_axs=None):
        '''Based on a time series sampled at t_pts with values v_pts, try to find a Gaussian pulse envelope.
        t_pts: a numpy array with the times of each sampled point. This must have the same dimensions as v_pts
        v_pts: a numpy array with the field at each sampled time
        lowpass_inc: in the case of a bad fit, a lowpass filter with increasingly high strength is applied. This specifies the increment at which the strength of the lowpass filter increases.
        scan_length: This is used when checking for noise
        noise_thresh: This is the ratio of peak frequency magnitude to high frequency components at which a signal is decided to be noisy.
        max_f0: This specifies the maximum allowable frequency at which the central frequency may occur. If the central frequency is found to be higher than this value, then a lowpass filter is applied.
        w0_type: there are multiple ways of finding the central frequency. This argument should be one of these strings, "odd", "avg", or "peak"
            -odd: select omega_0 such that arg(a(omega)) is close to an odd function. This helps to minimize the imaginary component of the extracted b(t)
            -avg: average over all magnitudes $|\tilde{E}(\omega)|$
            -peak: Take the value of $\omega$ which maximizes $|\tilde{E}(\omega)|$
        lowpass_center: if >0, then lowpass filters are centered around this frequency so as to avoid attenuating genuine signal
        '''
        self.t_pts = t_pts
        self.v_pts = v_pts
        #inferred values
        self.dt = t_pts[1]-t_pts[0]
        self.v_abs = np.abs(v_pts)
        self.t0_ind = self._guess_t0_ind()
        #take the fourier transform       
        self.freqs = fft.rfftfreq(len(v_pts), d=self.dt)
        max_f0_ind = len(self.freqs)//2
        #Shift the pulse so that the peak magnitude is at zero frequency. Taking the inverse fourier transform of this shifted magnitude gives pulse amplitude. The frequency at which the magnitude is maximized is the same as the unmodulated sine wave.
        vf0 = 2*fft.rfft(np.roll(v_pts, -self.t0_ind))
        self.vf = np.copy(vf0) 
        self.mags = np.abs(self.vf)
        self._low_stren = 0.0
        self._last_cent_f = 0.0
        while True:
            self.phi = np.angle(self.vf[0])
            self.f0_ind = np.argmax(self.mags)
            if w0_type == 'peak' or w0_type == 'odd':
                self.f0,self.f0_ind = self._est_f0(mode='peak')
            else:
                self.f0,self.f0_ind = self._est_f0(mode='avg')
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
        self._param_est(axs=param_axs)
        if verbose > 1:
            print("\tf0_corr: {}\n\tt0_corr: {}\n\tphi_corr: {}".format(self.f0, self._t0_corr, self._phi_corr))
        self.envelope = None
        self.envelope_asym = None
        self.asym_lowpass = DEF_LOW_SCALE
        self.asym_f0 = None
        self.peaks_arr = None
        self.t0_ind = np.argmax(self.v_abs)
        self.t0 = self.t_pts[self.t0_ind]
        if w0_type == 'odd':
            self.f0,self.f0_ind = self._est_f0(mode='odd')

    def _get_sym_amps(self, vec_vf, f0):
        #weirdness converting between scipy rfft and fft
        fsize = 2*self.mags.shape[0] - 2 + (self.v_pts.shape[0] % 2)
        ef = np.zeros(fsize)
        vef = np.zeros(fsize)
        df = self.freqs[1]-self.freqs[0]
        f0i = int(f0/df)
        lrp = f0 - f0i*df
        for i in range(min(f0i, self.freqs.shape[0]-f0i-1)):
            #find amplitude for \tilde{b}(\omega) (note the factor of two comes from equation 13
            vef[i] = 0.5*(vec_vf[f0i-i+1]-vec_vf[f0i-i]+vec_vf[f0i+i+1]-vec_vf[f0i+i])*lrp/df + 0.5*(vec_vf[f0i-i]+vec_vf[f0i+i])
            vef[-i] = vef[i]
            #transform from \tilde{b} to \tilde{a}
            ef[-i] = 2*np.pi*(f0-self.freqs[i])*vef[-i]
            ef[i] = 2*np.pi*(f0+self.freqs[i])*vef[i]
        return ef, vef

    def calc_envelope(self, f0=-1, sym_mode="sym"):
        '''Calculate a(t) and b(t)
        f0: the central frequency to use for calculations
        sym_mode: if set to "sym", then perform a symmeterization about f0 by taking 2b(w) = E(w+w0)/(w+w0) + E(w0-w)/(w0-w), if set to "sym-o", do the same, but choose the value for w0 that minimizes the square error. If set to "exp", then fit the amplitude to a generalized Gaussian (see _fit_sym_amps). Otherwise, no symmeterization is performed'''
        if f0 <= 0:
            f0 = self.f0
            f0i = self.f0_ind
        else:
            f0i = int(f0 / (self.freqs[1]-self.freqs[0]))
        pad_n =  (0, self.mags.shape[0]-2+(self.v_pts.shape[0]%2))
        if sym_mode is not None:
            if sym_mode=="exp":
                ef,vf,f0 = self._fit_sym_amps(DEF_DEG)
                f0i = int(f0 / (self.freqs[1]-self.freqs[0]))
                self.env_fourier = np.roll(np.pad(ef, pad_n), -f0i)
                self.vec_env_fourier = np.roll(np.pad(vf, pad_n), -f0i)
            else:
                #start by finding b_0(omega-omega_0)
                vec_vf = self.mags/(2*np.pi*self.freqs+np.exp(-DIV_S*self.freqs/self.f0))
                if sym_mode == "sym-o":
                    vw,ew = 1,1
                    def fp(f0, af):
                        ef,vef = af(vec_vf, f0[0])
                        rr = vw*np.sum( (np.abs(vef)-np.abs(self.vec_env_fourier))**2 ) + ew*np.sum( (np.abs(ef)-np.abs(self.env_fourier))**2 )
                        return rr 
                    res = opt.minimize(fp, f0, args=self._get_sym_amps, method='Nelder-Mead')
                    if verbose > 4:
                        print("optimized", res)
                    #account for phases
                    f0 = res.x[0]
                self.env_fourier, self.vec_env_fourier = self._get_sym_amps(vec_vf, f0)
            #get angles
            ang_fun = self.get_ang_func(f0)
            vp = ang_fun(fft.fftfreq(len(self.v_pts), d=self.dt))
            #note that we don't need to account for phi_corr here since ang_fun() already removed this part
            self.env_fourier = self.env_fourier*np.exp(1j*vp)
            self.vec_env_fourier = self.vec_env_fourier*np.exp(1j*vp)
        else:
            self.env_fourier = np.roll(np.pad(self.vf, pad_n), -f0i)*np.exp(-1j*self._phi_corr)
            self.vec_env_fourier = np.roll(np.pad(self.vf/(2*np.pi*self.freqs+np.exp(-DIV_S*self.freqs/self.f0)), pad_n), -f0i)*np.exp(-1j*self._phi_corr)
        #we want to ensure self consistency between envelopes and central frequencies
        self.f0 = f0
        self.f0_ind = f0i
        #take the inverse fourier transform and shift to the peak of the pulse
        self.envelope = 2*np.roll(fft.ifft(self.env_fourier), self.t0_ind)
        '''sig_direct = np.real(self.envelope*np.exp(1j*(2*np.pi*self.f0*(self.t_pts-self.t_pts[self.t0_ind]) + self._phi_corr)))
        t_ind_diff = np.argmax(self.v_pts) - np.argmax(np.real(sig_direct))
        self.envelope = np.roll(self.envelope, t_ind_diff)
        self.t0 = t_ind_diff*self.dt + self.t_pts[self.t0_ind]
        print("\tpeak shift: {}".format(np.sign(t_ind_diff)*self.t_pts[np.abs(t_ind_diff)]))'''
        #for odd numbers of time points, the inverse fourier transform will have one fewer points
        if len(self.t_pts) > len(self.envelope):
            self.envelope = np.pad( self.envelope, (0,len(self.t_pts)-len(self.envelope)) )

    def get_envelope(self):
        if self.envelope is None:
            self.calc_envelope()
        return self.envelope

    def get_mu(self, env=None):
        if env is None:
            env = np.abs( self.get_envelope() )
            env = env / (np.trapz(env)*self.dt)
        mu = np.trapz(env*self.t_pts)*self.dt
        return mu, env

    def get_sig(self, env=None):
        mu, env = self.get_mu(env=env)
        sig = np.sqrt( np.trapz(env*(self.t_pts-mu)**2)*self.dt )
        return mu, sig, env

    def get_skewness(self, env=None):
        #Calculate mu and sigma. Use these to get the skewness
        mu, sig, env = self.get_sig(env=env)
        skew = np.trapz(env*(self.t_pts - mu)**3)/(sig**3)
        return mu, sig, skew, env

    def plt_raw_fdom(self, axs):
        # Add the magnitude 
        axs.set_xlabel('$\omega$ (1/fs)') 
        axs.set_ylabel('|$E(\omega)$|', color = 'black') 
        axs.plot(2*np.pi*self.freqs, np.abs(self.vf), color = 'black', label='magnitude')
        #plot the bounds used for fitting
        frange = self._get_freq_fwhm()
        axs.axvline(2*np.pi*self.freqs[frange[0]], color='gray', linestyle=':')
        axs.axvline(2*np.pi*self.freqs[frange[1]], color='gray', linestyle=':')
        axs.axvline(2*np.pi*self.f0, color='gray')
        axs.tick_params(axis ='y', labelcolor = 'black')
        axs.set_ylim(0, max(self.mags)*1.1)
        # Adding Twin Axes
        ax2 = axs.twinx()
        ax2.set_ylim(-1, 1)
        ax2.set_ylabel("arg[$E(\omega)$]$/\pi$", color = 'green')
        ax2.tick_params(axis ='y', labelcolor = 'green')
        ax2.plot(2*np.pi*self.freqs, np.angle(self.vf)/np.pi, color='green', label='phase')
    def plt_raw_tdom(self, axs):
        peak_t = self.t_pts[self.t0_ind] - self._t0_corr
        #get the envelope and perform a fitting
        env = self.get_envelope()
        axs.plot(self.t_pts, self.v_pts, color='black', label='simulated data')
        axs.plot(self.t_pts, env, color='teal', label='envelope')
        axs.plot(self.t_pts, self.get_envelope(), color='red', label='envelope')
        axs.axvline(peak_t, color='teal', linestyle=':')
        axs.plot(self.t_pts, env*np.cos(2*np.pi*self.f0*(self.t_pts-peak_t)+self._phi_corr), color='orange', label='extracted pulse')

    def compare_envelopes(self, axs):
        env = self.get_envelope()
        axs.plot(self.t_pts, self.v_pts, color='black', label='measured E(t)')
        axs.fill_between(self.t_pts, np.real(env), -np.real(env), color='blue', label='Re$[a(t)]$', alpha=0.2)
        axs.fill_between(self.t_pts, np.imag(env), -np.imag(env), color='red', label='Im$[a(t)]$', alpha=0.2)
        #axs.legend(loc='upper right')

    def compare_signals(self, axs):
        print(self._phi_corr)
        env = self.get_envelope()
        max_field_t = self.t_pts[np.argmax(self.v_pts)]
        #get the signal from a(t)e^(i(\omega_0(t-t_0) + \phi)) + c.c.
        sig_direct_re = np.real(env)*np.cos(2*np.pi*self.f0*(self.t_pts-self.t0) + self._phi_corr)
        sig_direct_im = np.imag(env)*np.sin(2*np.pi*self.f0*(self.t_pts-self.t0) + self._phi_corr)
        axs.plot(self.t_pts, self.v_pts, color='black', label='measured')
        env_vec = 2*np.roll( fft.ifft(self.vec_env_fourier), int(self.t0/self.dt) )
        axs.fill_between(self.t_pts, np.real(env), -np.real(env), color='green', label='Re$[a(t)]$', alpha=0.2)
        axs.fill_between(self.t_pts, np.imag(env), -np.imag(env), color='blue', label='Im$[a(t)]$', alpha=0.2)
        #f0,_ = self._est_f0(mode='avg')
        f0 = self.f0
        fit_ser = np.diff( np.real(env_vec)*np.sin(2*np.pi*f0*(self.t_pts-self.t0)+self._phi_corr) )/self.dt
        fit_ser = np.roll(fit_ser, -np.argmax(fit_ser)+np.argmax(self.v_pts))
        axs.plot(self.t_pts[:-1], fit_ser*np.max(self.v_pts)/np.max(fit_ser), color='red', label='extracted', linestyle='-.')
        axs.set_xlim(self.t_pts[0], self.t_pts[-1])

    def compare_fspace(self, axs, sym_mode='sym-o', plt_raxs=True):
        full_freqs = fft.fftfreq(len(self.t_pts))
        #make unsymmeterized envelope plots
        self.calc_envelope(sym_mode=None)
        env_unsym = self.get_envelope()
        axs.plot(fft.fftshift(full_freqs), fft.fftshift(np.abs(self.env_fourier)), color='blue', label='|$a(\omega)$|')
        axs.plot(fft.fftshift(full_freqs), fft.fftshift(np.abs(self.vec_env_fourier)), color='orange', label='|$b(\omega)$|')
        #plot the symmeterized envelope for comparison
        self.calc_envelope(sym_mode=sym_mode)
        env = self.get_envelope()
        axs.plot(fft.fftshift(full_freqs), fft.fftshift(np.abs(self.env_fourier)), color='blue', linestyle='--')
        axs.plot(fft.fftshift(full_freqs), fft.fftshift(np.abs(self.vec_env_fourier)), color='orange', linestyle='--')
        ax2 = axs.twinx()
        ax2.tick_params(axis ='y', labelcolor = 'green')
        ext_angles = fft.fftshift(np.angle(self.vec_env_fourier))
        ax2.plot(fft.fftshift(full_freqs), [fix_angle(a)/np.pi for a in ext_angles], color='green', linestyle='--')
        #figure out the actual angles for comparison
        pad_n = (0, self.mags.shape[0]-2+self.v_pts.shape[0]%2)
        actual_angles = fft.fftshift( np.roll(np.pad(np.angle(self.vf)-self._phi_corr, pad_n), -self.f0_ind) )
        ax2.plot(fft.fftshift(full_freqs), [fix_angle(a)/np.pi for a in actual_angles], color='green', label='arg[$a(\omega)$]')
        ax2.set_ylim(-1, 1)
        if not plt_raxs:
            ax2.get_yaxis().set_visible(False)

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
    def __init__(self, fname, pass_alpha=1.0, slice_dir='x', prefix='.', keep_n=2.5, scan_length=5):
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
        err_2 = 0.02
        return np.array(self.f[clust][points[ind]]['time']['Re']), err_2

    def opt_pulse_full(self, t_pts, a_pts, err_sq, fig_name='', raw_ax=None, final_ax=None):
        if a_pts.shape != t_pts.shape:
            raise ValueError("t_pts and a_pts must have the same shape")
        if t_pts.shape[0] == 0:
            raise ValueError("empty time series supplied!")

        #do the signal processing to find the envelope and decompose it into a series of peaks
        psig = signal(t_pts, a_pts, lowpass_inc=2.0)
        w0 = 2*np.pi*psig.f0
        phi = psig._phi_corr
        #from the envelope find the amplitude and phase
        sig_env = psig.get_envelope()
        max_env_ind = np.argmax(sig_env)
        amp = sig_env[max_env_ind]
        #store the information
        res = extract_info(amp, w0, phi)
        if fig_name != '':
            res_name = fig_name+"_res.txt"
            write_reses(res_name, [res])
        if raw_ax is not None:
            psig.compare_fspace(raw_ax[0])
            psig.compare_signals(raw_ax[1])
        _,_,skew,_ = psig.get_skewness() 
        return res, skew

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
