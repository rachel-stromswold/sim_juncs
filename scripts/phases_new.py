import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.stats import linregress
from scipy.special import hermite
#from scipy.fft import ifft, fft, fftfreq
import scipy.fft as fft
import scipy.signal as ssig
import time
import pickle
import os.path
import matplotlib.pyplot as plt

'''     fitting paramaters       '''
HERM_N = 3 #use the first HERM_N even hermite polynomials (i.e. HERM_N=1 uses just H_e0(w), HERM_N=2 adds H_e2(w) and so on) 
HERMS = [hermite(n) for n in range(2*HERM_N+2)]
HERM_SCALE = np.sqrt(2)
HERM_SCALE=1
PLT_COLORS = [
    "#016876",
    "#cf4f4f",
    "#97b47d",
    "#e49c4f",
    "#a5538d",
    "#8dbcf9",
    "#8b4a5f",
    "#ffa8ff"
]
PHASE_STDS = 2

verbose = 5

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

def fix_angle_seq(angles, center_ind=0):
    '''Given a sequence of angles from -pi to pi look for -pi -> pi rollovers and perform adjustments such that each angle is the closest equivalent angle to the previous.'''
    shift = 0
    n = len(angles)
    #fix all angles right of center
    for i in range(center_ind+1,n):
        dist_raw = abs(angles[i] - angles[i-1])
        if abs(angles[i]+2*np.pi - angles[i-1]) < dist_raw:
            for j in range(i,n):
                angles[j] += 2*np.pi
        elif abs(angles[i]-2*np.pi - angles[i-1]) < dist_raw:
            for j in range(i,n):
                angles[j] -= 2*np.pi
    if center_ind <= 0:
        return angles
    #fix all angles left of center
    for i in range(center_ind-1, -1, -1):
        dist_raw = abs(angles[i] - angles[i+1])
        if abs(angles[i]+2*np.pi - angles[i+1]) < dist_raw:
            for j in range(i+1):
                angles[j] += 2*np.pi
        elif abs(angles[i]-2*np.pi - angles[i+1]) < dist_raw:
            for j in range(i+1):
                angles[j] -= 2*np.pi
    return angles

def get_fwhm(pulse, threshold):
    '''
    Given a strictly positive time or frequency series pulse, estimate its center and standard deviation
    pulse: the series for which we find information
    threshold: the widths are the smallest distance from the peak such that pulse[i] >= threshold*pulse[max]
    returns: a tuple with the peak index, lower cutoff index and higher cutoff index
    '''
    f0i = np.argmax(pulse)
    cut_amp = pulse[f0i]*threshold
    f_min, f_max = f0i, f0i
    bound = min(f0i, len(pulse)-f0i-1)
    for j in range(bound):
        if pulse[f0i - j] > cut_amp:
            f_min = f0i - j
        if pulse[f0i + j] > cut_amp:
            f_max = f0i + j
    return f0i, f_min, f_max

class signal:
    def __init__(self, t_pts, v_pts, fname=None, thin=1, phase_axs=None):
        self.dt = t_pts[1] - t_pts[0]
        self.t_pts = t_pts
        self.v_pts = v_pts
        #estimate the central time
        self.shift_i = np.argmax(np.abs(v_pts))
        self._fit_shift = t_pts[self.shift_i]
        #take the fourier transform       
        freqs = fft.rfftfreq(len(v_pts), d=self.dt)
        vf0 = 2*fft.rfft(np.roll(v_pts, -self.shift_i))
        self.vfm = np.abs(vf0)
        #estimate the central frequency and envelope width in frequency space
        guess_f0i, lo_fi, hi_fi = get_fwhm(self.vfm, threshold=np.exp(-PHASE_STDS))
        guess_f0 = freqs[guess_f0i]
        guess_sigma = (freqs[hi_fi]-freqs[lo_fi])/np.sqrt(2*PHASE_STDS)
        self.fit_bounds = np.array([lo_fi, hi_fi])
        self.vfa = fix_angle_seq(np.angle(vf0), center_ind=guess_f0i)
        #obtain a better estimate of f0 using averaging
        guess_f0 = np.trapz(self.vfm[:2*guess_f0i]*freqs[:2*guess_f0i])/np.trapz(self.vfm[:2*guess_f0i])
        guess_f0i = int( np.rint(guess_f0/(freqs[1]-freqs[0])) )
        #now using the selected frequency range, estimate the phase by performing linear regression
        lr_res = linregress(freqs[lo_fi:hi_fi], self.vfa[lo_fi:hi_fi])
        guess_t0 = -lr_res.slope
        guess_phi = fix_angle(lr_res.intercept)
        #now cancel out the f*t0 rotation
        #vf0 *= np.exp( 2j*np.pi*freqs*guess_t0)
        self.vfr = np.real(vf0)
        self.vfi = np.imag(vf0)
        #finally we perform optimization
        def residuals(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],...
            dd = x[3]*(freqs-x[2])
            emags = np.zeros(freqs.shape)
            for i in range(HERM_N):
                emags += x[i+4]*HERMS[2*i](dd)
            '''emags *= 0.5*freqs*np.exp(2j*np.pi*(x[0]-freqs*x[1]) - dd**2)
            #return np.sum( (np.exp(-dd**2)*emags - self.vfm)**2 + (x[0] - freqs*x[1] - self.vfa)**2 )
            return np.sum( (np.real(emags)-self.vfr)**2 + (np.imag(emags)-self.vfi)**2 )'''
            emags *= 0.5*freqs*np.exp(-dd**2)
            return np.sum( emags**2 + self.vfm**2 - 2*emags*self.vfm*np.cos(2*np.pi*(x[0]-freqs*x[1])-self.vfa) )
        def grad_res(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],...
            dd = x[3]*(freqs-x[2])
            emags = np.zeros(freqs.shape[0])
            eargs = 2*np.pi*(x[0] - freqs*x[1])
            demags = np.zeros(freqs.shape[0])
            cemags = np.zeros((HERM_N, freqs.shape[0]))
            for i in range(HERM_N):
                emags += freqs*x[i+4]*HERMS[2*i](dd)*np.exp(-dd**2)
                cemags[i,:] = freqs*HERMS[2*i](dd)*np.exp(-dd**2)
                demags += freqs*x[i+4]*HERMS[2*i+1](dd)*np.exp(-dd**2)
            ret = np.zeros(HERM_N+4)
            mag_as = emags*self.vfm*np.sin(eargs-self.vfa)
            mag_gs = emags - 2*self.vfm*np.cos(eargs-self.vfa)
            ret[0] = 2*np.pi*np.sum( mag_as )
            ret[1] = -2*np.pi*np.sum( freqs*mag_as )
            ret[2] = 0.5*np.sum( mag_gs*demags*x[3] )
            ret[3] = -0.5*np.sum( mag_gs*demags*(freqs-x[2]) )
            for i in range(HERM_N):
                ret[i+4] = 0.5*np.sum( mag_gs*cemags[i,:] )
            return ret
        '''def residuals(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],...
            dl,dr = x[3]*(freqs+x[2]), x[3]*(freqs-x[2])
            l,r = np.zeros(freqs.shape), np.zeros(freqs.shape)
            for i in range(HERM_N):
                l += x[i+4]*HERMS[2*i](HERM_SCALE*dl)
                r += x[i+4]*HERMS[2*i](HERM_SCALE*dr)
            l = 0.5*freqs*( r*np.exp(2j*np.pi*(x[0]-freqs*x[1]) - dr**2) - l*np.exp(-2j*np.pi*(x[0]+freqs*x[1]) - dl**2) )
            return np.sum( (np.real(l) - self.vfr)**2 + (np.imag(l) - self.vfi)**2 )'''
        x0 = np.append(np.array([guess_phi, guess_t0, guess_f0, 1/guess_sigma]), np.zeros(HERM_N))
        x0[4] = self.vfm[guess_f0i]/freqs[guess_f0i]
        '''herm_den = np.sqrt(np.pi)
        for i in range(HERM_N):
            dd = x0[3]*(freqs-x0[2])
            x0[i+4] = np.trapz(self.vfm*HERMS[2*i](dd)*np.exp(-dd**2), dx=freqs[1]-freqs[0])/herm_den
            herm_den *= 4*(i+1)*(i+2)'''
        #test that analytic gradients look correct
        if (verbose > 4):
            num_g = np.zeros(HERM_N+4)
            for i in range(HERM_N+4):
                ei = np.zeros(HERM_N+4)
                ei[i] = 0.001
                num_g[i] = ( residuals(x0 + ei) - residuals(x0 - ei) )/ei[i]/2
            print("gradients at x0:")
            print("\tnumeric:\t", num_g)
            print("\tanalytic:\t", grad_res(x0))
        opt_res = opt.minimize(residuals, x0, jac=grad_res)
        #opt_res.x = x0
        if verbose > 4:
            print("x0:\t", x0)
            print("f(x0):\t", residuals(x0))
            print("optimization result")
            print(opt_res)
        #use the results
        self.phi = opt_res.x[0]
        self.t0 = opt_res.x[1] + self._fit_shift
        self.f0 = opt_res.x[2]
        self.sigma = 1/opt_res.x[3]
        self.herm_coeffs = opt_res.x[4:]
        #generate plots
        if phase_axs is not None:
            dat_series = [fix_angle(ang)/np.pi for ang in self.vfa[lo_fi:hi_fi]]
            fit_series_0 = (guess_t0*freqs[lo_fi:hi_fi] + guess_phi)/np.pi
            fit_series_1 = (self.t0*freqs[lo_fi:hi_fi] + self.phi)/np.pi
            dat_omegas = 2*np.pi*freqs[lo_fi:hi_fi]
            phase_axs.scatter(dat_omegas, dat_series, color='black', label="simulation")
            phase_axs.plot(dat_omegas, fit_series_0, color=PLT_COLORS[0], linestyle='--', label="linear regression")
            phase_axs.plot(dat_omegas, fit_series_1, color=PLT_COLORS[1], linestyle='--', label="full optimization")
            phase_axs.annotate('$\\varphi = ${:.2f}, $t_0 = ${:.2f} fs\n$R^2 = ${:.2f}'.format(self.phi, self.t0, lr_res.rvalue**2), xy=(0,0), xytext=(0.2, 0.80), xycoords='figure fraction')
            phase_axs.set_ylim(-1,1)

    def get_fenv(self, freqs, field='E'):
        '''
        Get the frequency domain envelope for the electric field (or vector potential)
        freqs: the frequencies for which the field is computed
        field: if set to 'E' (default) compute the electric field, otherwise compute the envelope for the vector potential
        '''
        dd = freqs/self.sigma
        emags = np.zeros(freqs.shape, dtype=complex)
        for i in range(HERM_N):
            emags += self.herm_coeffs[i]*HERMS[2*i](dd)
        emags *= 0.5*np.exp(-2j*np.pi*freqs*self.t0 - dd**2)
        if field == 'E':
            return (freqs+self.f0)*emags
        return emags
        '''dd = (freqs - self.f0)/self.sigma
        #first compute the hermite polynomials
        emags = np.zeros(freqs.shape)
        for i in range(HERM_N):
            emags += self.herm_coeffs[i]*HERMS[2*i](HERM_SCALE*dd)
        #convert from vector potential to electric field
        if field == 'E':
            emags *= freqs
        #TODO: figure out imaginary component
        #reutrn the sum of the polynomials times the gaussian
        return emags*np.exp(2j*np.pi*(self.phi-freqs*self.t0) - dd**2)'''

    def get_tenv(self, ts, field='E'):
        '''
        Get the time domain envelope for the electric field (or vector potential)
        freqs: the frequencies for which the field is computed
        field: if set to 'E' (default) compute the electric field, otherwise compute the envelope for the vector potential
        '''
        freqs = fft.fftfreq(len(ts), d=ts[1]-ts[0])
        fenv = self.get_fenv(freqs, field=field)#*np.exp(2j*np.pi*(self.phi+freqs*self._fit_shift))
        return fft.ifft(fenv)
        '''
        dd = (freqs - self.f0)/self.sigma
        #first compute the hermite polynomials
        emags = np.zeros(freqs.shape)
        for i in range(HERM_N):
            emags += self.herm_coeffs[i]*HERMS[2*i](HERM_SCALE*dd)
        #convert from vector potential to electric field
        if field == 'E':
            emags *= (freqs + self.f0)
        #TODO: figure out imaginary component
        #reutrn the sum of the polynomials times the gaussian
        return np.roll( fft.ifft(emags*np.exp(-dd**2)), int((self.t0-ts[0])/(ts[1]-ts[0])) )'''

    def plt_raw_fdom(self, axs):
        freqs = fft.rfftfreq(len(self.v_pts), d=self.dt)
        # setup labels and annotations
        axs.set_xlabel('$\omega$ (1/fs)') 
        axs.set_ylabel('|$E(\omega)$|', color = 'black') 
        axs.axvline(2*np.pi*freqs[self.fit_bounds[0]], color='gray', linestyle=':')
        axs.axvline(2*np.pi*freqs[self.fit_bounds[1]], color='gray', linestyle=':')
        axs.axvline(2*np.pi*self.f0, color='gray')
        axs.tick_params(axis ='y', labelcolor = 'black')
        # Adding Twin Axes
        '''ax2 = axs.twinx()
        ax2.set_ylim(-1, 1)
        ax2.set_ylabel("arg[$E(\omega)$]$/\pi$", color='green')
        ax2.tick_params(axis ='y', labelcolor='green')
        #add the magnitude
        axs.plot(2*np.pi*freqs, np.abs(self.vfm), color = 'black', label='magnitude')
        axs.plot(2*np.pi*freqs, np.abs(self.get_fenv(freqs)), color=PLT_COLORS[0])
        ax2.plot(2*np.pi*freqs, self.vfa/np.pi, color='green', label='phase')'''
        #get the envelope function and multiply it by the appropriate rotation. By definition, the envelope is centered at zero frequency, so we have to roll the envelope back to its original position
        fit_field = self.get_fenv(freqs-self.f0)*np.exp(2j*np.pi*(self.phi + freqs*self._fit_shift))
        axs.plot(2*np.pi*freqs, self.vfm, color='black')
        axs.plot(2*np.pi*freqs, np.abs(fit_field), color='black', linestyle='-.')
        axs.plot(2*np.pi*freqs, self.vfr, color=PLT_COLORS[0], label='real')
        axs.plot(2*np.pi*freqs, self.vfi, color=PLT_COLORS[1], label='imaginary')
        axs.fill_between(2*np.pi*freqs, np.real(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[0], alpha=0.2)
        axs.fill_between(2*np.pi*freqs, np.imag(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[1], alpha=0.2)

    def plt_raw_tdom(self, axs):
        peak_t = self.t_pts[int(self.t0/self.dt)] - self.t0
        #get the envelope and perform a fitting
        axs.plot(self.t_pts, self.v_pts, color=PLT_COLORS[0], label='simulated data')
        '''freqs = fft.rfftfreq(len(self.v_pts), d=self.dt)
        fit_field = self.get_fenv(freqs-self.f0)*np.exp(2j*np.pi*(self.phi + freqs*self._fit_shift))
        axs.plot(self.t_pts, np.roll(fft.irfft(fit_field), self.shift_i), color=PLT_COLORS[1], label='extracted')'''
        env = self.get_tenv(self.t_pts)
        axs.plot(self.t_pts, env, color=PLT_COLORS[1], label='envelope')
        axs.axvline(peak_t, color='teal', linestyle=':')
        axs.plot(self.t_pts, env*np.cos(2*np.pi*self.f0*(self.t_pts-self.t0)+self.phi), color='orange', label='extracted pulse')

    def compare_fspace(self, axs, sym_mode='sym-o', plt_raxs=True):
        freqs = fft.fftfreq(len(self.v_pts), d=self.dt)
        #make unsymmeterized envelope plots
        axs.plot(fft.fftshift(freqs), fft.fftshift(np.abs(self.get_fenv(freqs, field='E'))), color=PLT_COLORS[0], label='|$a(\omega)$|')
        axs.plot(fft.fftshift(freqs), fft.fftshift(np.abs(self.get_fenv(freqs, field='A'))), color=PLT_COLORS[1], label='|$a(\omega)$|')
        #plot the symmeterized envelope for comparison
        '''self.calc_envelope(sym_mode=sym_mode)
        env = self.get_envelope()
        axs.plot(fft.fftshift(freqs), fft.fftshift(np.abs(self.env_fourier)), color='blue', linestyle='--')
        axs.plot(fft.fftshift(freqs), fft.fftshift(np.abs(self.vec_env_fourier)), color='orange', linestyle='--')'''
        '''ax2 = axs.twinx()
        ax2.tick_params(axis ='y', labelcolor='green')
        ext_angles = fft.fftshift(np.angle(self.vec_env_fourier))
        ax2.plot(fft.fftshift(freqs), [fix_angle(a)/np.pi for a in ext_angles], color='green', linestyle='--')
        #figure out the actual angles for comparison
        pad_n = (0, self.vfa.shape[0]-2+self.v_pts.shape[0]%2)
        actual_angles = fft.fftshift( np.roll(np.pad(np.angle(self.vf)-self.phi, pad_n), -self.f0) )
        ax2.plot(fft.fftshift(freqs), [fix_angle(a)/np.pi for a in actual_angles], color='green', label='arg[$a(\omega)$]')
        ax2.set_ylim(-1, 1)
        if not plt_raxs:
            ax2.get_yaxis().set_visible(False)'''

    def compare_tspace(self, axs):
        freqs = fft.fftfreq(len(self.v_pts), d=self.dt)
        env = self.get_fenv(freqs, field='E')
        max_field_t = self.t_pts[np.argmax(self.v_pts)]
        #get the signal from a(t)e^(i(\omega_0(t-t_0) + \phi)) + c.c.
        sig_direct_re = np.real(env)*np.cos(2*np.pi*self.f0*(self.t_pts-self.t0) + self.phi)
        sig_direct_im = np.imag(env)*np.sin(2*np.pi*self.f0*(self.t_pts-self.t0) + self.phi)
        axs.plot(self.t_pts, self.v_pts, color='black', label='simulated')
        #first plot the envelope
        env = self.get_tenv(self.t_pts, field='E')
        axs.fill_between(self.t_pts, np.real(env), -np.real(env), color=PLT_COLORS[0], label='Re$[a(t)]$', alpha=0.2)
        axs.fill_between(self.t_pts, np.imag(env), -np.imag(env), color=PLT_COLORS[1], label='Im$[a(t)]$', alpha=0.2)
        #now plot the wave
        env = self.get_tenv(self.t_pts, field='A')
        fit_ser = np.diff( np.real(env)*np.sin(self.f0*(self.t_pts-self.t0)+self.phi) )/self.dt
        #fit_ser = np.roll(fit_ser, self.shift_i)
        axs.plot(self.t_pts[:-1], fit_ser*np.max(self.v_pts)/np.max(fit_ser), color='red', label='extracted', linestyle='-.')
        axs.set_xlim(self.t_pts[0], self.t_pts[-1])

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
    def get_err_sq(self):
        return self.res_arr[2*N_RES_PARAMS]

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
        phi = psig.phi
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
            psig.compare_tspace(raw_ax[1])

        return res

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
            res = self.opt_pulse_full(self.t_pts, np.real(v_pts), err_2, fig_name=fig_name, raw_ax=raw_ax, final_ax=fin_ax)
            #save the point
            n_evals += 1
            good_js.append(j)
            ret.set_point(j, res, err_2)
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
