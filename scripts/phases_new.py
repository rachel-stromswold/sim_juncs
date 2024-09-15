import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.stats import linregress
from scipy.special import hermite
#from scipy.fft import ifft, fft, fftfreq
import scipy.fft as fft
import scipy.signal as ssig
import matplotlib.pyplot as plt

'''     fitting paramaters       '''
EPSILON = 0.01
HERM_N = 20 #use the first HERM_N even hermite polynomials (i.e. HERM_N=1 uses just H_e0(w), HERM_N=2 adds H_e2(w) and so on) 
HERMS = [hermite(n) for n in range(2*HERM_N+2)]
HERM_SCALE = np.sqrt(2)
HERM_SCALE=1
HERM_OFF = 4
POLY_FIT_DEVS = 8 #the number of gaussian standard deviations to take when fitting the angles
ANG_POLY_N = 1 #in addition to phi and t0, consider higher order odd terms in the polynomial for the angle. i.e. arg(a(w)) = \phi - \omega t_0 + \sum_{n=1}^{ANG_POLY_N} \omega^{2n+1} t_n

NOISY_N_PEAKS = 2

verbose = 5

'''
To find starting Hermite-Gaussian coefficients, we use the orthoganality condition $\int H_n(a(x-x0))H_m(a(x-x0))e^{-a^2(x-x0)^2{ dx = \sqrt{pi} 2^n n! \delta_{nm}/a$. This lets us expand $|E(w)| = \sum_n k_n H_n(a(w-w0))e^{-a^2(w-w0)^2}$. However, we need terms of the form $|E(w)| = \sum_n c_{2n} H_n(a(w-w0))e^{-a^2(w-w0)^2}$ to ensure that the DC component vanishes. This function constructs a matrix, A, that takes $\bm{k} = A\bm{c}$. A will always be singular, so a small epsilon*w0 correction is added to the odd diagonals so that A^-1 exists at the cost of no longer making the coefficients exact. This is a small price since we plug the resulting c coefficients into gradient descent.
a: the scale factor for the Hermite Gaussian. This should be a good guess in order to get sane coefficients
w0: the central frequency for the Hermite Gaussian
nn: the number of Hermite polynomials to include in the expansion
returns: A^-1
'''
def herm_conv(w0, a, nn):
    c_to_k = w0*np.identity(2*nn)
    for n in range(nn):
        c_to_k[2*n+1, 2*n+1] *= EPSILON
        if n > 0:
            c_to_k[2*n-1, 2*n] = 2*n/a
        c_to_k[2*n+1, 2*n] = 0.5/a
    return np.linalg.inv(c_to_k)

#returns the angle theta transformed to be in the range (-pi,pi]
def fix_angle(theta):
    try:
        for i, t in enumerate(theta):
            theta[i] = fix_angle(t)
    except TypeError:
        if theta > np.pi:
            theta -= 2*np.pi*np.floor(theta/(2*np.pi))
            if theta > np.pi:
                theta -= 2*np.pi
        elif theta <= -np.pi:
            theta += 2*np.pi*(np.floor(-theta/(2*np.pi)) + 1)
            if theta > np.pi:
                theta -= 2*np.pi
    return theta

def fix_angle_seq(angles, center_ind=0, scan_n=1):
    '''Given a sequence of angles from -pi to pi look for -pi -> pi rollovers and perform adjustments such that each angle is the closest equivalent angle to the previous.'''
    shift = 0
    n = len(angles)
    shiftrange = np.arange(-scan_n, scan_n+1)
    #fix all angles right of center
    for i in range(center_ind+1,n):
        #find the equivalent angle which minimizes the jump
        shift_n = shiftrange[np.argmin(np.abs(angles[i]+2*np.pi*shiftrange - angles[i-1]))]
        #shift by the minimizing angle
        if shift_n != 0:
            for j in range(i,n):
                angles[j] += 2*np.pi*shift_n
    if center_ind <= 0:
        return angles
    #fix all angles left of center
    for i in range(center_ind-1, -1, -1):
        shift_n = shiftrange[np.argmin(np.abs(angles[i]+2*np.pi*shiftrange - angles[i+1]))]
        #shift by the minimizing angle
        if shift_n != 0:
            for j in range(i+1):
                angles[j] += 2*np.pi*shift_n
    return angles

'''
Given a cost function cost_fn and an analytic gradient grad_fn, check that the two gradients approximately match
'''
def test_grads(cost_fn, grad_fn, x):
    if verbose < 1:
        return
    an_g = grad_fn(x)
    num_g = np.zeros(x.shape)
    for i in range(x.shape[0]):
        ei = np.zeros(x.shape[0])
        ei[i] = 0.001
        num_g[i] = ( cost_fn(x + ei) - cost_fn(x - ei) )/ei[i]/2
        if num_g[i] and np.abs( (num_g[i]-an_g[i])/num_g[i] ) > 0.01:
            print("Warning! disagreement between analytic and numeric gradients on axis", i)
    if (verbose > 4): 
        print("gradients at x0:")
        print("\tnumeric:\t", num_g)
        print("\tanalytic:\t", an_g)

class signal:
    @staticmethod
    def _fourier_env_formx(fs, x, herm_n=3, ang_n=1, calc_grads=False):
        dd = fs - abs(x[2])
        eargs = x[0] - fs*x[1]
        ang_off = HERM_OFF+herm_n
        for m in range(ang_n):
            eargs += x[ang_off+m]*dd**(2*m+3)
        #eargs *= 2*np.pi
        emags = np.zeros(fs.shape[0])
        mag_grads = None
        ang_grads = None
        if calc_grads:
            #calculate angle gradients
            ang_grads = np.zeros((x.shape[0], fs.shape[0]))
            ang_grads[0,:] = 1 
            ang_grads[1,:] = -fs
            for m in range(ang_n):
                ang_grads[ang_off+m,:] = dd**(2*m+3)
                ang_grads[2,:] -= np.sign(x[2])*(2*m+3)*x[m+ang_off]*dd**(2*m+2)
                #ang_grads[3,:] += (2*m+3)*x[m+ang_off]*(fs-abs(x[2]))*dd**(2*m+2)
            #calculate magnitude gradients
            dd *= x[3]
            mag_grads = np.zeros((x.shape[0], fs.shape[0]))
            for m in range(herm_n):
                emags += x[HERM_OFF+m]*HERMS[2*m](dd)
                mag_grads[2,:] -= fs*np.sign(x[2])*x[3]*x[HERM_OFF+m]*(dd*HERMS[2*m](dd) - HERMS[2*m+1](dd))*np.exp(-dd**2/2)
                mag_grads[3,:] += fs*(fs-abs(x[2]))*x[HERM_OFF+m]*(dd*HERMS[2*m](dd) - HERMS[2*m+1](dd))*np.exp(-dd**2/2)
                mag_grads[HERM_OFF+m,:] = fs*HERMS[2*m](dd)*np.exp(-dd**2/2)
        else:
            dd *= x[3]
            for m in range(herm_n):
                emags += x[HERM_OFF+m]*HERMS[2*m](dd)
        emags *= np.exp(-dd**2/2)
        return emags, eargs, mag_grads, ang_grads

    @staticmethod
    def _do_grad_tests(fs, herm_n, ang_n, n_tests=10):
        print("\n====================================\n")
        ang_off = HERM_OFF+herm_n
        def wrap_eargs(x):
            _, eargs, _, _ = signal._fourier_env_formx(fs, x)
            return np.sum(eargs)
        def grad_wrap_eargs(x):
            _, _, _, ang_grad = signal._fourier_env_formx(fs, x, calc_grads=True)
            return np.sum(ang_grad, axis=1)
        for i in range(n_tests):
            x0 = np.random.random(ang_off+ang_n)*2 - 1
            print("testing argument gradients at ", x0, "(fs =", fs, ")")
            test_grads(wrap_eargs, grad_wrap_eargs, x0)
        print("\n====================================\n")
        def wrap_emags(x):
            emags, _, _, _ = signal._fourier_env_formx(fs, x)
            return np.sum(fs*emags)
        def grad_wrap_emags(x):
            _, _, mag_grad, _ = signal._fourier_env_formx(fs, x, calc_grads=True)
            return np.sum(mag_grad, axis=1)
        for i in range(n_tests):
            x0 = np.random.random(ang_off+ang_n)*2 - 1
            print("testing magnitude gradients at ", x0, "(fs =", fs, ")")
            test_grads(wrap_emags, grad_wrap_emags, x0)
        print("\n====================================\n")

    @staticmethod
    def _guess_params(freqs, vfm, vfa, herm_n, ang_n, check_noise=True):
        ang_off = HERM_OFF + herm_n
        x0 = np.zeros(ang_off+ang_n)

        df = freqs[1] - freqs[0]
        #estimate the central frequency and envelope width in frequency space
        norm = 1/np.trapz(vfm[1:]/freqs[1:])
        x0[2] = np.trapz(vfm)*norm
        x0[3] = 2/np.sqrt(np.trapz(vfm*freqs)*norm - x0[2]**2)
        lo_fi, hi_fi = max(int((x0[2] - POLY_FIT_DEVS/x0[3])/df), 1), min(int((x0[2] + POLY_FIT_DEVS/x0[3])/df), len(freqs)-1)
        dd = x0[3]*(freqs-x0[2])
        #check whether noise was detected using the number of appreciable peaks (magnitude greater than a half of the largest)
        if check_noise:
            max_vfm = np.max(vfm)
            n_big_peaks = 0
            smallest_big = 0
            for i in range(1, len(vfm)-1):
                if vfm[i] > vfm[i-1] and vfm[i] > vfm[i+1]:
                    if vfm[i]/max_vfm > 0.5:
                        if n_big_peaks == 0:
                            smallest_big = freqs[i]
                        n_big_peaks += 1
            if n_big_peaks > NOISY_N_PEAKS:
                #apply a lowpass filter and try again
                vf = vfm*np.exp(1j*vfa)*np.sinc(x0[3]*(freqs-smallest_big))
                return signal._guess_params(freqs, np.abs(vf), np.angle(vf), herm_n, ang_n, check_noise=False)

        res = linregress(freqs[lo_fi:hi_fi], vfa[lo_fi:hi_fi])
        x0[0] = fix_angle(res.intercept)
        x0[1] = res.slope
        x0[HERM_OFF] = np.max(vfm)
        #return x0, lo_fi, hi_fi

        #fit the derivative of the phase to a polynomial of (f-f0)^2/sigma^2 so that only odd powers will appear when we integrate
        dd = freqs[lo_fi:hi_fi] - x0[2]
        vfa = fix_angle_seq(vfa[lo_fi:hi_fi], scan_n=1)
        vfm = vfm[lo_fi:hi_fi]
        def ang_cost(x):
            fit = x[0]*np.ones(dd.shape)
            for m in range(ang_n+1):
                fit += x[m+1]*dd**(2*m+1)
            return np.sum( (vfm*(vfa - fit))**2 )

        ang_opt_res = opt.minimize(ang_cost, np.zeros(2+ang_n))
        x0[1] = -ang_opt_res.x[1]*x0[3]
        x0[0] = ang_opt_res.x[0] + x0[1]*x0[2]
        for m in range(ang_n):
            x0[ang_off+m] = ang_opt_res.x[m+2]

        #now work out the Hermite-Gaussian expansion assuming w0 and alpha=1/sigma are held fixed at the values we just guessed
        #vfm = vfm*np.cos(vfa - ang_opt_res.x[0] - ang_opt_res.x[1]*dd)
        herm_den = 1
        k_ser = np.zeros(dd.shape)
        ks = np.zeros(2*herm_n)
        cs = np.zeros(2*herm_n)
        for m in range(2*herm_n):
            ks[m] = np.trapezoid(vfm*HERMS[m](dd)*np.exp(-dd**2/2), dx=df)/herm_den
            cs[m] = np.trapezoid(vfm*HERMS[m](dd)*np.exp(-dd**2/2)/freqs[lo_fi:hi_fi], dx=df)/herm_den
            k_ser += ks[m]*HERMS[m](dd)*np.exp(-dd**2/2)
            herm_den *= 2*(m+1)
        ks *= np.max(vfm)/np.max(k_ser)
        cs = np.dot(herm_conv(x0[2], x0[3], herm_n), ks)
        x0[HERM_OFF:HERM_OFF+herm_n] = cs[::2]
        #x0[HERM_OFF] = cs[0]

        #TODO: delete
        '''k_ser, c_ser = np.zeros(dd.shape), np.zeros(dd.shape)
        for i in range(herm_n):
            k_ser += (ks[2*i]*HERMS[2*i](dd) + ks[2*i+1]*HERMS[2*i+1](dd))*np.exp(-dd**2/2)
            c_ser += freqs[lo_fi:hi_fi]*cs[2*i]*HERMS[2*i](dd)*np.exp(-dd**2/2)
        plt.plot(dd, vfm, color='black')
        plt.plot(dd, vfm*np.sinc(dd))
        plt.plot(dd, k_ser*np.max(vfm)/np.max(k_ser))
        plt.plot(dd, c_ser*np.max(vfm)/np.max(c_ser))
        plt.show()'''

        return x0, lo_fi, hi_fi

    '''@staticmethod
    def _hess_mag(fs, x, dd, demags, cemags, herm_n, ang_n):
        ret = np.zeros((x.shape[0], x.shape[0], fs.shape[0]))
        ret[2,2,:] = -fs*de2mags*x[3]**2
        ret[3,3,:] = fs*de2mags*(fs-abs(x[2]))**2
        ret[2,3,:] = -fs*demags - dd*de2mags
        for m in range(herm_n):
            ret[HERM_OFF+m,2,:] = -fs*x[3]*(dd*HERMS[2*m](dd) - HERMS[2*m+1](dd))*np.exp(-dd**2/2)
            ret[2,HERM_OFF+m,:] = ret[HERM_OFF+m,2,:]
            ret[HERM_OFF+m,3,:] = fs*(fs-abs(x[2]))*(dd*HERMS[2*m](dd) - HERMS[2*m+1](dd))*np.exp(-dd**2/2)
            ret[3,HERM_OFF+m,:] = ret[HERM_OFF+m,3,:]
        for m in range(ang_n):
            ret[ang_off+m,2,:] = -np.sign(x[2])*x[3]*(2*m+3)*dd**(2*m+2)
            ret[ang_off+m,3,:] = (fs-abs(x[2]))*x[3]*(2*m+3)*dd**(2*m+2)
            ret[2,ang_off+m,:] = ret[ang_off+m,2,:]
            ret[3,ang_off+m,:] = ret[ang_off+m,3,:]
            ret[2,3,:] -= np.sign(x[2])*(2*m+3)**2*x[m+ang_off]*dd**(2*m+2)
            ret[3,2,:] = ret[2,3,:]
        ret[3,2,:] = ret[2,3,:]
        return ret'''

    def __init__(self, t_pts, v_pts, herm_n=3, ang_n=1, skip_opt=False, x0=None):
        self.herm_n = herm_n
        self.ang_n = ang_n
        ang_off = HERM_OFF+herm_n
        freqs = np.array([0,0])
        if t_pts is not None and v_pts is not None:
            self.dt = t_pts[1] - t_pts[0]
            #take the fourier transform       
            freqs = fft.rfftfreq(len(v_pts), d=self.dt)
            vf0 = fft.rfft(v_pts)
            self.vfm = np.abs(vf0)
            self.vfa = np.angle(vf0)
            res_norm = np.log(np.sum(self.vfm**2))
        #get guesses for initial parameters based on heuristics
        if x0 is not None:
            lo_fi, hi_fi = 0, 0
        else:
            x0, lo_fi, hi_fi = signal._guess_params(freqs, self.vfm, self.vfa, herm_n, ang_n)
        self.fit_bounds = np.array([lo_fi, hi_fi])
        #construct the cost function and analytic gradients
        def residuals(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],..., t_1=x[4]...
            emags, eargs, _, _ = signal._fourier_env_formx(freqs, x, herm_n=self.herm_n, ang_n=self.ang_n)
            emags *= freqs
            return np.sum(emags**2 + self.vfm**2 - 2*emags*self.vfm*np.cos(eargs-self.vfa))
            #return np.log( np.sum(emags**2 + self.vfm**2 - 2*emags*self.vfm*np.cos(eargs-self.vfa)) ) - res_norm
        def grad_res(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],..., t_1=x[4]...
            ret = np.zeros(x.shape)
            emags, eargs, grad_mag, grad_ang = signal._fourier_env_formx(freqs, x, herm_n=self.herm_n, ang_n=self.ang_n, calc_grads=True)
            emags *= freqs
            mag_as = np.reshape(np.tile(emags*self.vfm*np.sin(eargs-self.vfa), x.shape[0]), (x.shape[0], emags.shape[0]))
            mag_gs = np.reshape(np.tile(emags-self.vfm*np.cos(eargs-self.vfa), x.shape[0]), (x.shape[0], emags.shape[0]))
            return 2*np.sum(  mag_as*grad_ang + mag_gs*grad_mag, axis=1)
        def hess_res(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],..., t_1=x[4]...
            ret = np.zeros(x.shape)
            emags, eargs, dd, demags, cemags = signal._fourier_env_formx(freqs, x, herm_n=self.herm_n, ang_n=self.ang_n, calc_grads=True)
            emags *= freqs
            mag_as = np.reshape(np.tile(emags*self.vfm*np.sin(eargs-self.vfa), x.shape[0]), (x.shape[0], emags.shape[0]))
            mag_gs = np.reshape(np.tile(emags - self.vfm*np.cos(eargs-self.vfa), x.shape[0]), (x.shape[0], emags.shape[0]))
            grad_ang = signal._grad_ang(freqs, x, dd, self.herm_n, self.ang_n)
            grad_mag = signal._grad_mag(freqs, x, demags, cemags, self.herm_n, self.ang_n)
            return 2*np.sum(  mag_as*grad_ang + mag_gs*grad_mag, axis=1)
        if verbose > 4:
            test_grads(residuals, grad_res, x0)
        #finally we perform optimization if specified, otherwise just use the initial guess
        if skip_opt:
            xf = x0
        else:
            opt_res = opt.minimize(residuals, x0, jac=grad_res)
            xf = opt_res.x
            xf[2] = abs(xf[2])
            if verbose > 1:
                print("optimization message:", opt_res.message)
                print("x0:\t", x0)
                print("f(x0):\t", residuals(x0))
                print("x1:\t", xf)
                print("f(x1):\t", residuals(xf))
                if verbose > 4:
                    print("optimization result")
                    print(opt_res)
        #use the results
        xf[0] = fix_angle(xf[0])
        self.phi = xf[0]
        self.t0 = xf[1]/2/np.pi
        self.f0 = xf[2]
        self.sigma = 1/xf[3]
        self.herm_coeffs = xf[HERM_OFF:HERM_OFF+self.herm_n]
        self.poly_coeffs = xf[HERM_OFF+self.herm_n:HERM_OFF+self.herm_n+self.ang_n]
        self.x = xf
        if t_pts is not None and v_pts is not None:
            self.cost = residuals(xf)

    def get_fspace(self, freqs):
        x = np.array([self.phi, 2*np.pi*self.t0, self.f0, 1/self.sigma])
        x = np.append(x, self.herm_coeffs)
        x = np.append(x, self.poly_coeffs)
        emags, eargs, _, _ = signal._fourier_env_formx(freqs, x)
        return freqs*emags, eargs
    
    def __call__(self, ts):
        freqs = fft.rfftfreq(len(ts), d=ts[1]-ts[0])
        emags, eargs = self.get_fspace(freqs)
        return fft.irfft(emags*np.exp(1j*eargs))

    def get_fenv(self, freqs, field='E'):
        '''
        Get the frequency domain envelope for the electric field (or vector potential)
        freqs: the frequencies for which the field is computed
        field: if set to 'E' (default) compute the electric field, otherwise compute the envelope for the vector potential
        '''
        x = np.array([0, 0, 0, 1/self.sigma])
        x = np.append(x, self.herm_coeffs)
        x = np.append(x, self.poly_coeffs)
        emags, eargs, _, _ = signal._fourier_env_formx(freqs, x)
        if field == 'E':
            return 2*(freqs+self.f0)*emags*np.exp(1j*eargs)
        return 2*emags*np.exp(1j*eargs)

    def get_tenv(self, ts, field='E'):
        '''
        Get the time domain envelope for the electric field (or vector potential)
        freqs: the frequencies for which the field is computed
        field: if set to 'E' (default) compute the electric field, otherwise compute the envelope for the vector potential
        '''
        dt = ts[1]-ts[0]
        freqs = fft.fftfreq(len(ts), d=dt)
        fenv = self.get_fenv(freqs, field=field)*np.exp(-2j*np.pi*self.t0*freqs)
        return fft.ifft(fenv)
