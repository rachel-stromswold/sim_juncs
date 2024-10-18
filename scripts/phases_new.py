from functools import partial
import time
import numpy as np
import scipy.optimize as opt
from scipy.stats import linregress
from scipy.special import hermite, lambertw
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
N_FIXED = 4
F0_I, ALPHA_I = 2, 3
HERM_OFF = 4
POLY_FIT_DEVS = 2 #the number of gaussian standard deviations to take when fitting the angles
ANG_POLY_N = 1 #in addition to phi and t0, consider higher order odd terms in the polynomial for the angle. i.e. arg(a(w)) = \phi - \omega t_0 + \sum_{n=1}^{ANG_POLY_N} \omega^{2n+1} t_n

NOISY_N_PEAKS = 4
NOISY_N_STD = 2

verbose = 3

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
Given a cost function cost_fn and an analytic gradient grad_fn, check that the analytic and numeric gradients approximately match
'''
def test_grads(cost_fn, grad_fn, x):
    nx = x.shape[0]
    if verbose < 1:
        return None
    success = True
    an_g = grad_fn(x)
    num_g = np.zeros(an_g.shape)
    for i in range(nx):
        ei = np.zeros(nx)
        ei[i] = EPSILON**2
        num_g[i] = ( cost_fn(x + ei) - cost_fn(x - ei) )/ei[i]/2
        if num_g[i]**2 > 1e-8 and np.abs( (num_g[i]-an_g[i])/num_g[i] ) > EPSILON:
            print("Warning! disagreement between gradients on axis", i, num_g[i], an_g[i])
            success = False
    if (verbose > 4): 
        print("gradients at x0:")
        print("\tnumeric:\t", num_g)
        print("\tanalytic:\t", an_g)
    return success

'''
Given a cost function cost_fn and an analytic gradient grad_fn, check that the analytic and numeric hessians approximately match
'''
def test_hess(cost_fn, hess_fn, x):
    if verbose < 1:
        return None
    success = True
    n = x.shape[0]
    an_h = hess_fn(x)
    num_h = np.zeros((n, n))
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = EPSILON**2
        for j in range(i, n):
            ej = np.zeros(n)
            ej[j] = EPSILON**2
            num_h[i,j] = 0.25*( cost_fn(x+ei+ej) - cost_fn(x+ei-ej) - cost_fn(x-ei+ej) + cost_fn(x-ei-ej) )/ei[i]/ej[j]
            num_h[j,i] = num_h[i,j]
            if num_h[i,j]**2 > 1e-8 and np.abs( (num_h[i,j]-an_h[i,j])/num_h[i,j] ) > 0.01:
                print("Warning! disagreement between hessians on axes", i, j, num_h[i,j], an_h[i,j])
                success = False
    if (verbose > 4): 
        print("hessians at x0:")
        print("\tnumeric:\t", num_h)
        print("\tanalytic:\t", an_h)
    return success

def hess_mul(hess, lin, n):
    return hess*np.reshape(np.tile(lin, n**2), (n, n, lin.shape[0]))
def grad_out(g1, g2, n):
    return np.tile(g1, n).reshape(n, n, g1.shape[1])*np.transpose(np.tile(g2, n).reshape(n, n, g2.shape[1]), axes=(1,0,2))

'''
Do simulated annealing using the callable log_like
c_fn: a callable for the log likelihood that accepts a numpy array with shape (x_reg.shape[0])
x0: a starting location for the sampling
n_steps: the number of steps to take
step_size: the size of step to take (a numpy array with the same shape as x0
'''
def anneal(c_fn, x0, n_steps, step_size, end_beta=100):
    c_old = c_fn(x0)
    n_accepted = 0
    step_size = np.abs(step_size)
    for i in range(n_steps):
        beta = end_beta*(i+1)/n_steps
        x1 = x0 + np.random.normal(0,1)*step_size
        c_new, g_new = c_fn(x1)
        x1 -= 0.1*g_new*step_size
        u = np.log(np.random.random())
        if u < beta*(c_old - c_new):
            n_accepted += 1
            if verbose > 4:
                print("accepted move c = {} -> {}".format(c_old, c_new))
            x0 = x1
            c_old = c_new
    if verbose > 2:
        print("accept probability: ", n_accepted/n_steps)
    return x0

class signal:
    @staticmethod
    def _fourier_env_formx(x, freqs, herm_n, ang_n, calc_grads=False):
        nf = freqs.shape[0]
        fs = np.tile(freqs, 2).reshape((2,nf))
        sign_dw = np.meshgrid(np.ones(nf), np.array([-1,1]))[1]
        #calculate arguments
        n = x.shape[0]

        ang_off = HERM_OFF+herm_n
        dw = fs[:,:] + sign_dw*x[2]
        eargs = -sign_dw*x[0] - fs*x[1]
        for m in range(ang_n):
            eargs += x[ang_off+m]*dw**(2*m+3)
        #eargs *= 2*np.pi
        emags = np.zeros(fs.shape)
        if calc_grads:
            #calculate angle gradients
            ang_grads = np.zeros((2, n, freqs.shape[0]))
            ang_grads[0,0,:] = 1
            ang_grads[1,0,:] = -1
            ang_grads[:,1,:] = -fs
            ang_hess = np.zeros((2, n, n, nf))
            for m in range(ang_n):
                ang_grads[:,ang_off+m,:] = dw**(2*m+3)
                ang_grads[:,2,:] += sign_dw*(2*m+3)*x[m+ang_off]*dw**(2*m+2)
                ang_hess[:,ang_off+m,2,:] = sign_dw*(2*m+3)*dw**(2*m+2)
                ang_hess[:,2,ang_off+m,:] = ang_hess[:,ang_off+m,2,:]
                ang_hess[:,2,2,:] += (2*m+3)*(2*m+2)*x[m+ang_off]*dw**(2*m+1)
            #calculate magnitude gradients
            dw *= x[3]
            mag_grads = np.zeros((2, n, nf))
            mag_hess = np.zeros((2, n, n, nf))
            for m in range(herm_n):
                emags += x[HERM_OFF+m]*HERMS[2*m](dw)
                #gradients
                demag = sign_dw*(dw*HERMS[2*m](dw) - HERMS[2*m+1](dw))*np.exp(-dw**2/2)
                de2mag = sign_dw*( (1+dw**2)*HERMS[2*m](dw) - 2*dw*HERMS[2*m+1](dw) + HERMS[2*m+2](dw))*np.exp(-dw**2/2)
                mag_grads[:,2,:] += sign_dw*x[3]*x[HERM_OFF+m]*demag
                mag_grads[:,3,:] += dw*x[HERM_OFF+m]*demag/x[3]
                mag_grads[:,HERM_OFF+m,:] = sign_dw*HERMS[2*m](dw)*np.exp(-dw**2/2)
                #hessian
                mag_hess[:,2,2,:] += x[3]**2*x[HERM_OFF+m]*de2mag
                mag_hess[:,3,3,:] += dw**2*x[HERM_OFF+m]*de2mag/x[3]**2
                mag_hess[:,2,3,:] += sign_dw*(x[HERM_OFF+m]*demag + dw*x[HERM_OFF+m]*de2mag)
                mag_hess[:,3,2,:] = mag_hess[:,2,3,:]
                mag_hess[:,HERM_OFF+m,2,:] = sign_dw*x[3]*demag
                mag_hess[:,2,HERM_OFF+m,:] = mag_hess[:,HERM_OFF+m,2,:]
                mag_hess[:,HERM_OFF+m,3,:] += dw*demag/x[3]
                mag_hess[:,3,HERM_OFF+m,:] = mag_hess[:,HERM_OFF+m,3,:]
            emags *= sign_dw*np.exp(-dw**2/2)
            return emags, eargs, mag_grads, ang_grads, mag_hess, ang_hess
        else:
            dw *= x[3]
            for m in range(herm_n):
                emags += x[HERM_OFF+m]*HERMS[2*m](dw)
            emags *= sign_dw*np.exp(-dw**2/2)
        return emags, eargs

    #construct the cost function and analytic gradients
    def residuals(self, x, freqs):
        nx = 4 + self.herm_n + self.ang_n
        nf = freqs.shape[0]
        #np.tile(a, (2,1)).reshape((2,2,10))
        #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],..., t_1=x[4]...
        grad = np.zeros(x.shape)
        emags, eargs, grad_mag, grad_ang, _, _ = signal._fourier_env_formx(x, freqs, self.herm_n, self.ang_n, calc_grads=True)
        emags *= freqs
        grad_mag *= freqs
        if self.discard_negative:
            mag_gs = np.tile(emags[0,:]-self.vfm*np.cos(eargs[0,:]-self.vfa), nx).reshape(nx, freqs.shape[0])
            mag_as = np.tile(emags[0,:]*self.vfm*np.sin(eargs[0,:]-self.vfa), nx).reshape(nx, freqs.shape[0])
            cc = np.sum( self.vfm**2 + emags[0,:]**2 - 2*emags[0,:]*self.vfm*np.cos(eargs[0,:]-self.vfa) )
            gcc = 2*np.sum(mag_as*grad_ang[0,:,:] + mag_gs*grad_mag[0,:,:], axis=1)
        else:
            mag_gs = emags - self.vfm2*np.cos(eargs-self.vfa2) + emags[::-1,:]*np.cos(eargs[1,:]-eargs[0,:])
            mag_as = emags*self.vfm2*np.sin(eargs-self.vfa2) - self.sign_dw*np.prod(emags,axis=0)*np.sin(np.sum(self.sign_dw*eargs,axis=0))
            mag_as = np.repeat(mag_as, nx, axis=0).reshape(2, nx, nf)
            mag_gs = np.repeat(mag_gs, nx, axis=0).reshape(2, nx, nf)
            cc = np.sum( self.vfm**2
                        + np.sum(emags**2 - 2*emags*self.vfm2*np.cos(eargs-self.vfa2),axis=0)
                        + 2*np.prod(emags,axis=0)*np.cos(np.sum(self.sign_dw*eargs,axis=0))
                        )
            gcc = 2*np.sum(np.sum(mag_as*grad_ang + mag_gs*grad_mag, axis=0), axis=1)
        #return cc, gcc
        lp = self.log_prior(x)
        return np.log(cc) - lp[0], gcc/cc - lp[1]

    def hess_res(self, x, freqs):
        nx = 4 + self.herm_n + self.ang_n
        nf = freqs.shape[0]
        #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],..., t_1=x[4]...
        emags, eargs, grad_mag, grad_ang, hess_mag, hess_ang = signal._fourier_env_formx(x, freqs, self.herm_n, self.ang_n, calc_grads=True)
        emags *= freqs
        grad_mag *= freqs
        hess_mag *= freqs
        ret = np.zeros((nx, nx, nf))
        cc = self.res_norm
        for i in range(2 - self.discard_negative):
            phase_prod = np.tile(self.vfm*np.exp(1j*(eargs[i,:]-self.vfa)), nx**2).reshape((nx,nx,nf))
            ret += hess_mul(hess_mag[i,:,:], emags[i,:], nx)
            ret += grad_out(grad_mag[i,:,:], grad_mag[i,:], nx)
            ret += np.imag(phase_prod)*(grad_out(grad_mag[i,:,:],grad_ang[i,:,:],nx)
                                        +grad_out(grad_ang[i,:,:],grad_mag[i,:,:],nx)
                                        +hess_mul(hess_ang[i,:,:,:], emags[i,:],nx))
            ret += np.real(phase_prod)*(hess_mul(grad_out(grad_ang[i,:,:],grad_ang[i,:,:],nx),emags[i,:],nx)-hess_mag[i,:,:,:])
            cc += np.sum(emags[i,:]**2 - 2*emags[i,:]*self.vfm*np.cos(eargs[i,:]-self.vfa))
        if not self.discard_negative:
            c10, s10 = np.cos(eargs[1,:]-eargs[0,:]), np.sin(eargs[1,:]-eargs[0,:])
            cc += np.sum(2*np.prod(emags, axis=0)*np.cos(eargs[1,:]-eargs[0,:]))
            ret += hess_mul(hess_mag[0,:,:,:], emags[1,:]*c10, nx)
            ret += hess_mul(hess_mag[1,:,:,:], emags[0,:]*c10, nx)
            ret += hess_mul(grad_out(grad_mag[0,:,:], grad_mag[1,:,:], nx), c10, nx)
            ret += hess_mul(grad_out(grad_mag[1,:,:], grad_mag[0,:,:], nx), c10, nx)
            ret += hess_mul(grad_out(grad_ang[0,:,:], np.sum(self.sign_2dw*grad_ang, axis=0), nx), np.prod(emags,axis=0)*c10, nx)
            ret -= hess_mul(grad_out(grad_ang[1,:,:], np.sum(self.sign_2dw*grad_ang, axis=0), nx), np.prod(emags,axis=0)*c10, nx)
            ret -= hess_mul(grad_out(grad_mag[0,:,:], np.sum(self.sign_2dw*grad_ang, axis=0), nx), emags[1,:]*s10, nx)
            ret -= hess_mul(grad_out(grad_mag[1,:,:], np.sum(self.sign_2dw*grad_ang, axis=0), nx), emags[0,:]*s10, nx)
            ret -= hess_mul(grad_out(np.sum(self.sign_2dw*grad_ang, axis=0), grad_mag[0,:,:], nx), emags[1,:]*s10, nx)
            ret -= hess_mul(grad_out(np.sum(self.sign_2dw*grad_ang, axis=0), grad_mag[1,:,:], nx), emags[0,:]*s10, nx)
            ret -= hess_mul(hess_ang[1,:,:,:], emags[0,:]*emags[1,:]*s10, nx)
            ret += hess_mul(hess_ang[0,:,:,:], emags[0,:]*emags[1,:]*s10, nx)
        lp = self.log_prior(x)
        gr = self.residuals(x, freqs)[1] + lp[1]
        #return 2*np.sum(ret, axis=2)
        return 2*np.sum(ret, axis=2)/cc - np.outer(gr,gr) - lp[2]

    @staticmethod
    def _do_grad_tests(fs, herm_n, ang_n, n_tests=10):
        print("\n====================================\n")
        ang_off = HERM_OFF+herm_n
        def wrap_eargs(x, dim=0):
            _, eargs = signal._fourier_env_formx(x, fs, herm_n, ang_n)
            return np.sum(eargs[dim,:])
        def grad_wrap_eargs(x, dim=0):
            _, _, _, ang_grad, _, _ = signal._fourier_env_formx(x, fs, herm_n, ang_n, calc_grads=True)
            return np.sum(ang_grad[dim,:,:], axis=1)
        def hess_wrap_eargs(x, dim=0):
            _, _, _, _, _, ang_hess = signal._fourier_env_formx(x, fs, herm_n, ang_n, calc_grads=True)
            return np.sum(ang_hess[dim,:,:,:], axis=2)
        for i in range(n_tests):
            x0 = np.random.random(ang_off+ang_n)*2 - 1
            print("testing argument gradients at ", x0, "(fs =", fs, ")")
            success = test_grads(wrap_eargs, grad_wrap_eargs, x0)
            success *= test_hess(wrap_eargs, hess_wrap_eargs, x0)
            success *= test_grads(partial(wrap_eargs,dim=1), partial(grad_wrap_eargs,dim=1), x0)
            success *= test_hess(partial(wrap_eargs,dim=1), partial(hess_wrap_eargs,dim=1), x0)
            if success:
                print("\033[92mSuccess!\033[0m")
            else:
                print("\033[31mFailure!\033[0m")
        print("\n====================================\n")
        def wrap_emags(x, dim=0):
            emags, _ = signal._fourier_env_formx(x, fs, herm_n, ang_n)
            return np.sum(emags[dim,:])
        def grad_wrap_emags(x, dim=0):
            _, _, mag_grad, _, _, _ = signal._fourier_env_formx(x, fs, herm_n, ang_n, calc_grads=True)
            return np.sum(mag_grad[dim,:,:], axis=1)
        def hess_wrap_emags(x, dim=0):
            _, _, _, _, mag_hess, _ = signal._fourier_env_formx(x, fs, herm_n, ang_n, calc_grads=True)
            return np.sum(mag_hess[dim,:,:,:], axis=2)
        for i in range(n_tests):
            x0 = np.random.random(ang_off+ang_n)*2 - 1
            print("testing magnitude gradients at ", x0, "(fs =", fs, ")")
            success = test_grads(wrap_emags, grad_wrap_emags, x0)
            success *= test_hess(wrap_emags, hess_wrap_emags, x0)
            success *= test_grads(partial(wrap_emags,dim=1), partial(grad_wrap_emags,dim=1), x0)
            success *= test_hess(partial(wrap_emags,dim=1), partial(hess_wrap_emags,dim=1), x0)
            if success:
                print("\033[92mSuccess!\033[0m")
            else:
                print("\033[31mFailure!\033[0m")
        print("\n====================================\n")

    def _guess_params_opt(self, freqs, check_noise=False, rel_height=0.1):
        df = freqs[1] - freqs[0]
        div_sig = np.append(np.zeros(1), self.vfm[1:]/freqs[1:])
        mag_x = np.zeros(self.herm_n+2)
        norm = 1/np.trapz(div_sig)
        mag_x[0] = np.trapz(div_sig*freqs)*norm
        mag_x[1] = 1/np.sqrt(np.trapz(div_sig*freqs**2)*norm - mag_x[0]**2)
        mag_x[2] = -np.max(div_sig)

        def mag_res(x):
            emags, _, mag_grads, _, _, _ = signal._fourier_env_formx(np.append(np.zeros(2), x), freqs[1:], self.herm_n, 0, calc_grads=True)
            cc = np.sum( (emags[0,:]-div_sig[1:])**2 )
            gcc = mag_grads[0,2:]*np.tile(emags[0,:]-div_sig[1:], self.herm_n+2).reshape(self.herm_n+2, freqs.shape[0]-1)
            lp = self.log_prior(np.append(np.zeros(2), x))
            return np.log(cc) - lp[0], 2*np.sum(gcc, axis=1)/cc - lp[1][2:HERM_OFF+self.herm_n]

        #test_grads(lambda x: mag_res(x)[0], lambda x: mag_res(x)[1], mag_x)
        t_start = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
        opt_res = opt.minimize(mag_res, mag_x, jac=True)
        if verbose > 1:
            t_ns = time.clock_gettime_ns(time.CLOCK_MONOTONIC)-t_start
            print("initial optimization {} in {} ms".format("succeeded" if opt_res.success else "failed", t_ns*1e-6))
            if verbose > 3:
                print(opt_res)

        if check_noise:
            peaks, props = ssig.find_peaks(div_sig, height=-opt_res.x[2]*rel_height)
            #if len(peaks) > NOISY_N_PEAKS:
            if len(peaks) > 0 and freqs[peaks[-1]] - opt_res.x[0] > NOISY_N_STD/abs(opt_res.x[1]):
                if verbose > 2:
                    print("noisy signal detected applying bandpass filter around f={}, w={}".format(opt_res.x[0], opt_res.x[1]))
                old_vfm = self.vfm[:]
                self.vfm *= np.sinc(opt_res.x[1]*(freqs-opt_res.x[0]))
                x0, lo_fi, hi_fi = self._guess_params_opt(freqs, check_noise=False)
                self.vfm = old_vfm
                return x0, lo_fi, hi_fi

        #TODO:delete
        '''emags_before, _ = signal._fourier_env_formx(np.append(np.zeros(2), mag_x), freqs, self.herm_n, 0, calc_grads=False)
        emags_after, _  = signal._fourier_env_formx(np.append(np.zeros(2), opt_res.x), freqs, self.herm_n, 0, calc_grads=False)
        plt.plot(freqs, div_sig)
        plt.plot(freqs, emags_before[0,:])
        plt.plot(freqs, emags_after[0,:])
        plt.show()'''

        x0 = np.zeros(HERM_OFF+self.herm_n+self.ang_n)
        x0[2:HERM_OFF+self.herm_n] = opt_res.x
        x0[2] = abs(x0[2])
        lo_fi, hi_fi = max(int((x0[2] - POLY_FIT_DEVS/x0[3])/df), 1), min(int((x0[2] + POLY_FIT_DEVS/x0[3])/df), len(freqs)-1)
        #now guess angles
        if hi_fi <= lo_fi:
            return x0, 1, 2
        res = linregress(freqs[lo_fi:hi_fi], fix_angle_seq(self.vfa[lo_fi:hi_fi]))
        x0[0] = fix_angle(res.intercept)
        x0[1] = -res.slope
        return x0, lo_fi, hi_fi

    @staticmethod
    def _guess_params(freqs, vfm, vfa, herm_n, ang_n, check_noise=True, rel_height=0.1, f0_hint=-1):
        df = freqs[1] - freqs[0]
        ang_off = HERM_OFF + herm_n
        x0 = np.zeros(ang_off+ang_n)
        #estimate the central frequency and envelope width in frequency space
        div_sig = np.append(np.zeros(1), vfm[1:]/freqs[1:])
        norm = 1/np.trapz(div_sig)
        x0[2] = np.trapz(div_sig*freqs)*norm
        x0[3] = 1/np.sqrt(np.trapz(div_sig*freqs**2)*norm - x0[2]**2) 
        x0[HERM_OFF] = -np.max(div_sig)
        lo_fi, hi_fi = max(int((x0[2] - POLY_FIT_DEVS/x0[3])/df), 1), min(int((x0[2] + POLY_FIT_DEVS/x0[3])/df), len(freqs)-1)
        peaks, props = ssig.find_peaks(div_sig, height=np.max(div_sig)*rel_height)
        if check_noise and len(peaks) > NOISY_N_PEAKS:
            if f0_hint <= 0:
                f0_hint = freqs[peaks[0]]
            if verbose > 2:
                print("noisy signal detected applying bandpass filter around f={}, w={}".format(f0_hint, x0[3]))
            vf = vfm*np.exp(1j*vfa)*np.sinc(x0[3]*(freqs-f0_hint))
            return signal._guess_params(freqs, np.abs(vf), np.angle(vf), herm_n, ang_n, check_noise=False)
        '''elif len(peaks) > 1:
            #sort peaks in descending order of height
            x0[2] = np.sum(freqs[peaks]*props["peak_heights"])/np.sum(props["peak_heights"])
            dd = x0[3]*(freqs-x0[2])
            herm_den = np.sqrt(np.pi)
            for m in range(herm_n):
                x0[HERM_OFF+m] = np.trapz(div_sig*HERMS[2*m](dd)*np.exp(-dd**2/2), dx=df*x0[3])/herm_den
                herm_den *= 4*(m+1)*(m+2)'''

        #guess the phase and central time
        res = linregress(freqs[lo_fi:hi_fi], fix_angle_seq(vfa[lo_fi:hi_fi]))
        x0[0] = fix_angle(res.intercept)
        x0[1] = -res.slope

        return x0, lo_fi, hi_fi

    @staticmethod
    def _guess_params_old(freqs, vfm, vfa, herm_n, ang_n, check_noise=True):
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
                return signal._guess_params_old(freqs, np.abs(vf), np.angle(vf), herm_n, ang_n, check_noise=False)

        res = linregress(freqs[lo_fi:hi_fi], vfa[lo_fi:hi_fi])
        x0[0] = fix_angle(res.intercept)
        x0[1] = res.slope
        x0[HERM_OFF] = -np.max(vfm)
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
        x0[1] = -ang_opt_res.x[1]
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
            ks[m] = np.trapz(vfm*HERMS[m](dd)*np.exp(-dd**2/2), dx=df)/herm_den
            cs[m] = np.trapz(vfm*HERMS[m](dd)*np.exp(-dd**2/2)/freqs[lo_fi:hi_fi], dx=df)/herm_den
            k_ser += ks[m]*HERMS[m](dd)*np.exp(-dd**2/2)
            herm_den *= 2*(m+1)
        ks *= np.max(vfm)/np.max(k_ser)
        cs = np.dot(herm_conv(x0[2], x0[3], herm_n), ks)
        x0[HERM_OFF:HERM_OFF+herm_n] = -cs[::2]
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

    def __init__(self, t_pts, v_pts, herm_n=3, ang_n=1, skip_opt=False, x0=None, method="trust-exact", discard_negative=False, run_grad_tests=False, log_prior=None):
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
            self.res_norm = np.sum(self.vfm**2)
            self.vfm2 = np.tile(self.vfm, 2).reshape(2, freqs.shape[0])
            self.vfa2 = np.tile(self.vfa, 2).reshape(2, freqs.shape[0])
        nx,nf = HERM_OFF+herm_n+ang_n, freqs.shape[0]
        self.discard_negative = discard_negative
        self.ts_include = 2
        self.sign_dw = np.meshgrid(np.ones(nf), np.array([-1,1]))[1]
        self.sign_2dw = np.repeat(self.sign_dw, nx, axis=0).reshape(2,nx,nf)
        if discard_negative:
            ts_include = 1
        if log_prior is None:
            self.log_prior = lambda x: (0, np.zeros(nx), np.zeros((nx, nx)))
        else:
            self.log_prior = log_prior

        #get guesses for initial parameters based on heuristics
        if x0 is not None:
            lo_fi, hi_fi = 0, 0
        else:
            x0, lo_fi, hi_fi = self._guess_params_opt(freqs)
        x0[2] = abs(x0[2])

        self.fit_bounds = np.array([lo_fi, hi_fi])

        #finally we perform optimization if specified, otherwise just use the initial guess
        if skip_opt:
            xf = x0
        else:
            this_args = (freqs,)
            try:
                opt_res = opt.minimize(self.residuals, x0, jac=True, hess=self.hess_res, args=this_args, method=method)

                if not opt_res.success:
                    opt_res = opt.minimize(self.residuals, x0, jac=True, args=this_args)
                    print("used BFGS because of non-convergence")
                else:
                    print("used trust-exact")
                '''ssize = np.sqrt( np.abs(0.5*residuals(x0)/np.diag(hess_res(x0))) )
                opt_res = opt.basinhopping(residuals, x0, minimizer_kwargs={"method":method, "jac":True, hess:hess_res})'''
            except:
                #use bfgs as a fallback if trust-exact crashes
                opt_res = opt.minimize(self.residuals, x0, jac=True, args=this_args)
                print("used BFGS because of error")
            if not opt_res.success:
                print("\033[31mOptimization failed!\033[0m")
            xf = opt_res.x
            xf[2] = abs(xf[2])
            if verbose > 2:
                print("optimization message:", opt_res.message)
                print("x0:\t", x0)
                print("f(x0):\t", self.residuals(x0, freqs))
                print("x1:\t", xf)
                print("f(x1):\t", self.residuals(xf, freqs))
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
            self.cost = self.residuals(xf, freqs)[0] - np.log(self.res_norm)
            self.residual = self.cost + self.log_prior(xf)[0] #the residual should not include information from the prior
            if verbose > 0 and not skip_opt:
                print( "R^2 = ", 1-np.exp(self.residual), "n_it = ", opt_res.nit )
        if verbose > 1 and not skip_opt:
            print("nfev =", opt_res.nfev, "njev =", opt_res.njev, "nhev =", opt_res.njev)

    def get_fspace(self, freqs):
        #TODO: include negative frequency terms
        x = np.array([self.phi, 2*np.pi*self.t0, self.f0, 1/self.sigma])
        x = np.append(x, self.herm_coeffs)
        x = np.append(x, self.poly_coeffs)
        emags, eargs = signal._fourier_env_formx(x, freqs, self.herm_n, self.ang_n)
        if self.discard_negative:
            return freqs*emags[0,:]*np.exp(1j*eargs[0,:])
        return freqs*np.sum(emags*np.exp(1j*eargs), axis=0)
    
    def __call__(self, ts):
        freqs = fft.rfftfreq(len(ts), d=ts[1]-ts[0])
        return fft.irfft( self.get_fspace(freqs) )

    def get_fenv(self, freqs, field='E'):
        '''
        Get the frequency domain envelope for the electric field (or vector potential)
        freqs: the frequencies for which the field is computed
        field: if set to 'E' (default) compute the electric field, otherwise compute the envelope for the vector potential
        '''
        x = np.array([0, 0, 0, 1/self.sigma])
        x = np.append(x, self.herm_coeffs)
        x = np.append(x, self.poly_coeffs)
        emags, eargs = signal._fourier_env_formx(x, freqs, self.herm_n, self.ang_n)
        emags = emags[0,:]
        eargs = eargs[0,:]
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

    '''
    Extract the effective phase and central time by redefining t0 to be the average of the time envelope
    '''
    def get_eff_t0_phi(self, ts):
        tenv = np.abs(self.get_tenv(ts))
        t0 = np.trapz(tenv*ts)/np.trapz(tenv)
        phi = 2*np.pi*self.f0*(t0 - self.t0) + self.phi
        return t0, phi
