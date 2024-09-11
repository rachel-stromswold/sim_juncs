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
EPSILON = 0.01
HERM_N = 3 #use the first HERM_N even hermite polynomials (i.e. HERM_N=1 uses just H_e0(w), HERM_N=2 adds H_e2(w) and so on) 
HERMS = [hermite(n) for n in range(2*HERM_N+2)]
HERM_SCALE = np.sqrt(2)
HERM_SCALE=1
HERM_OFF = 4
POLY_FIT_DEVS = 2 #the number of gaussian standard deviations to take when fitting the angles
ANG_POLY_OFF = HERM_OFF + HERM_N
ANG_POLY_N = 1 #in addition to phi and t0, consider higher order odd terms in the polynomial for the angle. i.e. arg(a(w)) = \phi - \omega t_0 + \sum_{n=1}^{ANG_POLY_N} \omega^{2n+1} t_n

MAX_T0_GUESSES = 2
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
verbose = 2

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
    def _fourier_env_formx(freqs, x, herm_n=HERM_N, ang_n=ANG_POLY_N, calc_grads=False):
        dd = x[3]*(freqs - abs(x[2]))
        eargs = x[0] - freqs*x[1]
        for m in range(ang_n):
            eargs += x[HERM_OFF+herm_n+m]*dd**(2*m+3)
        #eargs *= 2*np.pi
        emags = np.zeros(freqs.shape[0])
        demags = None
        cemags = None
        if calc_grads:
            demags = np.zeros(freqs.shape[0])
            cemags = np.zeros((herm_n, freqs.shape[0]))
            for m in range(HERM_N):
                emags += x[HERM_OFF+m]*HERMS[2*m](dd)
                demags += x[HERM_OFF+m]*(dd*HERMS[2*m](dd) - HERMS[2*m+1](dd))
                cemags[m,:] = HERMS[2*m](dd)*np.exp(-dd**2/2)
            demags *= np.exp(-dd**2/2)
        else:
            for m in range(HERM_N):
                emags += x[HERM_OFF+m]*HERMS[2*m](dd)
        emags *= np.exp(-dd**2/2)
        return emags, eargs, dd, demags, cemags

    @staticmethod
    def _do_grad_tests(fs, n_tests=10):
        print("\n====================================\n")

        def wrap_eargs(x):
            _, eargs, _, _, _ = signal._fourier_env_formx(fs, x)
            return eargs[0]
        def grad_wrap_eargs(x):
            _, eargs, dd, _, _ = signal._fourier_env_formx(fs, x, calc_grads=True)
            ret = np.zeros(x.shape)
            ret[0] = 1
            ret[1] = -fs[0]
            for m in range(ANG_POLY_N):
                ret[ANG_POLY_OFF+m] = np.sum( dd**(2*m+3) )
                ret[2] -= np.sign(x[2])*(2*m+3)*x[m+ANG_POLY_OFF]*np.sum( x[3]*dd**(2*m+2) )
                ret[3] += (2*m+3)*x[m+ANG_POLY_OFF]*np.sum( (fs-abs(x[2]))*dd**(2*m+2) )
            return ret
        for i in range(n_tests):
            x0 = np.random.random(ANG_POLY_OFF+ANG_POLY_N)*2 - 1
            print("testing argument gradients at ", x0, "(fs =", fs, ")")
            test_grads(wrap_eargs, grad_wrap_eargs, x0)
        print("\n====================================\n")
        def wrap_emags(x):
            emags, _, _, _, _ = signal._fourier_env_formx(fs, x)
            return fs[0]*emags[0]
        def grad_wrap_emags(x):
            emags, eargs, dd, demags, cemags = signal._fourier_env_formx(fs, x, calc_grads=True)
            ret = np.zeros(x.shape)
            ret[2] = -np.sign(x[2])*np.sum( fs*demags*x[3] )
            ret[3] = np.sum( fs*demags*(fs-abs(x[2])) )
            for m in range(HERM_N):
                ret[HERM_OFF+m] = np.sum( fs*cemags[m,:] )
            return ret
        for i in range(n_tests):
            x0 = np.random.random(ANG_POLY_OFF+ANG_POLY_N)*2 - 1
            print("testing magnitude gradients at ", x0, "(fs =", fs, ")")
            test_grads(wrap_emags, grad_wrap_emags, x0)
        print("\n====================================\n")

    @staticmethod
    def _guess_params(freqs, vfm, vfa):
        x0 = np.zeros(ANG_POLY_OFF+ANG_POLY_N)

        df = freqs[1] - freqs[0]
        #estimate the central frequency and envelope width in frequency space
        norm = 1/np.trapz(vfm[1:]/freqs[1:])
        x0[2] = np.trapz(vfm)*norm
        x0[3] = 2/np.sqrt(np.trapz(vfm*freqs)*norm)
        lo_fi, hi_fi = max(int((x0[2] - POLY_FIT_DEVS/x0[3])/df), 0), min(int((x0[2] + POLY_FIT_DEVS/x0[3])/df), len(freqs)-1)
        dd = x0[3]*(freqs-x0[2])

        #now work out the Hermite-Gaussian expansion assuming w0 and alpha=1/sigma are held fixed at the values we just guessed
        herm_den = 1
        k_ser = np.zeros(freqs.shape)
        ks = np.zeros(2*HERM_N)
        for m in range(2*HERM_N):
            ks[m] = np.trapezoid(vfm*HERMS[m](dd)*np.exp(-dd**2/2), dx=df)/herm_den
            k_ser += ks[m]*HERMS[m](dd)*np.exp(-dd**2/2)
            herm_den *= 2*(m+1)
        ks *= np.max(vfm)/np.max(k_ser)
        cs = np.dot(herm_conv(x0[2], x0[3], HERM_N), ks)
        x0[HERM_OFF:HERM_OFF+HERM_N] = cs[::2]
        #x0[HERM_OFF] = cs[0]

        #TODO: delete
        '''k_ser, c_ser = np.zeros(freqs.shape), np.zeros(freqs.shape)
        for i in range(HERM_N):
            k_ser += (ks[2*i]*HERMS[2*i](dd) + ks[2*i+1]*HERMS[2*i+1](dd))*np.exp(-dd**2/2)
            c_ser += freqs*cs[2*i]*HERMS[2*i](dd)*np.exp(-dd**2/2)
        plt.plot(dd, vfm)
        plt.plot(dd, k_ser*np.max(vfm)/np.max(k_ser))
        plt.plot(dd, c_ser*np.max(vfm)/np.max(c_ser))
        plt.show()'''

        #fit the derivative of the phase to a polynomial of (f-f0)^2/sigma^2 so that only odd powers will appear when we integrate
        dd = dd[lo_fi:hi_fi]
        vfa = fix_angle_seq(vfa[lo_fi:hi_fi], scan_n=1)
        vfm = vfm[lo_fi:hi_fi]
        def ang_cost(x):
            fit = x[0]*np.ones(dd.shape)
            for m in range(ANG_POLY_N+1):
                fit += x[m+1]*dd**(2*m+1)
            return np.sum( (vfm*(vfa - fit))**2 )

        ang_opt_res = opt.minimize(ang_cost, np.zeros(2+ANG_POLY_N))

        x0[1] = -ang_opt_res.x[1]*x0[3]
        x0[0] = ang_opt_res.x[0] + x0[1]*x0[2]
        for m in range(ANG_POLY_N):
            x0[ANG_POLY_OFF+m] = ang_opt_res.x[m+2]

        return x0, lo_fi, hi_fi

    def __init__(self, t_pts, v_pts, skip_opt=False):
        self.dt = t_pts[1] - t_pts[0]
        self.v_pts = v_pts
        #take the fourier transform       
        freqs = fft.rfftfreq(len(v_pts), d=self.dt)
        vf0 = fft.rfft(v_pts)
        self.vfm = np.abs(vf0)
        self.vfa = np.angle(vf0)
        #get guesses for initial parameters based on heuristics
        x0, lo_fi, hi_fi = signal._guess_params(freqs, self.vfm, self.vfa)
        self.fit_bounds = np.array([lo_fi, hi_fi])
        #construct the cost function and analytic gradients
        def residuals(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],..., t_1=x[4]...
            emags, eargs, _, _, _ = signal._fourier_env_formx(freqs, x)
            emags *= freqs
            #return np.sum(emags**2 + self.vfm**2 - 2*emags*self.vfm*np.cos(eargs-self.vfa))
            return np.log( np.sum(emags**2 + self.vfm**2 - 2*emags*self.vfm*np.cos(eargs-self.vfa)) )
        def grad_res(x):
            #phi=x[0], t_0=x[1], f_0=x[2], 1/sigma=x[3], c_0=x[4],..., t_1=x[4]...
            ret = np.zeros(ANG_POLY_OFF+ANG_POLY_N)
            emags, eargs, dd, demags, cemags = signal._fourier_env_formx(freqs, x, calc_grads=True)
            emags *= freqs
            demags *= freqs
            mag_as = 2*emags*self.vfm*np.sin(eargs-self.vfa)
            mag_gs = 2*( emags - self.vfm*np.cos(eargs-self.vfa) )
            ret[0] = np.sum( mag_as )
            ret[1] = -np.sum( freqs*mag_as )
            ret[2] = -np.sum( mag_gs*demags*x[3] )
            ret[3] = np.sum( mag_gs*demags*(freqs-x[2]) )
            for i in range(HERM_N):
                ret[i+HERM_OFF] = np.sum( freqs*mag_gs*cemags[i,:] )
            for m in range(ANG_POLY_N):
                ret[m+ANG_POLY_OFF] = np.sum( mag_as*dd**(2*m+3) )
                ret[2] -= (2*m+3)*x[m+ANG_POLY_OFF]*np.sum( mag_as*x[3]*dd**(2*m+2) )
                ret[3] += (2*m+3)*x[m+ANG_POLY_OFF]*np.sum( mag_as*(freqs-x[2])*dd**(2*m+2) )
            if x[2] < 0:
                ret[2] *= -1
            #return ret
            return ret/np.sum(emags**2 + self.vfm**2 - 2*emags*self.vfm*np.cos(eargs-self.vfa))
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
                print("diffs:\t", x0/xf)
                if verbose > 4:
                    print("optimization result")
                    print(opt_res)
        #use the results
        self.phi = xf[0]
        self.t0 = xf[1]/2/np.pi
        self.f0 = xf[2]
        self.sigma = 1/xf[3]
        self.herm_coeffs = xf[HERM_OFF:HERM_OFF+HERM_N]
        self.poly_coeffs = xf[ANG_POLY_OFF:ANG_POLY_OFF+ANG_POLY_N]

    def get_fspace(self, freqs):
        x = np.array([self.phi, 2*np.pi*self.t0, self.f0, 1/self.sigma])
        x = np.append(x, self.herm_coeffs)
        x = np.append(x, self.poly_coeffs)
        emags, eargs, _, _, _ = signal._fourier_env_formx(freqs, x)
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
        emags, eargs, _, _, _ = signal._fourier_env_formx(freqs, x)
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
