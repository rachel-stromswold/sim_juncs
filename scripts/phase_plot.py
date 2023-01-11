import phases
import utils
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
import scipy.fft
import time
import pickle
import os.path

plt.rc('font', size=14)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

N_COLS = 1

AMP_RANGE = (0, 1.1)
OMG_RANGE = (0, 4.0)
SIG_RANGE = (0, 10.0)
PHI_RANGE = (-1.0, 1.0)
PAD_FACTOR = 1.05
#FREQ_0 = 6.56
FREQ_0 = 0.43
N_FREQ_COMPS=100
EPS_0 = 200.6
#EPS_0 = 2.5

#The highest waveguide mode (n) that is considered in fits sin((2n-1)pi/L)
HIGHEST_MODE=3

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Fit time a_pts data to a gaussian pulse envelope to perform CEP estimation.')
parser.add_argument('--fname', type=str, help='h5 file to read', default='field_samples.h5')
parser.add_argument('--gap-width', type=float, help='junction width', default=0.1)
parser.add_argument('--gap-thick', type=float, help='junction thickness', default=0.2)
parser.add_argument('--diel-const', type=float, help='dielectric constant of material', default=3.5)
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
parser.add_argument('--slice-dir', type=str, help='prefix to use when opening files', default='x')
parser.add_argument('--regression', action='store_true', help='If set to true, perform linear regression on phases across junction', default=False)
args = parser.parse_args()

#decide whether slices are along the x or z axis
slice_dir = args.slice_dir
slice_name = 'z'
if slice_dir == 'x':
    slice_name = 'z'

class phase_finder:
    '''
    Initialize a phase finder reading from the  specified by fname and a juction width and gap specified by width and height respectively.
    fname: h5 file to read time samples from
    width: the width of the junction
    height: the thickness of the junction
    pass_alpha: The scale of the low pass filter applied to the time series. This corresponds to taking a time average of duration 2/pass_alpha on either side
    '''
    def __init__(self, fname, width, height, pass_alpha=0.5):
        '''read metadata from the h5 file'''
        self.geom = utils.Geometry("params.conf", gap_width=width, gap_thick=height)
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
        #figure out the extent of the simulation points in the spanned direction
        z_min = min(self.geom.meep_len_to_um(self.geom.l_junc -self.geom.z_center), \
                self.geom.meep_len_to_um(self.f[self.clust_names[0]]['locations'][slice_dir][0]-self.geom.z_center))
        self.x_range = (PAD_FACTOR*z_min, -PAD_FACTOR*z_min)
        #read information to figure out time units and step sizes
        t_min = self.f['info']['time_bounds'][0]
        t_max = self.f['info']['time_bounds'][1]
        self.dt = (t_max-t_min)/self.n_t_pts
        self.t_pts = np.linspace(t_min, t_max, num=self.n_t_pts)
        #these are used for applying a low pass filter to the time data, the low-pass is a sinc function applied in frequency space
        self.f_pts = scipy.fft.fftfreq(self.n_t_pts, d=self.dt)
        self.low_filter = np.sin(self.f_pts*np.pi*pass_alpha) / (self.f_pts*np.pi*pass_alpha)
        self.low_filter[0] = 1
        #this normalizes the square errors to have reasonable units
        self.sq_er_fact = self.dt/(t_max-t_min)
        #figure out a range of incident frequencies
        avg_field_0 = self.f['info']['sources']['amplitude']
        avg_sig_0 = self.f['info']['sources']['width']
        avg_phase_0 = self.f['info']['sources']['phase']
        self.df = 6/(avg_sig_0*N_FREQ_COMPS)
        self.freq_range = np.linspace(-FREQ_0-3/avg_sig_0, -FREQ_0+3/avg_sig_0, num=N_FREQ_COMPS)
        self.freq_comps = avg_field_0*avg_sig_0*np.exp(-1j*avg_phase_0-((self.freq_range+FREQ_0)*avg_sig_0)**2/2)/np.sqrt(2*np.pi)

    def get_wavelen(self):
        return self.f['info']['sources']['wavelen']
    def get_width(self):
        return self.f['info']['sources']['width']
    def get_phase(self):
        return self.f['info']['sources']['phase']
    def get_start_time(self):
        return self.f['info']['sources']['start_time']
    def get_end_time(self):
        return self.f['info']['sources']['end_time']
    def get_amplitude(self):
        return self.f['info']['sources']['amplitude']

    def get_clust_location(self, clust):
        '''return the location of the cluster with the name clust along the propagation direction (z if slice_dir=x x if slice_dir=z)'''
        return self.geom.meep_len_to_um(self.f[clust]['locations'][slice_name][0] - self.geom.t_junc)

    def get_point_times(self, clust, ind, low_pass=False):
        #fetch a list of points and their associated coordinates
        points = list(self.f[clust].keys())[1:]
        v_pts = np.array(self.f[clust][points[ind]]['time']['Re']) + 1j*np.array(self.f[clust][points[ind]]['time']['Im'])
        # (dt*|dE/dt|)^2 evaluated at t=t_max (one factor of two comes from the imaginary part, the other from the gaussian. We do something incredibly hacky to account for the spatial part. We assume that it has an error identical to the time error, and add the two in quadrature hence the other factor of two.
        #err_2 = 8*np.max(np.diff(v_pts)**2)
        err_2 = 0.02
        if low_pass:
            vf_pts = scipy.fft.fft(v_pts)*self.low_filter
            v_pts = scipy.fft.ifft(vf_pts)
        return v_pts, err_2

    def read_cluster(self, clust):
        '''Read an h5 file and return three two dimensional numpy arrays of the form ([amplitudes, errors], [sigmas, errors], [phases, errors]
        '''
        #for performance metrics
        n_evals = 0
        t_start = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
        x = self.get_clust_location(clust)
        #fetch a list of points and their associated coordinates
        points = list(self.f[clust].keys())[1:]
        zs = (np.array(self.f[clust]['locations'][slice_dir])-self.geom.z_center)/self.geom.um_scale
        n_z_pts = len(points)
        if zs.shape[0] != n_z_pts:
            raise ValueError("points and locations do not have the same size")

        amp_arr = np.zeros((2,n_z_pts))
        t_0_arr = np.zeros((2,n_z_pts))
        sig_arr = np.zeros((2,n_z_pts))
        omega_arr = np.zeros((2,n_z_pts))
        phase_arr = np.zeros((2,n_z_pts))
        good_zs = []
        phase_cor = 0
        for j in range(len(points)):
            v_pts, err_2 = self.get_point_times(clust, j, low_pass=False)
            #before doing anything else, save a plot of just the time series
            plt.plot(self.t_pts, v_pts)
            plt.savefig("{}/fit_figs/t_series_{}_{}.pdf".format(args.prefix,clust,j))
            plt.clf()
            print("clust={}, j={}\n\tfield error^2={}".format(clust, j, err_2))
            try:
                res,res_env = phases.opt_pulse_env(self.t_pts, v_pts, a_sigmas_sq=err_2, keep_n=2.5, fig_name="{}/fit_figs/fit_{}_{}".format(args.prefix,clust,j))
            except:
                print("fitting to raw time series failed, applying low pass filter")
                v_pts, err_2 = self.get_point_times(clust, j, low_pass=True)
                res,res_env = phases.opt_pulse_env(self.t_pts, v_pts, a_sigmas_sq=err_2, keep_n=2.5, fig_name="{}/fit_figs/fit_{}_{}".format(args.prefix,clust,j))
            n_evals += 1
            print("\tsquare errors = {}, {}\n\tx={}\n\tdiag(H^-1)={}".format(res.fun, res_env.fun, res.x, np.diagonal(res.hess_inv)))
            #only include this point if the fit was good
            if res.fun*self.sq_er_fact < 50.0:
                good_zs.append(zs[j])
                jj = len(good_zs)-1
                amp, t_0, sig, omega, cep = phases.get_params(res)
                amp_arr[0,jj] = amp
                amp_arr[1,jj] = np.sqrt(res_env.hess_inv[0][0]/(err_2*self.n_t_pts))/2
                t_0_arr[0,jj] = t_0
                t_0_arr[1,jj] = np.sqrt(res.hess_inv[1][1]/(err_2*self.n_t_pts))
                sig_arr[0,jj] = sig
                sig_arr[1,jj] = np.sqrt(res.hess_inv[2][2]/(err_2*8*sig*self.n_t_pts))
                omega_arr[0,jj] = omega
                omega_arr[1,jj] = np.sqrt(res.hess_inv[3][3]/(err_2*self.n_t_pts))
                phase_arr[0,jj] = cep/np.pi
                phase_arr[1,jj] = np.sqrt(res.hess_inv[4][4]/(err_2*np.pi*self.n_t_pts))
            else:
                print("bad fit! clust={}, j={}".format(clust,j))

        #we want to use spatial symmetry to infer the points which have not been collected
        i_zer = len(good_zs)
        for ii, zz in enumerate(good_zs):
            if zz >= 0:
                i_zer = ii
                break
        #keep track of the last point that satisfies z<=0
        i_cent = i_zer-1
        new_size = 2*i_zer
        #if zero is included then this is a special case where there is an odd number of points
        if i_zer < len(good_zs) and good_zs[i_zer] == 0:
            i_cent = i_zer
            new_size += 1
        print(phase_arr[0])
        print(phase_arr[1])
        print("new_size={}, good_z.len={}, i_zer={}, i_cent={}".format(new_size, len(good_zs), i_zer,i_cent))
        if new_size > len(good_zs):
            new_good_zs = np.zeros(new_size)
            amp_arr = np.resize(amp_arr, (2, new_size))
            t_0_arr = np.resize(t_0_arr, (2, new_size))
            sig_arr = np.resize(sig_arr, (2, new_size))
            omega_arr = np.resize(omega_arr, (2, new_size))
            phase_arr = np.resize(phase_arr, (3, new_size))
            #average data points to the right of the center
            for ii in range(len(good_zs)-i_zer):
                amp_arr[0,i_cent-ii] = (amp_arr[0,i_cent-ii] + amp_arr[0,i_zer+ii])/2
                amp_arr[1,i_cent-ii] = np.sqrt(amp_arr[1,i_cent-ii]**2 + amp_arr[1,i_zer+ii]**2)
                t_0_arr[0,i_cent-ii] = (t_0_arr[0,i_cent-ii] + t_0_arr[0,i_zer+ii])/2
                t_0_arr[1,i_cent-ii] = np.sqrt(t_0_arr[1,i_cent-ii]**2 + t_0_arr[1,i_zer+ii]**2)
                sig_arr[0,i_cent-ii] = (sig_arr[0,i_cent-ii] + sig_arr[0,i_zer+ii])/2
                sig_arr[1,i_cent-ii] = np.sqrt(sig_arr[1,i_cent-ii]**2 + sig_arr[1,i_zer+ii]**2)
                omega_arr[0,i_cent-ii] = (omega_arr[0,i_cent-ii] + omega_arr[0,i_zer+ii])/2
                omega_arr[1,i_cent-ii] = np.sqrt(omega_arr[1,i_cent-ii]**2 + omega_arr[1,i_zer+ii]**2)
                phase_arr[0,i_zer-ii] = (phase_arr[0,i_cent-ii] + phase_arr[0,i_zer+ii])/2
                phase_arr[1,i_zer-ii] = np.sqrt(phase_arr[1,i_cent-ii]**2 + phase_arr[1,i_zer+ii]**2)
            for ii in range(new_size-i_zer):
                new_good_zs[i_cent-ii] = good_zs[i_cent-ii]
                new_good_zs[i_zer+ii] = -good_zs[i_cent-ii]
                amp_arr[0,i_zer+ii] = amp_arr[0,i_cent-ii]
                amp_arr[1,i_zer+ii] = amp_arr[1,i_cent-ii]
                t_0_arr[0,i_zer+ii] = t_0_arr[0,i_cent-ii]
                t_0_arr[1,i_zer+ii] = t_0_arr[1,i_cent-ii]
                sig_arr[0,i_zer+ii] = sig_arr[0,i_cent-ii]
                sig_arr[1,i_zer+ii] = sig_arr[1,i_cent-ii]
                omega_arr[0,i_zer+ii] = omega_arr[0,i_cent-ii]
                omega_arr[1,i_zer+ii] = omega_arr[1,i_cent-ii]
                phase_arr[0,i_zer+ii] = phase_arr[0,i_cent-ii]
                phase_arr[1,i_zer+ii] = phase_arr[1,i_cent-ii]
        else:
            new_good_zs = np.array(good_zs[:new_size])
            amp_arr = np.resize(amp_arr, (2, new_size))
            t_0_arr = np.resize(t_0_arr, (2, new_size))
            sig_arr = np.resize(sig_arr, (2, new_size))
            omega_arr = np.resize(omega_arr, (2, new_size))
            phase_arr = np.resize(phase_arr, (3, new_size))
        t_dif = time.clock_gettime_ns(time.CLOCK_MONOTONIC) - t_start
        print("Completed optimizations in {:.5E} ns, average time per eval: {:.5E} ns".format(t_dif, t_dif/n_evals))
        return new_good_zs, amp_arr, sig_arr, omega_arr, phase_arr

    def get_junc_bounds(self):
        return self.geom.meep_len_to_um(self.geom.l_junc - self.geom.z_center), self.geom.meep_len_to_um(self.geom.r_junc - self.geom.z_center)

    def lookup_fits(self, clust_name, recompute=False):
        '''Load the fits from the pickle file located at fname or perform the fits if it doesn't exist. Return the results'''
        data_name = '{}/dat_{}'.format(args.prefix, clust_name)
        if recompute or not os.path.exists(data_name):
            #figure out data by phase fitting
            dat_xs, amp_arr, sig_arr, omega_arr, phase_arr = self.read_cluster(clust_name)
            with open(data_name, 'wb') as fh:
                pickle.dump([dat_xs, amp_arr, sig_arr, omega_arr, phase_arr], fh)
            return dat_xs, amp_arr, sig_arr, omega_arr, phase_arr  
        else:
            with open(data_name, 'rb') as fh:
                vals = pickle.load(fh)
            return vals[0], vals[1], vals[2], vals[3], vals[4]

    def get_field_amps(self, x_pts, z, diel_const, vac_wavelen=0.7, n_modes=HIGHEST_MODE):
        length = self.geom.meep_len_to_um(self.geom.r_junc - self.geom.l_junc)
        #these are some useful constants that we define. Note that lc is just the wavelength inside the dielectric, kc_sq is the square of 2pi/lc and l_rat_sq is (lc/2l)^2 where l is the length of the sample. We derived the expression beta^2 = kc^2 (1 - l_rat_sq*n^2), so neg_beta_sq is just -beta^2.
        #lc = vac_wavelen/np.sqrt(np.real(diel_const))
        lc = vac_wavelen/np.real(diel_const)
        kc_sq = (2*np.pi / lc)**2
        l_rat_sq = 0.25*(lc/length)**2
        fields = 1j*np.zeros(len(x_pts))
        print("length={}, L_c={}, kc^2={}, ratio^2={}".format(length, lc, kc_sq, l_rat_sq))
        phaseless_comps = self.freq_comps*np.exp(1j*self.get_phase())
        for k in range(n_modes):
            n = 2*k + 1
            neg_beta_sq = kc_sq*(l_rat_sq*(n**2) - 1)
            print("beta^2={}".format(-neg_beta_sq))
            #real beta implies propagation, no decay
            if z > 0 and neg_beta_sq > 0:
                fields = fields + ( np.sin(n*np.pi*(x_pts+length/2)/length)/n )*np.sum( phaseless_comps*np.exp(-z*np.sqrt(neg_beta_sq)) )*self.df
            else:
                fields += np.sin(n*np.pi*(x_pts+length/2)/length)*np.sum(phaseless_comps)*self.df/n
        #the factor of 0.5 comes because only half of the field is propagating towards the junction, the other half is going away
        return 4*0.5*self.get_amplitude()*np.real(fields)/np.pi

    def get_field_phases(self, x_pts, z, diel_const, inc_phase, n_modes=HIGHEST_MODE):
        length = self.geom.meep_len_to_um(self.geom.r_junc - self.geom.l_junc)
        alphas = (self.freq_range**2+0j)/0.08988
        fours = np.exp(inc_phase)*np.ones(len(x_pts))
        for k in range(n_modes):
            n = 2*k + 1
            betas = np.sqrt(alphas*diel_const - (n*np.pi/length)**2)
            fours = fours + self.df*np.sum(np.exp(-1j*betas*z))*np.sin(n*np.pi*(x_pts+length/2)/length)/n
        return np.array([phases.fix_angle(x) for x in np.arctan2(np.imag(fours), np.real(fours))/np.pi])

    def get_n_component(self, x_pts, amps, n):
        length = self.geom.meep_len_to_um(self.geom.r_junc - self.geom.l_junc)
        inner = np.sin(n*np.pi*(x_pts+length/2)/length)*amps
        #use trapezoid rule to integrate
        ret = 0
        for i in range(1, len(x_pts)):
            ret += 0.5*(x_pts[i]+x_pts[i-1])*(inner[i]+inner[i-1])
        return 2*ret/length #2/L normalization factor since \int_0^L sin^2(npix/L)dx = L/2

    def get_field_amps_jac(self, x_pts, z, diel_const, n_modes=3):
        length = self.geom.meep_len_to_um(self.geom.r_junc - self.geom.l_junc)
        fields = 1j*np.zeros(len(x_pts))
        alphas = (self.freq_range**2+0j)/0.08988
        for k in range(n_modes):
            n = 2*k + 1
            beta = np.sqrt(alphas*diel_const - (n*np.pi/length)**2)
            #real beta implies propagation, no decay
            if z > 0:
                fields = fields + ( np.sin(n*np.pi*(x_pts+length/2)/length)/n )* \
                np.sum( 1j*z*self.freq_comps*alphas*np.exp(1j*z*beta)/beta )*self.df
        return 2*np.real(fields)/np.pi

    #def fit_eps(self, x_pts, amps, z):
    def fit_eps(self, inds=[]):
        import scipy.optimize as opt
        #ff = lambda x: np.sum( (self.get_field_amps(x_pts, z, x[0]) - amps)**2 )
        #select which clusters to use based on caller argument
        if inds == -1 or len(inds) == 0:
            inds = range(1, len(pf.clust_names))
        clusts = [pf.clust_names[ind] for ind in inds]

        junc_bounds = self.get_junc_bounds()
        length = self.geom.meep_len_to_um(self.geom.r_junc - self.geom.l_junc)
        slice_zs = []
        slice_xs = []
        slice_amps = []
        avg_eps_est = 0
        for clust in clusts:
            dat_xs, amp_arr, _, _ = self.lookup_fits(clust)
            #the TE mode expansion is only valid inside the junction, so we need to truncate the arrays
            lowest_valid = -1
            highest_valid = -1
            slice_zs.append(self.get_clust_location(clust))
            '''high_off = (junc_bounds[1]-junc_bounds[0])/(2*HIGHEST_MODE-1)
            fit_bounds = (junc_bounds[0]+high_off, junc_bounds[1]-high_off)'''
            fit_bounds = junc_bounds
            for j, x in enumerate(dat_xs):
                if lowest_valid < 0 and x >= fit_bounds[0]:
                    lowest_valid = j
                if x <= fit_bounds[1]:
                    highest_valid = j
            slice_xs.append(dat_xs[lowest_valid:highest_valid+1])
            slice_amps.append(amp_arr[:, lowest_valid:highest_valid+1])
            comp_1 = self.get_n_component(slice_xs[-1], slice_amps[-1][0], 1)
            comp_3 = self.get_n_component(slice_xs[-1], slice_amps[-1][0], 3)/3 #divide by component along this basis vector
            comp_5 = self.get_n_component(slice_xs[-1], slice_amps[-1][0], 5)/5
            r31 = (comp_3/comp_1)**2
            r53 = (comp_5/comp_3)**2
            print("n=1: %f" % comp_1)
            print("n=3: %f" % comp_3)
            print("n=5: %f" % comp_5)
            #estimate epsilon based on ratios between different components under varying assumptions of propagation or decay
            est_1p3d = (np.pi*0.2998/(FREQ_0*length))**2*(9-r31)
            est_1d3d = (np.pi*0.2998/(FREQ_0*length))**2*(r31-9)/(r31-1)
            est_3p5d = (np.pi*0.2998/(FREQ_0*length))**2*(9*r53-25)
            est_3d5d = (np.pi*0.2998/(FREQ_0*length))**2*(9*r53-25)/(r53-1)
            print("\tcomp_3/comp_1 = {}\n\t--> (1p3d) estimated eps = {}\n\t--> (1d3d) estimated eps = {}".format(np.sqrt(r31), est_1p3d, est_1d3d))
            print("\tcomp_5/comp_3 = {}\n\t--> (1p3d) estimated eps = {}\n\t--> (1d3d) estimated eps = {}".format(np.sqrt(r53), est_3p5d, est_3d5d))
            #based on which component is larger, average our estimates for epsilon
            if r31 < 1:
                avg_eps_est += est_1p3d
            else:
                avg_eps_est += est_1d3d
            if r53 <= 1:
                avg_eps_est += est_3p5d
            else:
                avg_eps_est += est_3d5d        
        if avg_eps_est <= 0:
            avg_eps_est = EPS_0
        else:
            avg_eps_est /= 2*len(clusts)
        print("avg_eps = %f" % avg_eps_est)

        #square error which we seek to minimize
        def ff(eps):
            ret = 0
            for z, dat_x, dat_amp in zip(slice_zs, slice_xs, slice_amps):
                ret += np.sum( (self.get_field_amps(dat_x, z, eps) - dat_amp[0])**2 )
            return ret
        def jac_ff(eps):
            ret = 0
            for z, dat_x, dat_amp in zip(slice_zs, slice_xs, slice_amps):
                ret += np.sum( 2*(self.get_field_amps(dat_x, z, eps) - dat_amp[0])*self.get_field_amps_jac(dat_x, z, eps) )
            return ret
        print((ff(avg_eps_est+0.001)-ff(avg_eps_est-0.001))/0.002, jac_ff(avg_eps_est))
        print("sq_err(x0)=%f" % ff(avg_eps_est))
        return opt.minimize(ff, avg_eps_est, jac=jac_ff), avg_eps_est

    def find_avg_vals(self, clust):
        dat_xs, amp_arr, sig_arr, phase_arr = self.lookup_fits(clust)
        field_0 = np.sum(amp_arr[0]) / amp_arr.shape[1]
        sig_0 = np.sum(sig_arr[0]) / sig_arr.shape[1]
        phase_0 = np.sum(phase_arr[0]) / sig_arr.shape[1]
        return field_0, sig_0, phase_0
    
def get_axis(axs_list, ind):
    #matplotlib is annoying and the axes it gives have a different type depending on the column
    if N_COLS == 1:
        return axs_list[ind]
    else:
        return axs_amp[ind//N_COLS,ind%N_COLS]

def make_theory_plots(pf):
    #figure out the edges of the junction
    junc_bounds = pf.get_junc_bounds()
    l_gold = [pf.x_range[0], junc_bounds[0]]
    r_gold = [junc_bounds[1], pf.x_range[-1]]
    #initialize plots
    fig_amp, axs_amp = plt.subplots(pf.n_clusts//N_COLS, N_COLS)
    x_pts = np.linspace(junc_bounds[0], junc_bounds[1], num=100)
    for i, clust in enumerate(pf.clust_names):
        z = pf.get_clust_location(clust)
        amps = pf.get_field_amps(x_pts, z, 0.5, 0.1, EPS_0)
        #actually make the plot
        tmp_axs = get_axis(axs_amp, i)
        tmp_axs.plot(x_pts, amps)
        tmp_axs.set_xlim(pf.x_range)
        tmp_axs.fill_between(l_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs.fill_between(r_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        #tmp_axs.set_ylim(AMP_RANGE)
    fig_amp.savefig(args.prefix+"/amps_theory.pdf")

def make_fits(pf, axs_mapping=None):
    #use a default set of axes if the user didn't supply one or the supplied axes were invalid
    if axs_mapping is None:
        axs_mapping = range(len(pf.clust_names))
    n_mapped = max(axs_mapping)+1
    if n_mapped > len(pf.clust_names):
        axs_mapping = range(len(pf.clust_names))
        n_mapped = pf.n_clusts
    #figure out the edges of the junction
    junc_bounds = pf.get_junc_bounds()
    l_gold = [pf.x_range[0], junc_bounds[0]]
    r_gold = [junc_bounds[1], pf.x_range[-1]]
    #initialize plots
    fig_amp, axs_amp = plt.subplots(n_mapped//N_COLS, N_COLS)
    fig_phs, axs_phs = plt.subplots(n_mapped//N_COLS, N_COLS)
    fig_omg, axs_omg = plt.subplots(n_mapped//N_COLS, N_COLS)
    #set up axes first
    for i in range(n_mapped):
        tmp_axs = get_axis(axs_amp, i)
        tmp_axs.set_xlim(pf.x_range)
        tmp_axs.set_ylim(AMP_RANGE)
        tmp_axs.fill_between(l_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs.fill_between(r_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        if i < len(pf.clust_names) - 1:
            tmp_axs.get_xaxis().set_visible(False)
        tmp_axs = get_axis(axs_phs, i)
        tmp_axs.set_xlim(pf.x_range)
        tmp_axs.set_ylim(PHI_RANGE)
        tmp_axs.fill_between(l_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs.fill_between(r_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs = get_axis(axs_omg, i)
        tmp_axs.set_xlim(pf.x_range)
        tmp_axs.set_ylim(OMG_RANGE)
        tmp_axs.fill_between(l_gold, OMG_RANGE[0], OMG_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs.fill_between(r_gold, OMG_RANGE[0], OMG_RANGE[1], color='yellow', alpha=0.3)
    #make a plot of average fits
    fig_amp.suptitle(r"amplitude across a junction $E_1=1/\sqrt{2}$, $E_2=0$")
    fig_phs.suptitle(r"phases across a junction $E_1=1/\sqrt{2}$, $E_2=0$")
    fig_fits, axs_fits = plt.subplots(2)
    fit_xs = np.zeros((3, n_mapped))
    x_cnt = np.linspace(junc_bounds[0], junc_bounds[1], num=100)
    axs_fits[0].scatter(fit_xs[0], fit_xs[1])
    axs_fits[1].scatter(fit_xs[0], fit_xs[2])
    axs_fits[0].set_ylim(AMP_RANGE)
    axs_fits[1].set_ylim([6, 7])
    #now perform a plot using the average of fits
    avg_fit = [np.sum(fit_xs[1][1:4])/fit_xs.shape[1]]
    for i, clust in zip(axs_mapping, pf.clust_names):
        dat_xs,amp_arr,sig_arr,omega_arr,phs_arr = pf.lookup_fits(clust)
        clust_z = pf.get_clust_location(clust)
        amps_fit = np.abs(pf.get_field_amps(x_cnt, clust_z, 10.48, vac_wavelen=0.7))
        phss_fit = np.abs(pf.get_field_phases(x_cnt, clust_z, 10.48, -np.pi/2))
        #plot amplitudes
        tmp_axs = get_axis(axs_amp, i)
        tmp_axs.plot(x_cnt, amps_fit, color='gray', linestyle=':')
        tmp_axs.scatter(dat_xs, amp_arr[0], s=3)
        tmp_axs.annotate(r"$z={}\mu$m".format(round(clust_z, 2)), (0.01, 0.68), xycoords='axes fraction')
        #plot phases
        tmp_axs = get_axis(axs_phs, i)
        #tmp_axs.plot(x_cnt, phss_fit)
        phase_th = -pf.get_phase()[0]/np.pi
        tmp_axs.plot([x_cnt[0], x_cnt[-1]], [phase_th, phase_th], color='gray', linestyle=':')
        tmp_axs.scatter(dat_xs, phs_arr[0], s=3)
        tmp_axs = get_axis(axs_omg, i)
        tmp_axs.annotate(r"$z={}\mu$m".format(round(clust_z, 2)), (0.01, 0.68), xycoords='axes fraction')
        omega_th = 2*np.pi*.299792458/pf.get_wavelen()[0] #2*pi*c/lambda
        tmp_axs.plot([x_cnt[0], x_cnt[-1]], [omega_th, omega_th], color='gray', linestyle=':')
        tmp_axs.scatter(dat_xs, omega_arr[0], s=3)

    #save the figures
    fig_amp.savefig(args.prefix+"/amps_theory.pdf")
    fig_phs.savefig(args.prefix+"/phases_theory.pdf")
    fig_omg.savefig(args.prefix+"/omega_sim.pdf")
    fig_fits.savefig(args.prefix+"/fit_plt.pdf")

def make_plots(pf):
    #figure out the edges of the junction
    junc_bounds = pf.get_junc_bounds()
    l_gold = [pf.x_range[0], junc_bounds[0]]
    r_gold = [junc_bounds[1], pf.x_range[-1]]
    #initialize plots
    fig_amp, axs_amp = plt.subplots(pf.n_clusts//N_COLS, N_COLS)
    fig_time, axs_time = plt.subplots(pf.n_clusts//N_COLS, N_COLS)
    fig_cep, axs_cep = plt.subplots(pf.n_clusts//N_COLS, N_COLS)
    for i, clust in enumerate(pf.clust_names):
        good_zs, amp_arr, sig_arr, phase_arr = pf.read_cluster(clust)
        #make two axis demo plot
        if i == 2:
            d_fig, d_ax1 = plt.subplots(figsize=(7,3))
            d_ax2 = d_ax1.twinx()
            #label things
            d_fig.suptitle("Phase and amplitude across a junction")
            d_ax1.set_ylabel(r"$\phi/\pi$", color='red')
            d_ax1.set_xlabel(r"${}$ ($\mu$m)".format(slice_dir))
            d_ax2.set_ylabel(r"Amplitude", color='blue')
            #phase plot
            d_ax1.scatter(good_zs, phase_arr[0], color='red')
            d_ax1.plot([good_zs[0], good_zs[-1]], [0, 0], color='black', linestyle=':')
            d_ax1.plot([good_zs[0], good_zs[-1]], [0.5, 0.5], color='gray', linestyle=':')
            d_ax1.plot([good_zs[0], good_zs[-1]], [-0.5, -0.5], color='gray', linestyle=':')
            d_ax1.set_xlim(pf.x_range)
            d_ax1.fill_between(l_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
            d_ax1.fill_between(r_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
            d_ax1.set_ylim(PHI_RANGE)
            d_ax1.set_yticks([-1, -0.5, 0, 0.5, 1])
            #amplitude plot
            d_ax2.scatter(good_zs, amp_arr[0], color='blue')
            d_ax2.set_xlim(pf.x_range)
            d_ax2.set_ylim(AMP_RANGE)
            d_fig.savefig(args.prefix+"/demo.pdf")

        tmp_axs_amp = get_axis(axs_amp, i)
        tmp_axs_time = get_axis(axs_time, i)
        tmp_axs_cep = get_axis(axs_cep, i)
        #make amplitude plot
        tmp_axs_amp.errorbar(good_zs, amp_arr[0], yerr=amp_arr[1], fmt='.', linestyle='')
        tmp_axs_amp.set_xlim(pf.x_range)
        tmp_axs_amp.fill_between(l_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_amp.fill_between(r_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_amp.set_ylim(AMP_RANGE)
        #make spread plot
        tmp_axs_time.errorbar(good_zs, sig_arr[0], yerr=sig_arr[1], label=r"$\sigma$", fmt='.', linestyle='')
        tmp_axs_time.set_xlim(pf.x_range)
        tmp_axs_time.fill_between(l_gold, SIG_RANGE[0], SIG_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_time.fill_between(r_gold, SIG_RANGE[0], SIG_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_time.set_ylim(SIG_RANGE)
        #make phase plot
        tmp_axs_cep.errorbar(good_zs, phase_arr[0], yerr=phase_arr[1], fmt='.', linestyle='')
        tmp_axs_cep.plot([good_zs[0], good_zs[-1]], [0, 0], color='black', linestyle=':')
        tmp_axs_cep.plot([good_zs[0], good_zs[-1]], [0.5, 0.5], color='gray', linestyle=':')
        tmp_axs_cep.plot([good_zs[0], good_zs[-1]], [-0.5, -0.5], color='gray', linestyle=':')
        tmp_axs_cep.set_xlim(pf.x_range)
        tmp_axs_cep.fill_between(l_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_cep.fill_between(r_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_cep.set_ylim(PHI_RANGE)
        tmp_axs_cep.set_yticks([-1, 1])
        #perform regression if desired
        if args.regression:
            if len(good_zs) > 1:
                cep_line = linregress(good_zs, phase_arr[0])
            else:
                cep_line = linregress([0,1],[0,0])#give bad dummy values
            print("x={}: phi(z)={}z+{}, phi(0.1)={}".format(x, cep_line.slope, cep_line.intercept, cep_line.slope*0.1+cep_line.intercept))
            tmp_axs_cep.plot(good_zs, cep_line.slope*np.array(good_zs) + cep_line.intercept, color='gray')
            tmp_axs_cep.annotate(r"$\phi={0:.2g}z+{1:.2g}$, $r^2={2:.2g}$".format(cep_line.slope, cep_line.intercept, cep_line.rvalue*cep_line.rvalue), (0.01, 0.84), xycoords='axes fraction')
        #hide x axis labels on everything
        if i < pf.n_clusts-N_COLS:
            tmp_axs_amp.get_xaxis().set_visible(False)
            tmp_axs_time.get_xaxis().set_visible(False)
            tmp_axs_cep.get_xaxis().set_visible(False)
        if i % N_COLS > 0:
            tmp_axs_amp.get_yaxis().set_visible(False)
            tmp_axs_time.get_yaxis().set_visible(False)
            tmp_axs_cep.get_yaxis().set_visible(False)
    fig_amp.suptitle("Amplitude as a function of depth")
    fig_time.suptitle("Pulse width as a function of depth")
    fig_cep.suptitle("Phase as a function of depth")
    fig_amp.supylabel(r"Amplitude (incident units)")
    fig_time.supylabel(r"$\sigma$ (fs)")
    fig_cep.supylabel(r"$\frac{\phi}{\pi}$")
    fig_amp.supxlabel(r"${}$ ($\mu$m)".format(slice_dir))
    fig_time.supxlabel(r"${}$ ($\mu$m)".format(slice_dir))
    fig_cep.supxlabel(r"${}$ ($\mu$m)".format(slice_dir))
    #fig_amp.tight_layout(pad=0.5)
    #fig_time.tight_layout(pad=0.5)
    #fig_cep.tight_layout(pad=0.5)
    fig_cep.savefig(args.prefix+"/phases.pdf")
    fig_amp.savefig(args.prefix+"/amps.pdf")
    fig_time.savefig(args.prefix+"/sigs.pdf")

pf = phase_finder(args.fname, args.gap_width, args.gap_thick)
#make_plots(pf)
#make_theory_plots(pf)
make_fits(pf, axs_mapping=[0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7])
