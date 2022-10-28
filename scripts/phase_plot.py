import phases
import utils
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
import time
import pickle
import os.path

plt.rc('font', size=14)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

N_COLS = 1

AMP_RANGE = (0, 0.55)
SIG_RANGE = (0, 10.0)
PHI_RANGE = (-1.0, 1.0)
PAD_FACTOR = 1.05

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
    def __init__(self, fname, width, height):
        '''read metadata from the h5 file'''
        self.geom = utils.Geometry("params.conf", gap_width=width, gap_thick=height)
        #open the h5 file and identify all the clusters
        self.f = h5py.File(fname, "r")
        self.clust_names = list(self.f.keys())[:-2]
        self.n_clusts = len(self.clust_names)
        #create a numpy array for time points, this will be shared across all points
        keylist = list(self.f['cluster_0'].keys())
        n_t_pts = len(self.f['cluster_0'][keylist[1]]['time'])
        self.n_post_skip = n_t_pts // phases.SKIP
        #figure out the extent of the simulation points in the spanned direction
        z_min = min(self.geom.meep_len_to_um(self.geom.l_junc -self.geom.z_center), \
                self.geom.meep_len_to_um(self.f['cluster_0']['locations'][slice_dir][0]-self.geom.z_center))
        self.x_range = (PAD_FACTOR*z_min, -PAD_FACTOR*z_min)
        #read information to figure out time units and step sizes
        t_min = self.f['info']['time_bounds'][0]
        t_max = self.f['info']['time_bounds'][1]
        self.dt = (t_max-t_min)/n_t_pts
        self.t_pts = np.linspace(t_min, t_max, num=n_t_pts)

    def get_clust_location(self, clust):
        '''return the location of the cluster with the name clust along the propagation direction (z if slice_dir=x x if slice_dir=z)'''
        return self.geom.meep_len_to_um(self.f[clust]['locations'][slice_name][0] - self.geom.t_junc)

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
        phase_arr = np.zeros((2,n_z_pts))
        good_zs = []
        phase_cor = 0
        for j, pt in enumerate(points):
            v_pts = np.array(self.f[clust][pt]['time']['Re'])
            # (dt*|dE/dt|)^2 evaluated at t=t_max (one factor of two comes from the imaginary part, the other from the gaussian. We do something incredibly hacky to account for the spatial part. We assume that it has an error identical to the time error, and add the two in quadrature hence the other factor of two.
            #err_2 = 8*np.max(np.diff(v_pts)**2)
            err_2 = 0.02
            print("clust={}, j={}\n\tfield error^2={}".format(clust, j, err_2))
            res,res_env = phases.opt_pulse_env(self.t_pts, v_pts, a_sigmas_sq=err_2, keep_n=2.5, fig_name="{}/fit_figs/fit_{}_{}".format(args.prefix,clust,j))
            n_evals += 1
            print("\tsquare errors = {}, {}\n\tx={}\n\tdiag(H^-1)={}".format(res.fun, res_env.fun, res.x, np.diagonal(res.hess_inv)))
            #only include this point if the fit was good
            if res.fun < 500.0:
                good_zs.append(zs[j])
                jj = len(good_zs)-1
                amp, t_0, sig, omega, cep = phases.get_params(res)
                amp_arr[0,jj] = amp
                amp_arr[1,jj] = np.sqrt(res_env.hess_inv[0][0]/(err_2*self.n_post_skip))/2
                t_0_arr[0,jj] = t_0
                t_0_arr[1,jj] = np.sqrt(res.hess_inv[1][1]/(err_2*self.n_post_skip))
                sig_arr[0,jj] = sig
                sig_arr[1,jj] = np.sqrt(res.hess_inv[2][2]/(err_2*8*sig*self.n_post_skip))
                phase_arr[0,jj] = cep/np.pi
                phase_arr[1,jj] = np.sqrt(res.hess_inv[4][4]/(err_2*np.pi*self.n_post_skip))
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
        if new_size > len(good_zs):
            new_good_zs = np.zeros(new_size)
            amp_arr = np.resize(amp_arr, (2, new_size))
            t_0_arr = np.resize(t_0_arr, (2, new_size))
            sig_arr = np.resize(sig_arr, (2, new_size))
            phase_arr = np.resize(phase_arr, (3, new_size))
            #average data points to the right of the center
            for ii in range(len(good_zs)-i_zer):
                amp_arr[0,i_cent-ii] = (amp_arr[0,i_cent-ii] + amp_arr[0,i_zer+ii])/2
                amp_arr[1,i_cent-ii] = np.sqrt(amp_arr[1,i_cent-ii]**2 + amp_arr[1,i_zer+ii]**2)
                t_0_arr[0,i_cent-ii] = (t_0_arr[0,i_cent-ii] + t_0_arr[0,i_zer+ii])/2
                t_0_arr[1,i_cent-ii] = np.sqrt(t_0_arr[1,i_cent-ii]**2 + t_0_arr[1,i_zer+ii]**2)
                sig_arr[0,i_cent-ii] = (sig_arr[0,i_cent-ii] + sig_arr[0,i_zer+ii])/2
                sig_arr[1,i_cent-ii] = np.sqrt(sig_arr[1,i_cent-ii]**2 + sig_arr[1,i_zer+ii]**2)
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
                phase_arr[0,i_zer+ii] = phase_arr[0,i_cent-ii]
                phase_arr[1,i_zer+ii] = phase_arr[1,i_cent-ii]
        else:
            new_good_zs = np.array(good_zs[:new_size])
            amp_arr = np.resize(amp_arr, (2, new_size))
            t_0_arr = np.resize(t_0_arr, (2, new_size))
            sig_arr = np.resize(sig_arr, (2, new_size))
            phase_arr = np.resize(phase_arr, (3, new_size))
        t_dif = time.clock_gettime_ns(time.CLOCK_MONOTONIC) - t_start
        print("Completed optimizations in {:.5E} ns, average time per eval: {:.5E} ns".format(t_dif, t_dif/n_evals))
        return new_good_zs, amp_arr, sig_arr, phase_arr

    def get_junc_bounds(self):
        return self.geom.meep_len_to_um(self.geom.l_junc - self.geom.z_center), self.geom.meep_len_to_um(self.geom.r_junc - self.geom.z_center)

    def lookup_fits(self, clust_name):
        '''Load the fits from the pickle file located at fname or perform the fits if it doesn't exist. Return the results'''
        data_name = '{}/dat_{}'.format(args.prefix, clust_name)
        if os.path.exists(data_name):
            with open(data_name, 'rb') as fh:
                vals = pickle.load(fh)
            return vals[0], vals[1], vals[2], vals[3]
        else:
            #figure out data by phase fitting
            dat_xs, amp_arr, sig_arr, phase_arr = pf.read_cluster(clust_name)
            with open(data_name, 'wb') as fh:
                pickle.dump([dat_xs, amp_arr, sig_arr, phase_arr], fh)
            return dat_xs, amp_arr, sig_arr, phase_arr

    def get_field_amps(self, x_pts, z, field_0, omega, diel_const, n_modes=3):
        #omega = 2*np.pi*.299792458 / vac_wavelen
        #kc = 2*np.pi / (vac_wavelen*np.sqrt(diel_const))
        #L_c = vac_wavelen*np.sqrt(diel_const)/2
        kc = omega*np.sqrt(diel_const)/0.2998
        length = self.geom.meep_len_to_um(self.geom.r_junc - self.geom.l_junc)
        #print("L_c={}, L={}".format(L_c, length))
        amps = np.zeros(len(x_pts))
        for k in range(n_modes):
            n = 2*k + 1
            #beta_sq = ((np.pi/L_c)**2)*(1 - (n*L_c/length)**2)
            beta_sq = kc**2 - (n*np.pi/length)**2
            #real beta implies propagation, no decay
            if beta_sq > 0 or z < 0:
                #print("n={}, z.beta_sq={}, kc={}".format(n, z*np.sqrt(beta_sq), kc))
                amps += np.sin(n*np.pi*(x_pts+length/2)/length)/n
            else:
                #print("n={}, z.beta_sq={}j, kc={}".format(n, z*np.sqrt(-beta_sq), kc))
                amps += np.sin(n*np.pi*(x_pts+length/2)/length)*np.exp(-z*np.sqrt(-beta_sq))/n
        return 4*field_0*np.real(amps)/np.pi

    def fit_eps(self, x_pts, amps, z, field_0):
        import scipy.optimize as opt
        ff = lambda x: np.sum( (self.get_field_amps(x_pts, z, x[0], x[1], 8.2) - amps)**2 )
        x0 = np.array([field_0, 6.56])
        print("sq_err(x0)=%f" % ff(x0))
        return opt.minimize(ff, x0)

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
        amps = pf.get_field_amps(x_pts, z, 0.5, 0.1, 3.5)
        #actually make the plot
        tmp_axs = get_axis(axs_amp, i)
        tmp_axs.plot(x_pts, amps)
        tmp_axs.set_xlim(pf.x_range)
        tmp_axs.fill_between(l_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs.fill_between(r_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        #tmp_axs.set_ylim(AMP_RANGE)
    fig_amp.savefig(args.prefix+"/amps_theory.pdf")

def make_fits(pf):
    #figure out the edges of the junction
    junc_bounds = pf.get_junc_bounds()
    l_gold = [pf.x_range[0], junc_bounds[0]]
    r_gold = [junc_bounds[1], pf.x_range[-1]]
    #initialize plots
    fig_amp, axs_amp = plt.subplots(pf.n_clusts//N_COLS, N_COLS)
    fig_fits, axs_fits = plt.subplots(2)
    fit_xs = np.zeros((3, pf.n_clusts))
    x_cnt = np.linspace(junc_bounds[0], junc_bounds[1], num=100)
    avg_field_0 = 0.0
    for i, clust in enumerate(pf.clust_names):
        dat_xs, amp_arr, sig_arr, phase_arr = pf.lookup_fits(clust)
        #the TE mode expansion is only valid inside the junction, so we need to truncate the arrays
        lowest_valid = -1
        highest_valid = -1
        for j, x in enumerate(dat_xs):
            if lowest_valid < 0 and x >= junc_bounds[0]:
                lowest_valid = j
            if x <= junc_bounds[1]:
                highest_valid = j
        dat_xs = dat_xs[lowest_valid:highest_valid+1]
        amp_arr = amp_arr[:, lowest_valid:highest_valid+1]
        #average the field
        if i == 0:
            avg_field_0 = np.sum(amp_arr[0]) / len(amp_arr[0])
        print(avg_field_0)
        z = pf.get_clust_location(clust)
        res = pf.fit_eps(dat_xs, amp_arr[0], z, avg_field_0)
        fit_xs[0, i] = z
        fit_xs[1, i] = res.x[0]
        fit_xs[2, i] = res.x[1]
        print(res)
        #amps_fit = pf.get_field_amps(x_cnt, z, 0.4, 6.56, 8.2)
        amps_fit_0 = pf.get_field_amps(x_cnt, z, avg_field_0, 6.56, 8.2)
        amps_fit = pf.get_field_amps(x_cnt, z, res.x[0], res.x[1], 8.2)
        #actually make the plot
        tmp_axs = get_axis(axs_amp, i)
        tmp_axs.plot(x_cnt, amps_fit_0, color='red')
        tmp_axs.plot(x_cnt, amps_fit)
        tmp_axs.scatter(dat_xs, amp_arr[0])
        tmp_axs.set_xlim(pf.x_range)
        tmp_axs.fill_between(l_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs.fill_between(r_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        #tmp_axs.set_ylim(AMP_RANGE)
    axs_fits[0].scatter(fit_xs[0], fit_xs[1])
    axs_fits[1].scatter(fit_xs[0], fit_xs[2])
    axs_fits[0].set_ylim(AMP_RANGE)
    axs_fits[1].set_ylim([6, 7])
    #now perform a plot using the average of fits
    avg_fit = [np.sum(fit_xs[1])/fit_xs.shape[1], np.sum(fit_xs[2])/fit_xs.shape[1]]
    print("average fit: {}".format(avg_fit))
    for i, clust in enumerate(pf.clust_names):
        amps_fit = pf.get_field_amps(x_cnt, z, avg_fit[0], avg_fit[1], 8.2)
        tmp_axs = get_axis(axs_amp, i)
        tmp_axs.plot(x_cnt, amps_fit, color='green')
    #save the figures
    fig_amp.savefig(args.prefix+"/amps_theory.pdf")
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
        good_zs, amp_arr, sig_arr, phase_arr = pf.lookup_fits(clust)
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
make_plots(pf)
#make_theory_plots(pf)
make_fits(pf)
