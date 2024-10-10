import phases_new as phases
from functools import partial
import argparse
import time
import pickle
import os.path
import numpy as np
import scipy.fft as fft
import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
#matplotlib.use('qtagg')
from mpl_toolkits.axes_grid1 import make_axes_locatable

np.random.seed(3141592)

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

AX_PAD = 0.05
C_RELW = 0.1

matplotlib.rc('figure', figsize=(20, 10))
plt.rc('font', size=18)
#plt.rc('xtick', labelsize=18)
#plt.rc('ytick', labelsize=18)

N_COLS = 1

AMP_RANGE = (0, 1.5)
OMG_RANGE = (0, 4.0)
SIG_RANGE = (0, 10.0)
PHI_RANGE = (-np.pi, np.pi)
N_FREQ_COMPS=100
EPS_0 = 200.6
#EPS_0 = 2.5

#The highest waveguide mode (n) that is considered in fits sin((2n-1)pi/L)
HIGHEST_MODE=3
PAD_FACTOR = 1.05

#parse arguments supplied via command line
def str_list(arg):
    return arg.split(',')
def flt_list(arg):
    return [float(s) for s in arg.split(',')]
def int_list(arg):
    return [int(s) for s in arg.split(',')]
parser = argparse.ArgumentParser(description='Fit time a_pts data to a gaussian pulse envelope to perform CEP estimation.')
parser.add_argument('--fname', type=str_list, help='h5 file to read', default=['field_samples_test.h5'])
parser.add_argument('--n-groups', type=int, help='for bowties there might be multiple groups of clusters which should be placed on the same axes', default=-1)
parser.add_argument('--recompute', action='store_true', help='If set, then phase estimation is recomputed. Otherwise, information is read from files saved in <prefix>', default=False)
parser.add_argument('--pad', action='store_true', help='If set, then plot the color bar for images', default=False)
parser.add_argument('--use-prior', action='store_true', help='If set, then plot then use priors for f0', default=False)
parser.add_argument('--save-fit-figs', action='store_true', help='If set, then intermediate plots of fitness are saved to <prefix>/fit_figs where <prefix> is specified by the --prefix flag.', default=False)
parser.add_argument('--herm-n', type=int, help='number of Hermite-Gaussian terms to include in fits', default=3)
parser.add_argument('--ang-n', type=int, help='half the degree of the polynomial used for fitting', default=3)
parser.add_argument('--gap-width', type=flt_list, help='junction width (in meep units)', default=[])
parser.add_argument('--gap-thick', type=flt_list, help='junction thickness (in meep units)', default=[])
parser.add_argument('--clust-range', type=int_list, help='If set, then plot the color bar for images', default=[])
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
parser.add_argument('--plot-pre-opt', action='store_true', help='If set include plots of the fitted pulse before full optimization. This is ignored if --point is not specified', default=False)
parser.add_argument('--point', type=str, help='if specified (as a comma separated tuple where the first index is the cluster number and the second is the point index) make a plot of fit parameters for a single point', default='')
parser.add_argument('--imtype', type=str, help='filetype to store for plotting defaults to pdf', default="pdf")
parser.add_argument('--opt-method', type=str, help='optimization method', default="trust-exact")
args = parser.parse_args()

NX = phases.HERM_OFF + args.herm_n + args.ang_n

def get_flims(vf, xlim, ylim):
    if xlim is None:
        xlim = [0, psig.f0 + 4*psig.sigma]
    if ylim is None:
        ylim = [0, np.max(np.abs(vf))*1.5]
        ylim[0] = -ylim[1]
    return xlim, ylim

def plot_fdom_e(ts, vs, psig, ax, xlim=None, ylim=None, ylabels=True):
    freqs = fft.rfftfreq(len(vs), d=ts[1]-ts[0])
    vf = fft.rfft(vs)
    xlim, ylim = get_flims(vf, xlim, ylim)

    #axs[0].annotate('a)', xy=(1,8), xytext=(0.07, 0.80), xycoords='figure fraction')
    #axs[1].annotate('b)', xy=(1,8), xytext=(0.07, 0.37), xycoords='figure fraction')
    #ax.annotate("R^2 = {:.3E}".format(psig.residual), xy=(0,0), xytext=(0.63, 0.64), xycoords='figure fraction')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('f (1/fs)')
    #ax.tick_params(axis ='y', labelcolor = 'black')
    ax.axvline(freqs[psig.fit_bounds[0]], color='gray', linestyle=':')
    ax.axvline(freqs[psig.fit_bounds[1]], color='gray', linestyle=':')
    ax.axvline(psig.f0, color='gray')


    #get the envelope function and multiply it by the appropriate rotation. By definition, the envelope is centered at zero frequency, so we have to roll the envelope back to its original position
    fit_field = psig.get_fspace(freqs)
    ax.plot(freqs, np.abs(vf), color='black', label='simulated magnitude')
    ax.plot(freqs, np.abs(fit_field), color='black', linestyle='-.', label='fitted magnitude')
    ax.plot(freqs, np.real(vf), color=PLT_COLORS[0], label='real')
    ax.plot(freqs, np.imag(vf), color=PLT_COLORS[1], label='imaginary')
    ax.fill_between(freqs, np.real(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[0], alpha=0.2)
    ax.fill_between(freqs, np.imag(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[1], alpha=0.2)
    ax.legend()
    #plot residuals
    residuals = np.abs(fit_field - vf)
    ax.scatter(freqs, residuals, color=PLT_COLORS[2])

def plot_fdom_vec(ts, vs, psig, ax, xlim=None, ylim=None, ylabels=True):
    #plot limits
    freqs = fft.rfftfreq(len(vs), d=ts[1]-ts[0])
    vf = fft.rfft(vs)
    xlim, ylim = get_flims(vf, xlim, ylim)
    ax.set_xlim([-xlim[1], xlim[1]])
    ax.set_ylim([8*ylim[0], 8*ylim[1]])
    ax.set_xlabel('f (1/fs)')
    #plot the vector potential envelope
    freq_center = freqs[len(freqs)//2]
    fit_field = psig.get_fenv(freqs-freq_center, field='A')
    ax.fill_between(freqs-freq_center, np.real(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[0], alpha=0.2)
    ax.fill_between(freqs-freq_center, np.imag(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[1], alpha=0.2)
    ax.plot(freqs-freq_center, np.abs(fit_field), color='black')
    ax.axvline(psig.f0-freq_center, color='gray')
    ax.axvline(-freq_center, color='gray')

def plot_raw_tdom(ts, vs, psig_opt, fname, psig_unopt=None):
    fig, axs = plt.subplots()
    #get the envelope and perform a fitting
    axs.plot(ts, vs, color='black', label='simulated data')
    freqs = fft.fftfreq(len(vs), d=ts[1]-ts[0])
    axs.plot(ts, np.real(psig_opt(ts)), color=PLT_COLORS[1], label="full optimization")
    if psig_unopt is not None:
        axs.plot(ts, np.real(psig_unopt(ts)), color=PLT_COLORS[0], label="initial guess")

    axs.axvline(psig_opt.t0, color='gray')
    #now use an envelope times a cosine
    tenv = np.abs(psig_opt.get_tenv(ts, field='E'))
    axs.fill_between(ts, -tenv, tenv, color=PLT_COLORS[1], alpha=0.2)
    axs.set_xlim([ts[0], ts[-1]])
    axs.set_ylim([-1, 1])
    axs.legend()
    fig.savefig(fname, bbox_inches='tight')

def log_prior(x, f0=0, scale=0, herm_n=1, ang_n=1):
    gd, hs = np.zeros(NX), np.zeros((NX,NX))
    gd[2] = scale*(f0 - x[2])
    hs[2,2] = -scale
    return -scale*(x[2] - f0)**2/2, gd, hs

class cluster_reader:
    '''
    Initialize a phase finder reading from the  specified by fname and a juction width and gap specified by width and height respectively.
    fname: h5 file to read time samples from
    width: the width of the junction
    height: the thickness of the junction
    '''
    def __init__(self, fname, herm_n, ang_n, use_prior=False, clust_range=[], prefix='.', recompute=False, save_fit_figs=False, pad=True):
        #otherwise, recalculate
        self.prefix = prefix
        self.herm_n = herm_n
        self.ang_n = ang_n
        #open the h5 file and identify all the clusters
        self.f = h5py.File(fname, "r")
        self.clust_names = []
        i = 0
        clust_len = 0
        #default to reading all points
        if len(clust_range) != 2:
            clust_range = [0, len(self.f.keys())]
        for key in self.f.keys():
            #make sure that the cluster has a valid name and it actually has points
            if 'cluster' in key and len(self.f[key]) > 1 and len(self.f[key]['locations']) > 0:
                if i == clust_range[1]:
                    break
                if i >= clust_range[0]:
                    self.clust_names.append(key)
                    this_clust_len = len(self.f[key]['locations'])
                    if clust_len != 0 and this_clust_len != clust_len:
                        print("Invalid file! clusters are not of the same length")
                        exit()
                    clust_len = this_clust_len
                i += 1
        self.n_clusts = len(self.clust_names)
        #read information to figure out time units and step sizes
        t_min = self.f['info']['time_bounds'][0]
        t_max = self.f['info']['time_bounds'][1]
        self.dt = self.f['info']['time_bounds'][2]
        self.n_t_pts = int(self.f['info']['n_time_points'][0])
        if pad:
            next_b2 = int(np.ceil(np.log(self.n_t_pts)/np.log(2)))
            self._pad_n = (1 << next_b2) - self.n_t_pts
            t_max += self.dt*self._pad_n
            self.n_t_pts += self._pad_n
        else:
            self._pad_n = 0
        self.t_pts = np.linspace(t_min, t_max, self.n_t_pts)

        #figure out the central frequency
        srcinfo = self.f['info']['sources']
        self.in_freq = .299792458 / srcinfo['wavelen'][0]
        self.freq_width = 0.5/np.pi/srcinfo['width'][0]
        if phases.verbose > 0:
            print("setting up paramater reader:")
            print("\ttime bounds = ({}, {}) in {} pts".format(t_min, t_max, self.n_t_pts))
            print("\tprior f0 = {}\\pm{}".format(self.in_freq, self.freq_width), flush=True)
        '''def log_prior(x):
            like = (x[2] - self.in_freq)**2/2/self.freq_width**2'''
        if use_prior:
            self.lp = partial(log_prior, f0=self.in_freq, scale=1/self.freq_width**2, herm_n=herm_n, ang_n=ang_n)
        else:
            self.lp = None

        #try loading precomputed data
        data_name = "{}/fit_data_{}".format(prefix, fname.split('/')[-1].split('.')[0])
        data_shape = (len(self.clust_names), clust_len, phases.HERM_OFF+herm_n+ang_n)
        if not recompute and os.path.exists(data_name):
            print("found cluster data ", data_name)
            with open(data_name, 'rb') as dat_f:
                pic = pickle.load(dat_f)
                if clust_range[0] >= pic[3][0] and clust_range[1] <= pic[3][1]:
                    ld, hd = clust_range[0] - pic[3][0], clust_range[1]
                    self._pts = pic[0][ld:hd,:,:]
                    self._raw_data = pic[1][ld:hd,:,:]
                    self.residuals = pic[2][ld:hd,:,:]
            return

        #otherwise, read all of the data
        self._pts = np.zeros( (len(self.clust_names), clust_len, 3) )
        self._raw_data = np.zeros( data_shape )
        self.residuals = np.zeros( (len(self.clust_names), clust_len, 2) )
        #find the start time for performance benchmarking
        if clust_len == 0:
            return
        t_start = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
        for i, clust in enumerate(self.clust_names):
            #set points
            points = list(self.f[clust].keys())[1:]
            self._pts[i,:,0] = np.array(self.f[clust]['locations']['x'])
            self._pts[i,:,1] = np.array(self.f[clust]['locations']['y'])
            self._pts[i,:,2] = np.array(self.f[clust]['locations']['z'])
            for j, point in enumerate(points):
                print("optimizing point", i, j, flush=True)
                v_pts = np.append(np.array(self.f[clust][point]['time']['Re']), np.zeros(self._pad_n))
                psig = phases.signal(self.t_pts, v_pts, herm_n=herm_n, ang_n=ang_n, log_prior=self.lp, method=args.opt_method)
                self._raw_data[i,j,:] = psig.x
                self.residuals[i,j,0] = psig.residual
                self.residuals[i,j,1] = psig.cost
                #save figures if specified
                if save_fit_figs:
                    fig, axs = plt.subplots()
                    plot_fdom_vec(self.t_pts, v_pts, psig, axs, xlim=[0,1])
                    fig.savefig("{}/fit_figs/ffit_{}_{}_vec.{}".format(self.prefix, clust, j, args.imtype))
                    fig, axs = plt.subplots()
                    plot_fdom_e(self.t_pts, v_pts, psig, axs, xlim=[0,1])
                    fig.savefig("{}/fit_figs/ffit_{}_{}.{}".format(self.prefix, clust, j, args.imtype))
                    plot_raw_tdom(self.t_pts, v_pts, psig, "{}/fit_figs/tfit_{}_{}_raw.{}".format(self.prefix, clust, j, args.imtype))
                    plt.close('all')
        t_dif = time.clock_gettime_ns(time.CLOCK_MONOTONIC) - t_start
        t_avg = t_dif/clust_len/len(self.clust_names)
        print("Completed optimizations in {:.5E} ns, average time per eval: {:.5E} ns".format(t_dif, t_avg))
        print("Average R^2: {:.5E}".format( np.mean(1-np.exp(self.residuals[:,:,0])) ))
        with open(data_name, 'wb') as dat_f:
            pickle.dump([self._pts, self._raw_data, self.residuals, clust_range], dat_f)

    def get_point_times(self, clust, ind):
        #fetch a list of points and their associated coordinates
        points = list(self.f[clust].keys())[1:]
        err_2 = 0.02
        v_pts = np.append(np.array(self.f[clust][points[ind]]['time']['Re']), np.zeros(self._pad_n))
        return v_pts, err_2

    def make_heatmap(self, ax, parameter, vlines=[], xlabels=[[],[]], ylabels=[[],[]], rng=[], cmap="magma", plot_cbar=False):
        imdat = None
        #if parameter == 'phase' or parameter == 't0' or parameter == 'amplitude':
        if parameter == 'amplitude':
            #these parameteters are special sinc they can't be directly read from the raw data but only accessed by a signal object
            imdat = np.zeros((self._raw_data.shape[0], self._raw_data.shape[1]))
            for i in range(self._raw_data.shape[0]):
                for j in range(self._raw_data.shape[1]):
                    sig = phases.signal(None, None, herm_n=self.herm_n, ang_n=self.ang_n, skip_opt=True, x0=self._raw_data[i,j,:])
                    if parameter == 'amplitude':
                        imdat[i,j] = np.max(np.abs(sig.get_tenv(self.t_pts)))
                    else:
                        t0, phi = sig.get_eff_t0_phi(self.t_pts)
                        if parameter == 't0':
                            imdat[i,j] = t0
                        else:
                            imdat[i,j] = phases.fix_angle(phi)
        elif parameter == 'R^2':
            imdat = 1-np.exp(self.residuals[:,:,0])
        elif parameter == 'residual':
            imdat = self.residuals[:,:,0]
        elif parameter == 'cost':
            imdat = self.residuals[:,:,1]
        elif parameter == 'phase':
            imdat = self._raw_data[:,:,0]
        elif parameter == 't0':
            imdat = self._raw_data[:,:,1]/2/np.pi
        elif parameter == 'f0':
            imdat = self._raw_data[:,:,2]
        elif parameter == 'sigma':
            imdat = 1/self._raw_data[:,:,2]
        elif parameter[:17] == 'Hermite amplitude':
            ind = int( parameter[17:] )
            if ind % 2 == 1:
                print("Error: only even Hermite Gaussian modes are allowed!")
                return
            ind = ind // 2
            imdat = self._raw_data[:,:,phases.HERM_OFF+ind]
        elif parameter[:15] == 'angle magnitude':
            ind = int( parameter[15:] )
            if ind % 2 == 0:
                print("Error: angle can only contain odd powers of dd!")
                return
            ind = ind // 2
            imdat = self._raw_data[:,:,phases.HERM_OFF+self.herm_n+ind]
        else:
            print("Error: unrecognized parameter!")
        imdat = np.append(imdat, np.flip(imdat, axis=1), axis=1)

        if len(xlabels[0]) == 0:
            ax.get_xaxis().set_visible(False)
        else:
            ax.set_xticks(xlabels[0])
            ax.set_xticklabels(xlabels[1])

        if len(ylabels[0]) == 0:
            ax.get_yaxis().set_visible(False)
        else:
            ax.set_yticks(ylabels[0])
            ax.set_yticklabels(ylabels[1])

        for xx in vlines:
            ax.axvline(x=xx, color='gray')
        #do plots
        if len(rng) == 2:
            im = ax.imshow(imdat, vmin=rng[0], vmax=rng[1], cmap=cmap)
        else:
            im = ax.imshow(imdat, cmap=cmap)

        if plot_cbar:
            pos = ax.get_position()
            cax = ax.get_figure().add_axes([pos.x1+pos. width*AX_PAD,pos.y0, pos.width*C_RELW, pos.height])
            ax.get_figure().colorbar(im, cax=cax)
            cax.set_ylabel(parameter) 
        return im

def fetch_ax(axs, i):
    if type(axs) == np.ndarray:
        return axs[i]
    return axs

if args.point == '':
    n_names = len(args.fname)

    phs_colors = [(0.18,0.22,0.07), (0.36,0.22,0.61), (1,1,1), (0.74,0.22,0.31), (0.18,0.22,0.07)]
    phase_cmap = LinearSegmentedColormap.from_list('phase_colors', phs_colors, N=100)

    am_fig, am_axs = plt.subplots(1, n_names, gridspec_kw={'wspace': AX_PAD})
    t0_fig, t0_axs = plt.subplots(1, n_names, gridspec_kw={'wspace': AX_PAD})
    ph_fig, ph_axs = plt.subplots(1, n_names, gridspec_kw={'wspace': AX_PAD})
    f0_fig, f0_axs = plt.subplots(1, n_names, gridspec_kw={'wspace': AX_PAD})
    rs_fig, rs_axs = plt.subplots(1, n_names, gridspec_kw={'wspace': AX_PAD})

    for i, name in enumerate(args.fname):
        print("loading ", name)
        cr = cluster_reader(name, args.herm_n, args.ang_n, use_prior=args.use_prior, clust_range=args.clust_range, prefix=args.prefix, recompute=args.recompute, save_fit_figs=args.save_fit_figs)

        if len(args.gap_width) == n_names:
            dx = cr._pts[0,1,0] - cr._pts[0,0,0]
            x_ext = int( (cr._pts[0,-1,0] - cr._pts[0,0,0])/dx )
            vlines = [x_ext - args.gap_width[i]/2/dx, x_ext + args.gap_width[i]/2/dx + 1]
        else:
            vlines = []
        f0_rng = [max(0,cr.in_freq-2*cr.freq_width),cr.in_freq+2*cr.freq_width]
        plot_cbar = (i == n_names-1)

        cr.make_heatmap(fetch_ax(am_axs, i), "amplitude", vlines=vlines, rng=[0,1], plot_cbar=plot_cbar)
        cr.make_heatmap(fetch_ax(t0_axs, i), "t0", vlines=vlines, rng=[0,cr.t_pts[-1]], plot_cbar=plot_cbar)
        cr.make_heatmap(fetch_ax(ph_axs, i), "phase", vlines=vlines, cmap=phase_cmap, rng=[-np.pi,np.pi], plot_cbar=plot_cbar)
        cr.make_heatmap(fetch_ax(f0_axs, i), "f0", vlines=vlines, rng=f0_rng, plot_cbar=plot_cbar)
        cr.make_heatmap(fetch_ax(rs_axs, i), "R^2", vlines=vlines, rng=[0,1], plot_cbar=plot_cbar)

    am_fig.savefig(os.path.join(args.prefix, "htmp_amp."+args.imtype), bbox_inches='tight')
    t0_fig.savefig(os.path.join(args.prefix, "htmp_t0."+args.imtype), bbox_inches='tight')
    ph_fig.savefig(os.path.join(args.prefix, "htmp_phase."+args.imtype), bbox_inches='tight')
    f0_fig.savefig(os.path.join(args.prefix, "htmp_f0."+args.imtype), bbox_inches='tight')
    rs_fig.savefig(os.path.join(args.prefix, "htmp_residuals."+args.imtype), bbox_inches='tight')

    plt.close(am_fig)
    plt.close(t0_fig)
    plt.close(ph_fig)
    plt.close(f0_fig)
    plt.close(rs_fig)
else:
    cr = cluster_reader(args.fname[0], args.herm_n, args.ang_n, use_prior=args.use_prior, prefix=args.prefix, recompute=args.recompute, save_fit_figs=args.save_fit_figs)
    #phases.signal._do_grad_tests(np.array([0.4]), 3, 2)
    point_arr = args.point.split(",")
    clust = "cluster_"+point_arr[0]
    j = int(point_arr[1])
    fig_name = "{}/{}_{}".format(args.prefix, clust, j)
    #read the data and set up the signal analyzer
    v_pts, _ = cr.get_point_times(clust, j)
    freqs = fft.rfftfreq(len(cr.t_pts), d=cr.t_pts[1]-cr.t_pts[0])
    psig_after  = phases.signal(cr.t_pts, v_pts, skip_opt=False, herm_n=args.herm_n, ang_n=args.ang_n, log_prior=cr.lp, method=args.opt_method)
    fit_series_1 = psig_after.get_fspace(freqs)

    print("Testing gradients:", end=" ")
    x0, _, _ = psig_after._guess_params_opt(freqs)
    t0 = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    success = phases.test_grads(lambda x: psig_after.residuals(x,freqs)[0], lambda x: psig_after.residuals(x,freqs)[1], x0)
    t1 = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    success *= phases.test_hess(lambda x: psig_after.residuals(x,freqs)[0], lambda x: psig_after.hess_res(x,freqs), x0)
    t2 = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    if success:
        print("\033[92mSuccess!\033[0m", end=" ")
    else:
        print("\033[31mFailure!\033[0m", end=" ")
    print("gradients took {} ms, Hessians took {} ms".format((t1-t0)*1e-6, (t2-t1)*1e-6))

    #set up the plot
    fig, axs = plt.subplots()
    axs.axvline(freqs[psig_after.fit_bounds[0]], color='gray', linestyle=':')
    axs.axvline(freqs[psig_after.fit_bounds[1]], color='gray', linestyle=':')
    axs.axvline(psig_after.f0, color=PLT_COLORS[0])
    axs.axvline(psig_after.f0, color=PLT_COLORS[1])
    axs.set_xlim(0, freqs[-1]/2)
    axs.set_ylim(-np.pi, np.pi)
    #plot each series and annotate
    axs.scatter(freqs, phases.fix_angle(psig_after.vfa), color='black', label="simulation")
    #get signals
    axs.plot(freqs, np.angle(fit_series_1), color=PLT_COLORS[1], label="full optimization")
    axs.annotate('$\\varphi = ${:.2f}, $t_0 = ${:.2f} fs'.format(psig_after.phi, psig_after.t0), xy=(0,10), xytext=(0.2, 0.80), xycoords='figure fraction')
    if args.plot_pre_opt:
        psig_before = phases.signal(cr.t_pts, v_pts, skip_opt=True, herm_n=args.herm_n, ang_n=args.ang_n, log_prior=cr.lp, run_grad_tests=True)
        fit_series_0 = psig_before.get_fspace(freqs)
        axs.plot(freqs, np.angle(fit_series_0), color=PLT_COLORS[0], label="initial guess")
        axs.legend()
    else:
        psig_before = None
    fig.savefig("{}_angle.{}".format(fig_name, args.imtype), bbox_inches='tight')
    plt.close(fig)

    #plot frequency space 
    if args.plot_pre_opt:
        fig,axs = plt.subplots(2)
        plot_fdom_vec(cr.t_pts, v_pts, psig_before, axs[0], xlim=[0,1])
        plot_fdom_vec(cr.t_pts, v_pts, psig_after, axs[1], xlim=[0,1], ylabels=False)
        fig.savefig("{}_fdom_vec.{}".format(fig_name, args.imtype), bbox_inches='tight')
        fig,axs = plt.subplots(2)
        plot_fdom_e(cr.t_pts, v_pts, psig_before, axs[0], xlim=[0,1])
        plot_fdom_e(cr.t_pts, v_pts, psig_after, axs[1], xlim=[0,1], ylabels=False)
        fig.savefig("{}_fdom.{}".format(fig_name, args.imtype), bbox_inches='tight')
    else:
        fig,axs = plt.subplots()
        plot_fdom_vec(cr.t_pts, v_pts, psig_after, axs, xlim=[0,1], ylabels=False)
        fig.savefig("{}_fdom_vec.{}".format(fig_name, args.imtype), bbox_inches='tight')
        fig,axs = plt.subplots()
        plot_fdom_e(cr.t_pts, v_pts, psig_after, axs, xlim=[0,1], ylabels=False)
        fig.savefig("{}_fdom.{}".format(fig_name, args.imtype), bbox_inches='tight')
    #plot time space 
    if args.plot_pre_opt:
        plot_raw_tdom(cr.t_pts, v_pts, psig_after, "{}_tdom.{}".format(fig_name, args.imtype), psig_unopt=psig_before)
    else:
        plot_raw_tdom(cr.t_pts, v_pts, psig_after, "{}_tdom.{}".format(fig_name, args.imtype))
