import phases_new as phases
import utils
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
matplotlib.use('qtagg')
from scipy.stats import linregress
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

matplotlib.rc('figure', figsize=(20, 10))
plt.rc('font', size=18)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

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
parser = argparse.ArgumentParser(description='Fit time a_pts data to a gaussian pulse envelope to perform CEP estimation.')
parser.add_argument('--fname', type=str, help='h5 file to read', default='field_samples.h5')
parser.add_argument('--n-groups', type=int, help='for bowties there might be multiple groups of clusters which should be placed on the same axes', default=-1)
parser.add_argument('--recompute', action='store_true', help='If set, then phase estimation is recomputed. Otherwise, information is read from files saved in <prefix>', default=False)
parser.add_argument('--plot-cbar', action='store_true', help='If set, then plot the color bar for images', default=False)
parser.add_argument('--plot-y-labels', action='store_true', help='If set, then plot the y axis labels', default=False)
parser.add_argument('--plot-x-labels', action='store_true', help='If set, then plot the y axis labels', default=False)
parser.add_argument('--plot-legend', action='store_true', help='If set, then plot the y axis labels', default=False)
parser.add_argument('--save-fit-figs', action='store_true', help='If set, then intermediate plots of fitness are saved to <prefix>/fit_figs where <prefix> is specified by the --prefix flag.', default=False)
parser.add_argument('--sym-mode', help='Type of symmeterization mode to use', default="sym-o")
parser.add_argument('--f0-mode', help='Type of central frequency estimate to use', default="avg")
parser.add_argument('--gap-width', type=float, help='junction width', default=0.1)
parser.add_argument('--gap-thick', type=float, help='junction thickness', default=0.2)
parser.add_argument('--diel-const', type=float, help='dielectric constant of material', default=3.5)
parser.add_argument('--lowpass', type=float, help='If optimizations fail, then the phase finder tries them again after applying a lowpass filter to the timeseries. This parameter specifies the strength of the lowpass filter. If this script is crashing, then try increasing this parameter.', default=1.0)
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
parser.add_argument('--point', type=str, help='if specified (as a comma separated tuple where the first index is the cluster number and the second is the point index) make a plot of fit parameters for a single point', default='')
parser.add_argument('--slice-dir', type=str, help='prefix to use when opening files', default='x')
args = parser.parse_args()

def plot_raw_fdom(ts, vs, psig, axs, xlim=None, ylim=None, ylabels=True):
    freqs = fft.rfftfreq(len(vs), d=ts[1]-ts[0])
    vf = fft.rfft(vs)
    if xlim is None:
        xlim = [0, psig.f0 + 4*psig.sigma]
    if ylim is None:
        ylim = [0, np.max(vf)*1.5]
        ylim[0] = -ylim[1]
    # setup labels and annotations 
    if ylabels:
        axs[0].set_ylabel('A(f)', color = 'black')
        axs[1].set_ylabel('E(f)', color = 'black')
    else:
        axs[0].get_yaxis().set_visible(False)
        axs[1].get_yaxis().set_visible(False)
    axs[0].set_xlim([-xlim[1], xlim[1]])
    axs[0].set_ylim([8*ylim[0], 8*ylim[1]])
    axs[1].set_xlim(xlim)
    axs[1].set_ylim(ylim)
    axs[1].set_xlabel('f (1/fs)')
    axs[1].tick_params(axis ='y', labelcolor = 'black')
    axs[1].axvline(freqs[psig.fit_bounds[0]], color='gray', linestyle=':')
    axs[1].axvline(freqs[psig.fit_bounds[1]], color='gray', linestyle=':')
    axs[1].axvline(psig.f0, color='gray')
    #plot the vector potential envelope
    freq_center = freqs[len(freqs)//2]
    fit_field = psig.get_fenv(freqs-freq_center, field='A')
    axs[0].fill_between(freqs-freq_center, np.real(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[0], alpha=0.2)
    axs[0].fill_between(freqs-freq_center, np.imag(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[1], alpha=0.2)
    axs[0].plot(freqs-freq_center, np.abs(fit_field), color='black')
    axs[0].axvline(psig.f0-freq_center, color='gray')
    axs[0].axvline(-freq_center, color='gray')
    #get the envelope function and multiply it by the appropriate rotation. By definition, the envelope is centered at zero frequency, so we have to roll the envelope back to its original position
    emags, eargs = psig.get_fspace(freqs)
    fit_field = emags*np.exp(1j*eargs)
    axs[1].plot(freqs, np.abs(vf), color='black')
    axs[1].plot(freqs, np.abs(fit_field), color='black', linestyle='-.')
    axs[1].plot(freqs, np.real(vf), color=PLT_COLORS[0], label='real')
    axs[1].plot(freqs, np.imag(vf), color=PLT_COLORS[1], label='imaginary')
    axs[1].fill_between(freqs, np.real(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[0], alpha=0.2)
    axs[1].fill_between(freqs, np.imag(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[1], alpha=0.2)
    #plot residuals
    r_ax = axs[1]
    '''r_ax = axs[1].twinx()
    r_ax.set_ylabel("residuals", color=PLT_COLORS[2])
    r_ax.set_xlim(xlim)
    r_ax.set_ylim([0, 2*ylim[1]])
    r_ax.tick_params(axis ='y', labelcolor = 'green')'''
    residuals = np.sqrt( emags**2 + np.abs(vf)**2 - 2*emags*np.abs(vf)*np.cos(eargs-np.angle(vf)) )
    r_ax.scatter(freqs, residuals, color=PLT_COLORS[2])

def plot_raw_tdom(ts, vs, psig_unopt, psig_opt, fname):
    fig, axs = plt.subplots()
    #get the envelope and perform a fitting
    axs.plot(ts, vs, color='black', label='simulated data')
    freqs = fft.fftfreq(len(vs), d=ts[1]-ts[0])
    if psig_unopt is not None:
        axs.plot(ts, np.real(psig_unopt(ts)), color=PLT_COLORS[0], label="initial guess")
    axs.plot(ts, np.real(psig_opt(ts)), color=PLT_COLORS[1], label="full optimization")
    axs.axvline(psig_opt.t0, color='gray')
    #now use an envelope times a cosine
    fit_pulse = np.real( psig_opt.get_tenv(ts, field='A') )
    tenv = np.abs(psig_opt.get_tenv(ts, field='E'))
    axs.fill_between(ts, -tenv, tenv, color=PLT_COLORS[1], alpha=0.2)
    #axs.fill_between(ts, np.zeros(ts.shape), fit_pulse, color=PLT_COLORS[1], alpha=0.2)
    fit_pulse *= np.sqrt(2)*np.sin(2*np.pi*psig_opt.f0*(ts - psig_opt.t0) + psig_opt.phi)
    axs.plot(ts[1:], np.diff(fit_pulse), color=PLT_COLORS[2])
    axs.set_xlim([ts[0], ts[-1]])
    axs.set_ylim([-1, 1])
    axs.legend()
    fig.savefig(fname, bbox_inches='tight')

class cluster_reader:
    '''
    Initialize a phase finder reading from the  specified by fname and a juction width and gap specified by width and height respectively.
    fname: h5 file to read time samples from
    width: the width of the junction
    height: the thickness of the junction
    '''
    def __init__(self, fname, clust_range=[], prefix='.', herm_n=3, ang_n=1, recompute=False, save_fit_figs=False):
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
        c_span = clust_range[1] - clust_range[0]
        for key in self.f.keys():
            #make sure that the cluster has a valid name and it actually has points
            if 'cluster' in key and len(self.f[key]) > 1 and len(self.f[key]['locations']) > 0:
                if i >= clust_range[0]:
                    #only include the specified range of clusters
                    if i >= c_span:
                        break
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
        self.n_t_pts = self.f['info']['n_time_points'][0]
        self.t_pts = np.linspace(t_min, t_max, self.n_t_pts)
        #figure out the central frequency
        self.in_freq = .299792458 / self.f['info']['sources']['wavelen'][0]
        self.freq_width = 0.5/np.pi/self.f['info']['sources']['width'][0]
        #try loading precomputed data
        data_name = "{}/dat_{}".format(prefix, fname.split('.')[0])
        data_shape = (len(self.clust_names), clust_len, phases.HERM_OFF+herm_n+ang_n)
        if not recompute and os.path.exists(data_name):
            with open(data_name, 'rb') as dat_f:
                pic = pickle.load(dat_f)
                self._pts = pic[0]
                self._raw_data = pic[1]
                self.residuals = pic[2]
                if self._raw_data.shape == data_shape:
                    return

        #otherwise, read all of the data
        self._pts = np.zeros( (len(self.clust_names), clust_len, 3) )
        self._raw_data = np.zeros( data_shape )
        self.residuals = np.zeros( (len(self.clust_names), clust_len) )
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
                print("optimizing point", i, j)
                v_pts = np.array(self.f[clust][point]['time']['Re'])
                psig = phases.signal(self.t_pts, v_pts, herm_n=herm_n, ang_n=ang_n)
                self._raw_data[i,j,:] = psig.x
                self.residuals[i,j] = psig.cost
                #save figures if specified
                if save_fit_figs:
                    fig, axs = plt.subplots(2)
                    plot_raw_fdom(self.t_pts, v_pts, psig, axs, xlim=[0,1])
                    fig.savefig("{}/fit_figs/ffit_{}_{}.svg".format(self.prefix, clust, j))
                    plot_raw_tdom(self.t_pts, v_pts, None, psig, "{}/fit_figs/tfit_{}_{}_raw.svg".format(self.prefix, clust, j))
        t_dif = time.clock_gettime_ns(time.CLOCK_MONOTONIC) - t_start
        t_avg = t_dif/clust_len/len(self.clust_names)
        print("Completed optimizations in {:.5E} ns, average time per eval: {:.5E} ns".format(t_dif, t_avg))
        with open(data_name, 'wb') as dat_f:
            pickle.dump([self._pts, self._raw_data, self.residuals], dat_f)

    def get_point_times(self, clust, ind):
        #fetch a list of points and their associated coordinates
        points = list(self.f[clust].keys())[1:]
        err_2 = 0.02
        return np.array(self.f[clust][points[ind]]['time']['Re']), err_2

    def make_heatmap(self, ax, parameter, vlines=[], xlabels=[[],[]], ylabels=[[],[]], rng=[], cmap="viridis", plot_cbar=True):
        imdat = None
        if parameter == 'phase' or parameter == 't0' or parameter == 'amplitude':
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
                            imdat[i,j] = phi
        elif parameter == 'residual':
            imdat = self.residuals
        elif parameter == 't0':
            sig = phases.signal(None, None, herm_n=self.herm_n, ang_n=self.ang_n, skip_opt=True, x0=self._raw_data[i,j,:])
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
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            cax.set_ylabel(parameter)
        return im
    
def plot_before_after_phase(ts, vs, fname):
    freqs = fft.rfftfreq(len(ts), d=ts[1]-ts[0])
    #get signals
    psig_unopt = phases.signal(ts, vs, skip_opt=True, ang_n=2)
    psig_opt = phases.signal(ts, vs, skip_opt=False, ang_n=2)
    dat_series = phases.fix_angle_seq(psig_unopt.vfa, scan_n=4)
    _, fit_series_0 = psig_unopt.get_fspace(freqs)
    _, fit_series_1 = psig_opt.get_fspace(freqs)
    #set up the plot
    fig, axs = plt.subplots()
    axs.axvline(freqs[psig_unopt.fit_bounds[0]], color='gray', linestyle=':')
    axs.axvline(freqs[psig_unopt.fit_bounds[1]], color='gray', linestyle=':')
    axs.axvline(psig_unopt.f0, color=PLT_COLORS[0])
    axs.axvline(psig_opt.f0, color=PLT_COLORS[1])
    scale = np.max(psig_unopt.vfa)-np.min(psig_unopt.vfa)
    axs.set_ylim(-scale, scale)
    #plot each series and annotate
    axs.plot(freqs, dat_series, color='black', label="simulation")
    axs.plot(freqs, fit_series_0, color=PLT_COLORS[0], label="initial guess")
    axs.plot(freqs, fit_series_1, color=PLT_COLORS[1], label="full optimization")
    axs.annotate('$\\varphi = ${:.2f}, $t_0 = ${:.2f} fs'.format(psig_opt.phi, psig_opt.t0), xy=(0,10), xytext=(0.2, 0.80), xycoords='figure fraction')
    axs.legend()
    fig.savefig(fname, bbox_inches='tight')
    return psig_unopt, psig_opt

if args.point == '':
    cr = cluster_reader(args.fname, prefix=args.prefix, recompute=args.recompute, save_fit_figs=args.save_fit_figs)
    dx = cr._pts[0,1,0] - cr._pts[0,0,0]
    x_ext = int( (cr._pts[0,-1,0] - cr._pts[0,0,0])/dx )
    vlines = [x_ext - args.gap_width/dx, x_ext + args.gap_width/dx]

    phs_colors = [(0.18,0.22,0.07), (0.36,0.22,0.61), (1,1,1), (0.74,0.22,0.31), (0.18,0.22,0.07)]
    cmap = LinearSegmentedColormap.from_list('phase_colors', phs_colors, N=100)

    fig, axs = plt.subplots(2)
    cr.make_heatmap(axs[0], "amplitude", vlines=vlines, cmap='magma')
    cr.make_heatmap(axs[1], "phase", vlines=vlines, cmap=cmap)
    #fig.tight_layout()
    fig.savefig(args.prefix+"/phase_amp_{}.svg".format(args.fname.split('.')[0]), bbox_inches='tight')
    plt.close(fig)
    fig, axs = plt.subplots()
    cr.make_heatmap(axs, "residual", vlines=vlines, cmap='magma', rng=[-5,0])
    fig.savefig(args.prefix+"/residual_{}.svg".format(args.fname.split('.')[0]), bbox_inches='tight')
    plt.close(fig)
else:
    #phases.signal._do_grad_tests(np.array([0.4]), 3, 2)
    point_arr = args.point.split(",")
    clust = "cluster_"+point_arr[0]
    j = int(point_arr[1])
    fig_name = "{}/fit_{}_{}".format(args.prefix, clust, j)
    cr = cluster_reader(args.fname, clust_range=[0,0], prefix=args.prefix, recompute=args.recompute, save_fit_figs=args.save_fit_figs)
    #read the data and set up the signal analyzer
    v_pts, _ = cr.get_point_times(clust, j)
    psig_before, psig_after = plot_before_after_phase(cr.t_pts, v_pts, "{}_param_est.svg".format(fig_name))
    #plot frequency space
    fig,axs = plt.subplots(2,2) 
    plot_raw_fdom(cr.t_pts, v_pts, psig_before, axs[:,0], xlim=[0,1]) 
    plot_raw_fdom(cr.t_pts, v_pts, psig_after, axs[:,1], xlim=[0,1], ylabels=False)
    fig.savefig("{}_raw_fdom.svg".format(fig_name), bbox_inches='tight')
    #plot time space 
    plot_raw_tdom(cr.t_pts, v_pts, psig_before, psig_after, "{}_raw_tdom.svg".format(fig_name))
