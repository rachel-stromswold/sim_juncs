import phases_new as phases
import utils
import argparse
import numpy as np
import scipy.fft as fft
import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
matplotlib.use('qtagg')
from scipy.stats import linregress
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
PHI_RANGE = (-1.0, 1.0)
SKW_RANGE = (-10.0, 50.0)
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

class waveguide_est:
    def __init__(self, width, height, pf):
        self.geom = utils.Geometry("params.conf", gap_width=width, gap_thick=height)
        #read source parameters from the phase finder
        self.vac_wavelen = pf.get_src_wavelen()[0]
        self.vac_pulse_width = pf.get_src_width()[0]
        self.vac_phase = pf.get_src_phase()[0]
        self.vac_amp = pf.get_src_amplitude()[0]
        #set object variables
        self.freq_0 = .299792458 / self.vac_wavelen
        self.df = 6/(self.vac_pulse_width*N_FREQ_COMPS)
        #create a linear range of frequencies
        self.freq_range = np.linspace(-self.freq_0-3/self.vac_pulse_width, -self.freq_0+3/self.vac_pulse_width, num=N_FREQ_COMPS)
        scale = self.vac_amp / np.sqrt(2*np.pi)
        self.freq_comps = scale*self.vac_pulse_width*np.exp(-1j*self.vac_phase-((self.freq_range+self.freq_0)*self.vac_pulse_width)**2/2)
        #get x ranges for plots
        z_min = min(self.geom.meep_len_to_um(self.geom.l_junc -self.geom.z_center), \
                self.geom.meep_len_to_um(pf.get_clust_span(0)[0]-self.geom.z_center))
        self.x_range = (PAD_FACTOR*z_min, -PAD_FACTOR*z_min)

    def reflect_pts(self, cr):
        '''mirror the cluster_result object cr and label the x axis appropriately
        '''
        return cr.mirror_pts(self.geom.z_center, self.geom.um_scale)

    def get_junc_bounds(self):
        return self.geom.meep_len_to_um(self.geom.l_junc - self.geom.z_center), self.geom.meep_len_to_um(self.geom.r_junc - self.geom.z_center)

    def setup_axes(self, ax, rng, label_x=True):
        '''Setup the pyplot axes ax
        ax: the axes to use
        rng: a tuple with at least two elements specifying upper and lower bounds for the axes'''
        #set axes ranges
        ax.set_xlim(self.x_range)
        ax.set_ylim(rng)
        #draw gold leads
        junc_bounds = self.get_junc_bounds()
        l_gold = [self.x_range[0], junc_bounds[0]]
        r_gold = [junc_bounds[1], self.x_range[-1]]
        ax.fill_between(l_gold, rng[0], rng[1], color='yellow', alpha=0.3)
        ax.fill_between(r_gold, rng[0], rng[1], color='yellow', alpha=0.3)
        ax.get_xaxis().set_visible(label_x)

    
def get_axis(axs_list, ind):
    #matplotlib is annoying and the axes it gives have a different type depending on the column
    if N_COLS == 1:
        return axs_list[ind]
    else:
        return axs_amp[ind//N_COLS,ind%N_COLS]

def find_grouping(pf, n_groups):
    '''returns: a tuple with a group length, number of groups and axs mapping giving the specified number of groups. This is found by dividing the clusters in pf into n_groups sets of subclusters. Each cluster is assigned an index in axs_mapping'''
    grp_len = pf.n_clusts
    if n_groups >= 1:
        grp_len = int(np.ceil(pf.n_clusts/n_groups))
    else:
        n_groups = 1
    axs_mapping = [i%grp_len for i in range(pf.n_clusts)]
    return grp_len, n_groups, axs_mapping

def make_heatmap(fg, ax, imdat, title, label, rng=None, cmap='viridis', vlines=[], xlabels=[[],[]], ylabels=[[],[]], w=-1):
    ax.set_title(title)
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
    if rng is None:
        im = ax.imshow(imdat, cmap=cmap)
    else:
        im = ax.imshow(imdat, vmin=rng[0], vmax=rng[1], cmap=cmap)
    if args.plot_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        cax.set_ylabel(label)
    fg.tight_layout()
    #adjust the size of the plot
    if w > 0:
        l = ax.figure.subplotpars.left
        r = ax.figure.subplotpars.right
        t = ax.figure.subplotpars.top
        b = ax.figure.subplotpars.bottom
        figw = float(w)/(r-l)
        figh = figw*(t-b)/(r-l)
        ax.figure.set_size_inches(figw, figh)

def make_heatmaps(pf, n_groups=-1):
    grp_len,n_groups,_ = find_grouping(pf, n_groups)
    cl_xs = []
    cl_amp = []
    cl_phs = []
    cl_ampr = []
    cl_skew = []
    cl_err = []
    wg = waveguide_est(args.gap_width, args.gap_thick, pf)
    for clust in pf.clust_names:
        cr = wg.reflect_pts(pf.lookup_fits(clust, recompute=args.recompute, save_fit_figs=args.save_fit_figs))
        #save x points, amplitudes and phases so that an average may be computed
        cl_xs.append(cr.xs)
        cl_amp.append(cr.get_amp())
        cl_phs.append(cr.get_phase())
        cl_ampr.append(cr.get_amp_ref())
        cl_skew.append(cr.get_skew())
        cl_err.append(cr.get_err_sq())
    cl_xs = np.array(cl_xs)
    cl_amp = np.array(cl_amp)
    cl_phs = np.array(cl_phs)
    cl_ampr = np.array(cl_ampr)
    cl_skew = np.array(cl_skew)
    cl_err = np.array(cl_err)

    #figure out the extent in the x and z directions
    x_min = cl_xs[0][0]
    x_max = cl_xs[0][-1]
    z_min = pf.get_clust_location(pf.clust_names[0])
    z_max = pf.get_clust_location(pf.clust_names[-1])
    dx = cl_xs[0][1] - cl_xs[0][0]
    xl = np.floor( (x_max - args.gap_width - x_min)/(2*dx) )+0.5
    xm = np.floor( (x_max - x_min)/(2*dx) )+0.5
    xr = np.floor( (x_max + args.gap_width - x_min)/(2*dx) )+0.5
    zt = 0
    zb = (pf.get_clust_location(pf.clust_names[-1]) - pf.get_clust_location(pf.clust_names[0]))

    phs_colors = [(0.18,0.22,0.07), (0.36,0.22,0.61), (1,1,1), (0.74,0.22,0.31), (0.18,0.22,0.07)]
    cmap = LinearSegmentedColormap.from_list('phase_colors', phs_colors, N=100)

    #set plot parameters which are reused
    title = "L={} μm λ={} μm".format(args.gap_thick, wg.vac_wavelen)
    n_cbar = 1
    wrs = [1]
    vlines = [xl, xr]
    xlabels = [[xl, xr], ["{}  ".format(-int(args.gap_width*500)), "  {}".format(int(args.gap_width*500))]]
    ylabels = [[], []]
    if args.plot_cbar:
        n_cbar = 2
        wrs = [32,1]
    if args.plot_y_labels:
        ylabels = [[1, len(cl_xs)-1], ["0", "{}".format(int(args.gap_thick*500))]]

    #save heatmaps of amplitudes and phases
    for i in range(n_groups):
        heat_fig, heat_ax = plt.subplots(2, gridspec_kw={'hspace':0.2})
        make_heatmap(heat_fig, heat_ax[0], 2*cl_amp[i*grp_len:(i+1)*grp_len], "", "amp.", rng=AMP_RANGE, cmap='magma', vlines=vlines, ylabels=ylabels, w=2.2)
        make_heatmap(heat_fig, heat_ax[1], cl_phs[i*grp_len:(i+1)*grp_len], "", r"φ/2π", rng=PHI_RANGE, cmap=cmap, vlines=vlines, xlabels=xlabels, ylabels=ylabels, w=2.2)
        heat_fig.savefig(args.prefix+"/heatmap_grp{}.svg".format(i), bbox_inches='tight')
        plt.close(heat_fig)
        #plot skews
        nfig = plt.figure()
        nax = plt.gca()
        make_heatmap(nfig, nax, cl_skew[i*grp_len:(i+1)*grp_len], "", "skewness", cmap='magma', rng=SKW_RANGE, vlines=vlines, xlabels=xlabels, ylabels=ylabels)
        nfig.savefig(args.prefix+"/heatmap_skew_grp{}.svg".format(i), bbox_inches='tight')
        plt.close(nfig)
        #plot errors
        heat_fig = plt.figure()
        heat_ax = plt.gca()
        make_heatmap(heat_fig, heat_ax, cl_err[i*grp_len:(i+1)*grp_len], "", "fit error (arb. units)", vlines=vlines, xlabels=xlabels, ylabels=ylabels)
        heat_fig.savefig(args.prefix+"/heatmap_err_grp{}.svg".format(i))
        plt.close(heat_fig)
        #plot reflected phases
        heat_fig, heat_ax = plt.subplots(2)
        make_heatmap(heat_fig, heat_ax[0], cl_phs[i*grp_len:(i+1)*grp_len], "Phases", r"φ/2π$", rng=PHI_RANGE, cmap='twilight_shifted')
        make_heatmap(heat_fig, heat_ax[1], cl_ampr[i*grp_len:(i+1)*grp_len], "Reflected amplitudes", r"φ/2π$", rng=AMP_RANGE, cmap='magma')
        heat_fig.savefig(args.prefix+"/heatmap_phs_grp{}.svg".format(i))
        plt.close(heat_fig)

    #save a figure of the average phase
    n_x_pts = len(cl_xs[0])
    avg_phs = np.zeros((n_groups, n_x_pts))
    tot_amp = np.zeros((n_groups, n_x_pts))
    for i in range(n_groups):
        for j in range(grp_len):
            k = i*grp_len + j
            if len(cl_xs[k]) == n_x_pts:
                tot_amp[i] = tot_amp[i] + cl_amp[k]
                avg_phs[i] = avg_phs[i] + cl_amp[k]*cl_phs[k]
    avg_phs = avg_phs / tot_amp
    tot_amp = tot_amp / grp_len
    avg_fig, avg_ax = plt.subplots(2)
    #setup the axes with the gold leads and a dashed line with incident phase
    wg.setup_axes(avg_ax[0], PHI_RANGE)
    wg.setup_axes(avg_ax[1], AMP_RANGE)
    phs_th = -pf.get_src_phase()[0]/np.pi
    avg_ax[0].plot([cl_xs[0][0], cl_xs[0][-1]], [phs_th, phs_th], color='gray', linestyle=':')
    avg_ax[0].set_ylabel(r"<φ>/π")
    avg_ax[1].set_ylabel(r"<E_0> (arb. units)")
    avg_ax[1].set_xlabel(r"x μm")
    #plot averages
    for i in range(n_groups):
        avg_ax[0].scatter(cl_xs[0], avg_phs[i])
        avg_ax[1].scatter(cl_xs[0], tot_amp[i])
    avg_fig.savefig(args.prefix+"/avgs.pdf")

def config_axs(axs, xlab, ylab, leg_loc="upper right", leg_ncol=-1, leg_alpha=0.8, w=-1, h=-1):
    if args.plot_x_labels:
        axs.set_xlabel(xlab)
    else:
        axs.get_xaxis().set_visible(False)
    if args.plot_y_labels:
        axs.set_ylabel(ylab)
    else:
        axs.get_yaxis().set_visible(False)
    if args.plot_legend:
        if leg_ncol > 0:
            axs.legend(loc=leg_loc, ncol=leg_ncol, framealpha=leg_alpha)
        else:
            axs.legend(loc=leg_loc)
    #adjust the size of the plot
    if w > 0 and h > 0:
        l = axs.figure.subplotpars.left
        r = axs.figure.subplotpars.right
        t = axs.figure.subplotpars.top
        b = axs.figure.subplotpars.bottom
        figw = float(w)/(r-l)
        figh = float(h)/(t-b)
        axs.figure.set_size_inches(figw, figh)

def plot_raw_fdom(psig, axs, xlim=None, ylim=None, ylabels=True):
    freqs = fft.rfftfreq(len(psig.v_pts), d=psig.dt)
    vf = fft.rfft(psig.v_pts)
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
    #fit_field = psig.get_fenv(freqs-psig.f0)*np.exp( 1j*(psig.phi - 2*np.pi*psig.t0*freqs) )
    emags, eargs = psig.get_fspace(freqs)
    fit_field = emags*np.exp(1j*eargs)
    axs[1].plot(freqs, np.abs(vf), color='black')
    axs[1].plot(freqs, np.abs(fit_field), color='black', linestyle='-.')
    axs[1].plot(freqs, np.real(vf), color=PLT_COLORS[0], label='real')
    axs[1].plot(freqs, np.imag(vf), color=PLT_COLORS[1], label='imaginary')
    axs[1].fill_between(freqs, np.real(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[0], alpha=0.2)
    axs[1].fill_between(freqs, np.imag(fit_field), np.zeros(freqs.shape), color=PLT_COLORS[1], alpha=0.2)
    axs[1].scatter(freqs, emags**2 + np.abs(vf)**2 - 2*emags*np.abs(vf)*np.cos(eargs-np.angle(vf)))

def plot_raw_tdom(ts, psig_unopt, psig_opt, fname):
    fig, axs = plt.subplots()
    #get the envelope and perform a fitting
    axs.plot(ts, psig_unopt.v_pts, color='black', label='simulated data')
    freqs = fft.fftfreq(len(psig_unopt.v_pts), d=ts[1]-ts[0])
    axs.plot(ts, np.real(psig_unopt(ts)), color=PLT_COLORS[0], label="initial guess")
    axs.plot(ts, np.real(psig_opt(ts)), color=PLT_COLORS[1], label="full optimization")
    #now use an envelope times a cosine
    fit_pulse = np.real( psig_opt.get_tenv(ts, field='A') )
    #axs.fill_between(ts, np.zeros(ts.shape), fit_pulse, color=PLT_COLORS[1], alpha=0.2)
    fit_pulse *= np.sqrt(2)*np.sin(2*np.pi*psig_opt.f0*(ts - psig_opt.t0) + psig_opt.phi)
    axs.plot(ts[1:], np.diff(fit_pulse), color=PLT_COLORS[2])
    axs.set_xlim([ts[0], ts[-1]])
    axs.set_ylim([-1, 1])
    axs.legend()
    fig.savefig(fname, bbox_inches='tight')

def plot_before_after_phase(ts, vs, fname):
    freqs = fft.rfftfreq(len(ts), d=ts[1]-ts[0])
    #get signals
    psig_unopt = phases.signal(ts, vs, skip_opt=True)
    psig_opt = phases.signal(ts, vs, skip_opt=False)
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

pf = phases.phase_finder(args.fname, prefix=args.prefix, pass_alpha=args.lowpass, keep_n=-1)
if args.point == '':
    make_heatmaps(pf, n_groups=args.n_groups)
else:
    #phases.signal._do_grad_tests( np.array([0.4]) )

    point_arr = args.point.split(",")
    clust = "cluster_"+point_arr[0]
    j = int(point_arr[1])
    fig_name = "{}/fit_{}_{}".format(pf.prefix,clust,j)
    #read the data and set up the signal analyzer
    v_pts, err_2 = pf.get_point_times(clust, j, low_pass=False)
    psig_before, psig_after = plot_before_after_phase(pf.t_pts, v_pts, "{}_param_est.svg".format(fig_name))
    #plot frequency space
    fig,axs = plt.subplots(2,2)
    plot_raw_fdom(psig_before, axs[:,0], xlim=[0,1])
    plot_raw_fdom(psig_after, axs[:,1], xlim=[0,1], ylabels=False)
    fig.savefig("{}_raw_fdom.svg".format(fig_name), bbox_inches='tight')
    #plot time space
    plot_raw_tdom(pf.t_pts, psig_before, psig_after, "{}_raw_tdom.svg".format(fig_name))
    #plot the frequency domain envelope information
    '''fig, axs = plt.subplots(2)
    psig.compare_fspace(axs[0], sym_mode=args.sym_mode, plt_raxs=not args.plot_y_labels)
    axs[0].set_xlim(-0.05, 0.05)
    psig.compare_tspace(axs[1])
    config_axs(axs[0], "f (1/fs)", "$E(ω)$", leg_loc="upper left")
    config_axs(axs[1], "t (fs)", "$E(t)$", leg_ncol=2)
    #fig.set_size_inches(3.5, 2)
    fig.savefig("{}_fdom.svg".format(fig_name), bbox_inches='tight')
    #plot the time domain envelope information
    fig, axs = plt.subplots()
    psig.compare_tspace(axs)
    #config_axs(axs, "$t$ (fs)", "$E(t)$", leg_ncol=2, w=2, h=1.5)
    fig.savefig("{}_tdom.svg".format(fig_name), bbox_inches='tight')'''
