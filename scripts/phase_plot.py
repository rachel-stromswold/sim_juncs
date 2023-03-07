import phases
import utils
import argparse
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import linregress

plt.rc('font', size=14)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

N_COLS = 1

AMP_RANGE = (0, 1.5)
OMG_RANGE = (0, 4.0)
SIG_RANGE = (0, 10.0)
PHI_RANGE = (-1.0, 1.0)
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
parser.add_argument('--save-fit-figs', action='store_true', help='If set, then intermediate plots of fitness are saved to <prefix>/fit_figs where <prefix> is specified by the --prefix flag.', default=False)
parser.add_argument('--gap-width', type=float, help='junction width', default=0.1)
parser.add_argument('--gap-thick', type=float, help='junction thickness', default=0.2)
parser.add_argument('--diel-const', type=float, help='dielectric constant of material', default=3.5)
parser.add_argument('--lowpass', type=float, help='If optimizations fail, then the phase finder tries them again after applying a lowpass filter to the timeseries. This parameter specifies the strength of the lowpass filter. If this script is crashing, then try increasing this parameter.', default=1.0)
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
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

    def get_loc_rel(self, z, ref_z=-1):
        '''return the location of the cluster with the name clust along the propagation direction (z if slice_dir=x x if slice_dir=z)'''
        if ref_z < 0:
            ref_z = self.geom.t_junc

    def convert_x_units(self, xs):
        '''given an np array of x data points, return that array relative to the center and in micrometers instead of meep units
        '''
        return (xs - self.geom.z_center)/self.geom.um_scale

    def reflect_pts(self, cr):
        '''mirror the cluster_result object cr and label the x axis appropriately
        '''
        return cr.mirror_pts(self.geom.z_center, self.geom.um_scale)

    def get_junc_bounds(self):
        return self.geom.meep_len_to_um(self.geom.l_junc - self.geom.z_center), self.geom.meep_len_to_um(self.geom.r_junc - self.geom.z_center)

    def get_field_amps(self, x_pts, z, diel_const, n_modes=HIGHEST_MODE):
        length = self.geom.meep_len_to_um(self.geom.r_junc - self.geom.l_junc)
        #these are some useful constants that we define. Note that lc is just the wavelength inside the dielectric, kc_sq is the square of 2pi/lc and l_rat_sq is (lc/2l)^2 where l is the length of the sample. We derived the expression beta^2 = kc^2 (1 - l_rat_sq*n^2), so neg_beta_sq is just -beta^2.
        #lc = vac_wavelen/np.sqrt(np.real(diel_const))
        lc = self.vac_wavelen/np.real(diel_const)
        kc_sq = (2*np.pi / lc)**2
        l_rat_sq = 0.25*(lc/length)**2
        fields = 1j*np.zeros(len(x_pts))
        print("length={}, L_c={}, kc^2={}, ratio^2={}".format(length, lc, kc_sq, l_rat_sq))
        phaseless_comps = self.freq_comps*np.exp(1j*self.vac_phase)
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
        return 4*0.5*self.vac_amp*np.real(fields)/np.pi

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

    #def fit_eps(self, x_pts, amps, z):
    def fit_eps(self, pf, inds=[]):
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
            cr = self.reflect_pts(pf.lookup_fits(clust, recompute=True))
            amp_arr = cr.get_amp()
            #the TE mode expansion is only valid inside the junction, so we need to truncate the arrays
            lowest_valid = -1
            highest_valid = -1
            slice_zs.append(self.get_clust_location(clust))
            '''high_off = (junc_bounds[1]-junc_bounds[0])/(2*HIGHEST_MODE-1)
            fit_bounds = (junc_bounds[0]+high_off, junc_bounds[1]-high_off)'''
            fit_bounds = junc_bounds
            for j, x in enumerate(cr.xs):
                if lowest_valid < 0 and x >= fit_bounds[0]:
                    lowest_valid = j
                if x <= fit_bounds[1]:
                    highest_valid = j
            slice_xs.append(cr.xs[lowest_valid:highest_valid+1])
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

def make_fits(pf, n_groups=-1, recompute=False):
    _,_,axs_mapping = find_grouping(pf, n_groups)
    #use a default set of axes if the user didn't supply one or the supplied axes were invalid
    if axs_mapping is None:
        axs_mapping = range(len(pf.clust_names))
    n_mapped = max(axs_mapping)+1
    if n_mapped > len(pf.clust_names):
        axs_mapping = range(len(pf.clust_names))
        n_mapped = pf.n_clusts

    wg = waveguide_est(args.gap_width, args.gap_thick, pf)
    junc_bounds = wg.get_junc_bounds()

    #initialize plots
    fig_amp, axs_amp = plt.subplots(n_mapped//N_COLS, N_COLS)
    fig_phs, axs_phs = plt.subplots(n_mapped//N_COLS, N_COLS)
    fig_omg, axs_omg = plt.subplots(n_mapped//N_COLS, N_COLS)
    #set up axes first
    for i in range(n_mapped):
        label_x = (i == len(pf.clust_names) - 1)
        wg.setup_axes(get_axis(axs_amp, i), AMP_RANGE, label_x=label_x)
        wg.setup_axes(get_axis(axs_phs, i), PHI_RANGE, label_x=label_x)
        wg.setup_axes(get_axis(axs_omg, i), OMG_RANGE, label_x=label_x)

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
        cr = wg.reflect_pts(pf.lookup_fits(clust, recompute=recompute, save_fit_figs=args.save_fit_figs))
        #save x points, amplitudes and phases so that an average may be computed
        #figure out expected amplitudes from waveguide modes
        clust_z = pf.get_clust_location(clust)
        amps_fit = np.abs(wg.get_field_amps(x_cnt, clust_z, args.diel_const))
        phss_fit = np.abs(wg.get_field_phases(x_cnt, clust_z, args.diel_const, -np.pi/2))
        #plot amplitudes
        tmp_axs = get_axis(axs_amp, i)
        tmp_axs.plot(x_cnt, amps_fit, color='gray', linestyle=':')
        tmp_axs.scatter(cr.xs, cr.get_amp(), s=3)
        tmp_axs.annotate(r"$z={}\mu$m".format(round(clust_z, 3)), (0.01, 0.68), xycoords='axes fraction')
        #plot phases
        tmp_axs = get_axis(axs_phs, i)
        #tmp_axs.plot(x_cnt, phss_fit)
        phase_th = -pf.get_src_phase()[0]/np.pi
        tmp_axs.plot([cr.xs[0], cr.xs[-1]], [phase_th, phase_th], color='gray', linestyle=':')
        tmp_axs.scatter(cr.xs, cr.get_phase(), s=3)
        #plot frequencies
        tmp_axs = get_axis(axs_omg, i)
        tmp_axs.annotate(r"$z={}\mu$m".format(round(clust_z, 3)), (0.01, 0.68), xycoords='axes fraction')
        omega_th = 2*np.pi*.299792458/pf.get_src_wavelen()[0] #2*pi*c/lambda
        tmp_axs.plot([x_cnt[0], x_cnt[-1]], [omega_th, omega_th], color='gray', linestyle=':')
        tmp_axs.scatter(cr.xs, cr.get_omega(), s=3)

    #save the figures
    fig_amp.savefig(args.prefix+"/amps.pdf")
    fig_phs.savefig(args.prefix+"/phases.pdf")
    fig_omg.savefig(args.prefix+"/sigs.pdf")
    fig_fits.savefig(args.prefix+"/fits.pdf")

def make_heatmap(fg, ax, imdat, title, label, rng=None, cmap='viridis'):
    ax[0].set_title(title)
    ax[0].get_xaxis().set_visible(False)
    ax[0].get_yaxis().set_visible(False)
    #do plots
    if rng is None:
        im = ax[0].imshow(imdat, cmap=cmap)
    else:
        im = ax[0].imshow(imdat, vmin=rng[0], vmax=rng[1], cmap=cmap)
    cbar = fg.colorbar(im, cax=ax[1])
    cbar.ax.set_ylabel(label)
    fg.tight_layout()

def plot_average_phase(pf, n_groups=-1):
    grp_len,n_groups,_ = find_grouping(pf, n_groups)
    cl_xs = []
    cl_amp = []
    cl_phs = []
    cl_ampr = []
    cl_err = []
    wg = waveguide_est(args.gap_width, args.gap_thick, pf)
    for clust in pf.clust_names:
        cr = wg.reflect_pts(pf.lookup_fits(clust, recompute=False))
        #save x points, amplitudes and phases so that an average may be computed
        cl_xs.append(cr.xs)
        cl_amp.append(cr.get_amp())
        cl_phs.append(cr.get_phase())
        cl_ampr.append(cr.get_amp_ref())
        cl_err.append(cr.get_err_sq())
    cl_xs = np.array(cl_xs)
    cl_amp = np.array(cl_amp)
    cl_phs = np.array(cl_phs)
    cl_ampr = np.array(cl_ampr)
    cl_err = np.array(cl_err)

    #figure out the extent in the x and z directions
    x_min = cl_xs[0][0]
    x_max = cl_xs[0][-1]
    z_min = pf.get_clust_location(pf.clust_names[0])
    z_max = pf.get_clust_location(pf.clust_names[-1])

    #phs_colors = [(0.18,0.22,0.07), (0.36,0.22,0.61), (1,1,1), (0.74,0.22,0.31), (0.28,0.32,0.17)]
    phs_colors = [(0.18,0.22,0.07), (0.36,0.22,0.61), (1,1,1), (0.74,0.22,0.31), (0.18,0.22,0.07)]
    cmap = LinearSegmentedColormap.from_list('phase_colors', phs_colors, N=100)

    #save heatmaps of amplitudes and phases
    for i in range(n_groups):
        heat_fig, heat_ax = plt.subplots(2,2, gridspec_kw={'width_ratios':[24,1]})
        make_heatmap(heat_fig, heat_ax[0], 2*cl_amp[i*grp_len:(i+1)*grp_len], r"Amplitudes L={} $\mu$m $\lambda$={} $\mu$m".format(args.gap_thick, wg.vac_wavelen), "amplitude (arb. units)", rng=AMP_RANGE, cmap='Reds')
        #make_heatmap(heat_fig, heat_ax[1], cl_phs[i*grp_len:(i+1)*grp_len], "Phases", r"$\phi/2\pi$", rng=PHI_RANGE, cmap='twilight_shifted')
        make_heatmap(heat_fig, heat_ax[1], cl_phs[i*grp_len:(i+1)*grp_len], "Phases", r"$\phi/2\pi$", rng=PHI_RANGE, cmap=cmap)
        heat_fig.savefig(args.prefix+"/heatmap_grp{}.pdf".format(i))
        plt.close(heat_fig)
        #plot errors
        heat_fig, heat_ax = plt.subplots(1,2, gridspec_kw={'width_ratios':[24,1]})
        make_heatmap(heat_fig, heat_ax, cl_err[i*grp_len:(i+1)*grp_len], "Square errors", "fit error (arb. units)", rng=[0,1])
        heat_fig.savefig(args.prefix+"/heatmap_err_grp{}.pdf".format(i))
        plt.close(heat_fig)
        #plot reflected phases
        heat_fig, heat_ax = plt.subplots(2,2, gridspec_kw={'width_ratios':[24,1]})
        make_heatmap(heat_fig, heat_ax[0], cl_phs[i*grp_len:(i+1)*grp_len], "Phases", r"$\phi/2\pi$", rng=PHI_RANGE, cmap='twilight_shifted')
        make_heatmap(heat_fig, heat_ax[1], cl_ampr[i*grp_len:(i+1)*grp_len], "Reflected amplitudes", r"$\phi/2\pi$", rng=AMP_RANGE, cmap='Reds')
        heat_fig.savefig(args.prefix+"/heatmap_phs_grp{}.pdf".format(i))
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
    avg_ax[0].set_ylabel(r"$<\phi>/\pi$")
    avg_ax[1].set_ylabel(r"$<E_0>$ (arb. units)")
    avg_ax[1].set_xlabel(r"x $\mu$m")
    #plot averages
    for i in range(n_groups):
        avg_ax[0].scatter(cl_xs[0], avg_phs[i])
        avg_ax[1].scatter(cl_xs[0], tot_amp[i])
    avg_fig.savefig(args.prefix+"/avgs.pdf")

pf = phases.phase_finder(args.fname, args.gap_width, args.gap_thick, prefix=args.prefix, pass_alpha=args.lowpass)
make_fits(pf, n_groups=args.n_groups, recompute=args.recompute)
plot_average_phase(pf, n_groups=args.n_groups)
