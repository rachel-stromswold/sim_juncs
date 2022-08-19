import numpy as np
import argparse
import h5py
import glob
import utils

import matplotlib

#plotting parameters

#compare the analytic and simulated electric fields at an arbitrary set of time steps
field_range = (-1.0, 1.0)
sq_er_range = (0.0, 0.002)
#the field times we are interested in
#field_times = ["1.", "2.", "3.", "4."]
#field_times = ["16.", "20.", "24.", "28.", "32.", "36."]
field_times = ["2.", "4.", "6.", "8.", "10.", "12."]
N_ROWS = 2
N_COLS = len(field_times)//N_ROWS

EPSILON = 0.1

MOVIE_FIELDS = 1
MOVIE_CROSS_SECS = 2
MOVIE_FRAME_SKIP = 1

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Query sky-surveys for redshift data corresponding to a gravitational-wave detection.')
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
parser.add_argument('--compare-prefix', type=str, help='prefix to use when opening files for comparison', default='')
parser.add_argument('--compare-scaling', type=float, help='Apply a scaling factor to fields obtained from the comparison prefix', default=1)
parser.add_argument('--grid-res', type=float, help='Size of the grid used for simulations', default='20.0')
parser.add_argument('--movie', type=str, help='If set to true, save pngs for an animation of the electic field', default='')
parser.add_argument('--fitting', action='store_true', help='If set to true, generate posterior probabilities for the speed of light', default=False)
parser.add_argument('--ene-flux', action='store_true', help='If set to true, calculate the flux of energies through points', default=False)
parser.add_argument('--show-plots', action='store_true', help='If set to true, show plots instead of saving them', default=False)
parser.add_argument('--gap-width', type=float, help='The width of the junction gap', default=-1.0)
parser.add_argument('--gap-thick', type=float, help='The thickness of the junction', default=-1.0)
args = parser.parse_args()
show_plt = args.show_plots
#disable display if we aren't showing plots
if not show_plt:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 6]

#figure out what type of movies to make
movie_type = 0
if args.movie != '':
    movie_type = MOVIE_CROSS_SECS
    if (args.movie == 'fields'):
        movie_type = MOVIE_FIELDS

#in theory we're working in c=1 units, but the simulation seems to behave better with a slightly lower value
conf_level = 0.95
C_PRIOR_RANGE = (0.999, 1.001)

PHASE_PRIOR_RANGE = (-0.1, 0.1)
light_speed = 1.0
field_label = "ex-"
pml_exclude_fact = 2.0
epsilon = 0.01

#initialize information classes
geom = utils.Geometry("params.conf", gap_width=args.gap_width, gap_thick=args.gap_thick)

#specify which z values we want for cross sections
CROSS_SEC_ZS = [geom.meep_len_to_um(2*geom.pml_thick), geom.meep_len_to_um(geom.t_junc + EPSILON), geom.meep_len_to_um(geom.z_center)]

#plot the fourier components of the Gaussian envelope for debugging purposes
if args.show_plots:
    plt.plot(geom.src.pulse_freqs, np.real(geom.eps_2s))
    plt.plot(geom.src.pulse_freqs, np.imag(geom.eps_2s))
    plt.show()

#list all h5 files and figure out the appropriate time intervals and spatial resolutions
ex_lst = utils.get_h5_list(field_label[:-1], args.prefix)
min_time = float(ex_lst[0][ex_lst[0].find(field_label)+len(field_label):-3])
max_time = float(ex_lst[-1][ex_lst[-1].find(field_label)+len(field_label):-3])
time_step = (max_time - min_time) / (len(ex_lst)-1)
time_lst = np.arange(min_time, max_time+time_step*1.5, time_step)
n_z_pts = 0
z_step = 0.0
with h5py.File(ex_lst[0], "r") as f:
    n_z_pts = f["ex.r"].shape[0]
z_step = geom.tot_len / (n_z_pts - 1)

utils.plot_h5_diel(args.prefix)

#initizlize the posterior distribution for the speed of light
p_pts = np.linspace(PHASE_PRIOR_RANGE[0], PHASE_PRIOR_RANGE[1])
p_posterior = np.ones(p_pts.shape[0])
c_pts = np.linspace(C_PRIOR_RANGE[0], C_PRIOR_RANGE[1])
c_posterior = np.ones(c_pts.shape[0])

#matplotlib doesn't like making columns with more than 5 rows
if geom.n_dims == 1:
    field_times = field_times[1:]
    fig_nd, axs_nd = plt.subplots(len(field_times))
    fig_xs, axs_xs = plt.subplots(2)
    er_fig, er_axs = plt.subplots(len(field_times))
else:
    fig_nd, axs_nd = plt.subplots(N_ROWS, N_COLS)
    fig_xs, axs_xs = plt.subplots(len(field_times), len(CROSS_SEC_ZS))
    er_fig, er_axs = plt.subplots(N_ROWS, N_COLS)
    for i in range(len(field_times)):
        axs_xs[i, 0].set_xlabel("above junction")
        axs_xs[i, 1].set_xlabel("middle junction")
        axs_xs[i, 2].set_xlabel("below junction")
fig_nd.suptitle("Field strength as a function of position (resolution={}, scale ratio={})".format(args.grid_res, geom.um_scale))

sq_err_int_space = 0.0
worst_sq_err_int_space = 0.0
#if the user supplied a comparison prefix, perform the comparison and check for errors
if (args.compare_prefix != ''):
    ex_lst_2 = utils.get_h5_list(field_label[:-1], args.compare_prefix)
    for i, ft in enumerate(field_times):
        fname_1,time = utils.get_closest_time(ft, args.prefix)
        fname_2,_ = utils.get_closest_time(ft, args.compare_prefix)
        field_list_1,_ = utils.get_fields_from_file(fname_1, max_z=geom.tot_len)
        field_list_2,_ = utils.get_fields_from_file(fname_2, max_z=geom.tot_len)
        diff_sq = np.sum((field_list_1 - args.compare_scaling*field_list_2)**2)
        print("Square error at time t={}: {}".format(time, diff_sq))

#generate plots comparing analytic and simulated results
for i, ft in enumerate(field_times):
    #plt.subplot(5, 1, i+1)
    grid_res = 20.0

    fname, time = utils.get_closest_time(ft, args.prefix)
    print("now reading {} using time t={}".format(fname, time))

    #layouts are different for 1d and 2d cases
    cur_axs = None
    cur_er_axs = None
    if geom.n_dims == 1:
        cur_axs = axs_nd[i]
        cur_er_axs = er_axs[i]
    else:
        cur_axs = axs_nd[i//N_COLS, i%N_COLS]
        cur_er_axs = er_axs[i//N_COLS, i%N_COLS]

    _,tmp_err = geom.plot_h5_fields(fname, False, time=time, axs=axs_nd[i//N_COLS, i%N_COLS], er_axs=cur_er_axs)
    geom.plot_cross_section(fname, CROSS_SEC_ZS, axs=axs_xs[i, :])
    for ax in axs_xs[i, :]:
        ax.set_ylabel("t={}".format(ft))
    sq_err_int_space += tmp_err / len(field_times)
    if tmp_err > worst_sq_err_int_space:
        worst_sq_err_int_space = tmp_err
    print("ft={}: err={}".format(ft, tmp_err))
    
if show_plt:
    plt.show()
else:
    fig_nd.tight_layout()
    fig_xs.tight_layout()
    er_fig.tight_layout()
    fig_nd.savefig(args.prefix+"/space_plot.pdf")
    fig_xs.savefig(args.prefix+"/cross_plot.pdf")
    er_fig.savefig(args.prefix+"/error_plot.pdf")
    plt.clf()

#make an animation of the wave if the user asked for one
if movie_type > 0:
    for i, fname in enumerate(ex_lst):
        if i % MOVIE_FRAME_SKIP == 0:
            outname = args.prefix+"/im_{}.png".format(i//MOVIE_FRAME_SKIP)
            #plot_h5_fields(fname, False)
            if movie_type == MOVIE_FIELDS:
                tmp_fig,er = geom.plot_h5_fields(fname, False)
                plt.savefig(outname)
                plt.close(tmp_fig)
            else:
                tmp_fig, tmp_axs = plt.subplots(len(CROSS_SEC_ZS))
                geom.plot_cross_section(fname, CROSS_SEC_ZS, tmp_axs)
                tmp_fig.savefig(outname)
                plt.close(tmp_fig)

if args.fitting:
    #normalize the posterior distribution and find mode and expectation value
    int_fact = (PHASE_PRIOR_RANGE[1]-PHASE_PRIOR_RANGE[0]) / p_pts.shape[0]
    p_posterior = p_posterior / ( np.sum(p_posterior)*int_fact )
    c_max_ind = np.argmax(p_posterior)
    p_ex = np.sum( p_pts*p_posterior )*int_fact
    p_vr = np.sum( p_posterior*(p_pts-p_ex)**2 )*int_fact
    #print useful information about the posterior distribution
    print("maximum posterior probability (phase={}, p(phase|D)={})".format(p_pts[c_max_ind], p_posterior[c_max_ind]))
    print("posterior expectation value phase={}+-{}".format(p_ex, np.sqrt(p_vr)))
    plt.plot(p_pts, p_posterior)
    plt.axvline(p_ex)
    plt.savefig(args.prefix+"/posterior_plot.pdf")
    plt.clf()
