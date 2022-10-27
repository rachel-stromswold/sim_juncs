import phases
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
import time

plt.rc('font', size=14)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

N_COLS = 1

AMP_RANGE = (0, 0.5)
SIG_RANGE = (0, 10.0)
PHI_RANGE = (-1.0, 1.0)

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Fit time a_pts data to a gaussian pulse envelope to perform CEP estimation.')
parser.add_argument('--fname', type=str, help='h5 file to read', default='field_samples.h5')
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
parser.add_argument('--slice-dir', type=str, help='prefix to use when opening files', default='x')
parser.add_argument('--regression', action='store_true', help='If set to true, perform linear regression on phases across junction', default=False)
args = parser.parse_args()

#decide whether slices are along the x or z axis
slice_dir = args.slice_dir
slice_name = 'z'
if slice_dir == 'x':
    slice_name = 'z'

with h5py.File(args.fname, "r") as f:
    clust_names = list(f.keys())[:-2]
    fig_amp, axs_amp = plt.subplots(len(clust_names)//N_COLS, N_COLS)
    fig_time, axs_time = plt.subplots(len(clust_names)//N_COLS, N_COLS)
    fig_cep, axs_cep = plt.subplots(len(clust_names)//N_COLS, N_COLS)
    #create a numpy array for time points, this will be shared across all points
    keylist = list(f['cluster_0'].keys())
    n_t_pts = len(f['cluster_0'][keylist[1]]['time'])
    n_post_skip = n_t_pts // phases.SKIP
    #read information to figure out time units and step sizes
    t_min = f['info']['time_bounds'][0]
    t_max = f['info']['time_bounds'][1]
    dt = (t_max-t_min)/n_t_pts
    t_pts = np.linspace(t_min, t_max, num=n_t_pts)
    #iterate over each cluster (x slice)
    t_start = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    n_evals = 0

    #DEBUGGING =====================================================================
    '''test_i = 4
    test_j = 9
    tmp_pts = list(f[clust_names[test_i]].keys())
    tmp_v_pts = np.array(f[clust_names[test_i]][tmp_pts[test_j]]['time']['Re'])
    res,res_env,est_omega,est_phi = phases.opt_pulse_env(t_pts, tmp_v_pts, keep_n=2.5, fig_name="{}/fit_figs/blah".format(args.prefix))
    amp, t_0, sig, omega, cep = phases.get_params(res)
    print("cep=%f cep/pi=%f" % (cep, cep/np.pi))'''

    x_0 = (9.0-f[clust_names[0]]['locations'][slice_name][0])/16

    for i, clust in enumerate(clust_names):
        x = (9.0-f[clust]['locations'][slice_name][0])/16 - x_0
        #fetch a list of points and their associated coordinates
        points = list(f[clust].keys())[1:]
        zs = (np.array(f[clust]['locations'][slice_dir])-9)/24
        x_range = (1.05*zs[0], -1.05*zs[0])
        l_gold = [x_range[0], zs[5]]
        r_gold = [-zs[5], x_range[-1]]

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
            v_pts = np.array(f[clust][pt]['time']['Re'])
            # (dt*|dE/dt|)^2 evaluated at t=t_max (one factor of two comes from the imaginary part, the other from the gaussian. We do something incredibly hacky to account for the spatial part. We assume that it has an error identical to the time error, and add the two in quadrature hence the other factor of two.
            #err_2 = 8*np.max(np.diff(v_pts)**2)
            err_2 = 0.02
            print("i={}, j={}\n\tfield error^2={}".format(i, j, err_2))
            res,res_env = phases.opt_pulse_env(t_pts, v_pts, a_sigmas_sq=err_2, keep_n=2.5, fig_name="{}/fit_figs/fit_{}_{}".format(args.prefix,i,j))
            n_evals += 1
            #res_notrunc,res_env_notrunc,_,_ = phases.opt_pulse_env(t_pts, v_pts, keep_n=-1)
            print("\tsquare errors = {}, {}\n\tx={}\n\tdiag(H^-1)={}".format(res.fun, res_env.fun, res.x, np.diagonal(res.hess_inv)))
            #only include this point if the fit was good
            if res.fun < 500.0:
                good_zs.append(zs[j])
                jj = len(good_zs)-1
                amp, t_0, sig, omega, cep = phases.get_params(res)
                '''rel_phase = phase_cor + cep/np.pi
                rel_phase_est = phase_cor + cep/np.pi
                #if the point is closer when translated by 2pi, translate it
                if jj > 0:
                    if np.abs(rel_phase - phase_arr[0,jj-1]) > np.abs(rel_phase+2-phase_arr[0,jj-1]):
                        rel_phase += 2
                        phase_cor += 2
                    elif np.abs(rel_phase - phase_arr[0,jj-1]) > np.abs(rel_phase-2-phase_arr[0,jj-1]):
                        rel_phase -= 2
                        phase_cor -= 2'''
                amp_arr[0,jj] = amp
                amp_arr[1,jj] = np.sqrt(res_env.hess_inv[0][0]/(err_2*n_post_skip))/2
                t_0_arr[0,jj] = t_0
                t_0_arr[1,jj] = np.sqrt(res.hess_inv[1][1]/(err_2*n_post_skip))
                sig_arr[0,jj] = sig
                sig_arr[1,jj] = np.sqrt(res.hess_inv[2][2]/(err_2*8*sig*n_post_skip))
                phase_arr[0,jj] = cep/np.pi
                phase_arr[1,jj] = np.sqrt(res.hess_inv[4][4]/(err_2*np.pi*n_post_skip))
            else:
                print("bad fit! i={}, j={}".format(i,j))

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
            d_ax1.scatter(new_good_zs, phase_arr[0], color='red')
            d_ax1.plot([new_good_zs[0], new_good_zs[-1]], [0, 0], color='black', linestyle=':')
            d_ax1.plot([new_good_zs[0], new_good_zs[-1]], [0.5, 0.5], color='gray', linestyle=':')
            d_ax1.plot([new_good_zs[0], new_good_zs[-1]], [-0.5, -0.5], color='gray', linestyle=':')
            d_ax1.set_xlim(x_range)
            d_ax1.fill_between(l_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
            d_ax1.fill_between(r_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
            d_ax1.set_ylim(PHI_RANGE)
            d_ax1.set_yticks([-1, -0.5, 0, 0.5, 1])
            #amplitude plot
            d_ax2.scatter(new_good_zs, amp_arr[0], color='blue')
            d_ax2.set_xlim(x_range)
            d_ax2.set_ylim(AMP_RANGE)
            d_fig.savefig(args.prefix+"/demo.pdf")
        #matplotlib is annoying and the axes it gives have a different type depending on the column
        if N_COLS == 1:
            tmp_axs_amp = axs_amp[i]
            tmp_axs_time = axs_time[i]
            tmp_axs_cep = axs_cep[i]
        else:
            tmp_axs_amp = axs_amp[i//N_COLS,i%N_COLS]
            tmp_axs_time = axs_time[i//N_COLS,i%N_COLS]
            tmp_axs_cep = axs_cep[i//N_COLS,i%N_COLS]
        #make amplitude plot
        tmp_axs_amp.errorbar(new_good_zs, amp_arr[0], yerr=amp_arr[1], fmt='.', linestyle='')
        #tmp_axs_amp.annotate(r"${0}={1:.2g} \mu m$".format(slice_name, x), (0.01, 1.05), xycoords='axes fraction')
        tmp_axs_amp.set_xlim(x_range)
        tmp_axs_amp.fill_between(l_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_amp.fill_between(r_gold, AMP_RANGE[0], AMP_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_amp.set_ylim(AMP_RANGE)
        #make spread plot
        tmp_axs_time.errorbar(new_good_zs, sig_arr[0], yerr=sig_arr[1], label=r"$\sigma$", fmt='.', linestyle='')
        #tmp_axs_time.annotate(r"${0}={1:.2g} \mu m$".format(slice_name, x), (0.01, 1.05), xycoords='axes fraction')
        tmp_axs_time.set_xlim(x_range)
        tmp_axs_time.fill_between(l_gold, SIG_RANGE[0], SIG_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_time.fill_between(r_gold, SIG_RANGE[0], SIG_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_time.set_ylim(SIG_RANGE)
        #make phase plot
        #tmp_axs_cep.annotate(r"${0}={1:.2g}$ $\mu$m".format(slice_name, x), (0.01, 1.05), xycoords='axes fraction')
        tmp_axs_cep.errorbar(new_good_zs, phase_arr[0], yerr=phase_arr[1], fmt='.', linestyle='')
        tmp_axs_cep.plot([new_good_zs[0], new_good_zs[-1]], [0, 0], color='black', linestyle=':')
        tmp_axs_cep.plot([new_good_zs[0], new_good_zs[-1]], [0.5, 0.5], color='gray', linestyle=':')
        tmp_axs_cep.plot([new_good_zs[0], new_good_zs[-1]], [-0.5, -0.5], color='gray', linestyle=':')
        tmp_axs_cep.set_xlim(x_range)
        tmp_axs_cep.fill_between(l_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_cep.fill_between(r_gold, PHI_RANGE[0], PHI_RANGE[1], color='yellow', alpha=0.3)
        tmp_axs_cep.set_ylim(PHI_RANGE)
        tmp_axs_cep.set_yticks([-1, 1])
        #perform regression if desired
        if args.regression:
            if len(good_zs) > 1:
                cep_line = linregress(new_good_zs, phase_arr[0])
            else:
                cep_line = linregress([0,1],[0,0])#give bad dummy values
            print("x={}: phi(z)={}z+{}, phi(0.1)={}".format(x, cep_line.slope, cep_line.intercept, cep_line.slope*0.1+cep_line.intercept))
            tmp_axs_cep.plot(new_good_zs, cep_line.slope*np.array(new_good_zs) + cep_line.intercept, color='gray')
            tmp_axs_cep.annotate(r"$\phi={0:.2g}z+{1:.2g}$, $r^2={2:.2g}$".format(cep_line.slope, cep_line.intercept, cep_line.rvalue*cep_line.rvalue), (0.01, 0.84), xycoords='axes fraction')
        #hide x axis labels on everything
        if i < len(clust_names)-N_COLS:
            tmp_axs_amp.get_xaxis().set_visible(False)
            tmp_axs_time.get_xaxis().set_visible(False)
            tmp_axs_cep.get_xaxis().set_visible(False)
        if i % N_COLS > 0:
            tmp_axs_amp.get_yaxis().set_visible(False)
            tmp_axs_time.get_yaxis().set_visible(False)
            tmp_axs_cep.get_yaxis().set_visible(False)
    t_dif = time.clock_gettime_ns(time.CLOCK_MONOTONIC) - t_start
    print("Completed optimizations in {:.5E} ns, average time per eval: {:.5E} ns".format(t_dif, t_dif/n_evals))
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
