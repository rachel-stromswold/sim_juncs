import phases
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
import time

DEBUG_FITS = False

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Fit time a_pts data to a gaussian pulse envelope to perform CEP estimation.')
parser.add_argument('--fname', type=str, help='h5 file to read', default='field_samples.h5')
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
args = parser.parse_args()

n_fits = 1
with h5py.File(args.fname, "r") as f:
    clust_names = list(f.keys())[:-2]
    fig_amp, axs_amp = plt.subplots(len(clust_names), 1)
    fig_time, axs_time = plt.subplots(len(clust_names), 1)
    fig_cep, axs_cep = plt.subplots(len(clust_names), 1)
    #create a numpy array for time points, this will be shared across all points
    n_t_pts = len(f['cluster_0']['point_00']['time'])
    #read information to figure out time units and step sizes
    t_min = f['info']['time_bounds'][0]
    t_max = f['info']['time_bounds'][1]
    dt = (t_max-t_min)/n_t_pts
    t_pts = np.linspace(t_min, t_max, num=n_t_pts)
    #iterate over each cluster (x slice)
    t_start = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    n_evals = 0
    for i, clust in enumerate(clust_names):
        x = (9.0-f[clust]['locations']['x'][0])/16
        #fetch a list of points and their associated coordinates
        points = list(f[clust].keys())[1:]
        zs = (np.array(f[clust]['locations']['z'])-7.4)/16
        n_z_pts = len(points)
        if zs.shape[0] != n_z_pts:
            raise ValueError("points and locations do not have the same size")

        amp_arr = np.zeros((2,n_z_pts))
        t_0_arr = np.zeros((2,n_z_pts))
        sig_arr = np.zeros((2,n_z_pts))
        phase_arr = np.zeros((3,n_z_pts))
        good_zs = []
        for j, pt in enumerate(points):
            v_pts = np.array(f[clust][pt]['time']['Re'])
            #err_2 = 2*np.sum( np.diff(v_pts)**2 ) / (n_t_pts-1) #=avg( ((dv/dt)*(time step))^2 )
            err_2 = 0.1
            print("field error=%f" % err_2)
            res,res_env,est_omega,est_phi = phases.opt_pulse_env(t_pts, v_pts, keep_n=2)
            n_evals += 1
            if DEBUG_FITS:
                res_notrunc,res_env_notrunc,_,_ = phases.opt_pulse_env(t_pts, v_pts, keep_n=-1)
                fig = plt.figure()
                ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
                ax.plot(t_pts, v_pts, color='black')
                ax.plot(t_pts, phases.gauss_env(res_env.x, t_pts), linestyle=':', color='blue')
                ax.plot(t_pts, phases.gauss_env(res_env_notrunc.x, t_pts), linestyle=':', color='red')
                ax.plot(t_pts, phases.gauss_series(res.x, t_pts), color='blue')
                ax.plot(t_pts, phases.gauss_series(res_notrunc.x, t_pts), color='red')
                fig.savefig("{}/fit_{}.pdf".format(args.prefix, n_fits))
                n_fits += 1
            print(res)
            print(res_env)
            #only include this point if the fit was good
            if res.x[0] > 1e-2:
                good_zs.append(zs[j])
                jj = len(good_zs)-1
                amp, t_0, sig, omega, cep = phases.get_params(res)
                amp_arr[0,jj] = amp
                amp_arr[1,jj] = np.sqrt(res.hess_inv[0][0]/err_2)
                t_0_arr[0,jj] = t_0
                t_0_arr[1,jj] = np.sqrt(res.hess_inv[1][1]/err_2)
                sig_arr[0,jj] = sig
                sig_arr[1,jj] = np.sqrt(res.hess_inv[2][2]/(err_2*8*sig))
                phase_arr[0,jj] = cep/np.pi
                phase_arr[1,jj] = est_phi/np.pi
                phase_arr[2,jj] = np.sqrt(res.hess_inv[4][4]/err_2)
        max_i = len(good_zs)
        axs_amp[i].scatter(good_zs, amp_arr[0][:max_i])
        axs_amp[i].annotate(r"$x={0:.2g} \mu m$".format(x), (0.01, 0.85), xycoords='axes fraction')
        axs_amp[i].set_ylim((0, 1.0))
        axs_time[i].scatter(good_zs, sig_arr[0][:max_i], label=r"$\sigma$")
        axs_time[i].annotate(r"$x={0:.2g} \mu m$".format(x), (0.01, 0.85), xycoords='axes fraction')
        axs_time[i].set_ylim((0, 50))
        cep_line = linregress(good_zs, phase_arr[1][:max_i])
        #axs_cep[i].scatter(good_zs, phase_arr[1][:max_i])
        axs_cep[i].scatter(good_zs, phase_arr[0][:max_i])
        axs_cep[i].plot(good_zs, cep_line.slope*np.array(good_zs) + cep_line.intercept, color='gray')
        axs_cep[i].annotate(r"$x={0:.2g} \mu m$".format(x), (0.01, 0.85), xycoords='axes fraction')
        axs_cep[i].annotate(r"$r^2={0:.4g}$".format(cep_line.rvalue*cep_line.rvalue), (0.01, 0.73), xycoords='axes fraction')
        axs_cep[i].set_ylim((-1, 1))
        #hide x axis labels on everything
        if i < len(clust_names)-1:
            axs_amp[i].get_xaxis().set_visible(False)
            axs_time[i].get_xaxis().set_visible(False)
            axs_cep[i].get_xaxis().set_visible(False)

    t_dif = time.clock_gettime_ns(time.CLOCK_MONOTONIC) - t_start
    print("Completed optimizations in {} ns, average time per eval: {} ns".format(t_dif, t_dif/n_evals))
    fig_amp.suptitle("Amplitude as a function of depth")
    fig_time.suptitle("Pulse width as a function of depth")
    fig_cep.suptitle("Phase as a function of depth")
    axs_amp[1].set_ylabel(r"Amplitude (incident units)")
    axs_time[1].set_ylabel(r"$\sigma$ (fs)")
    axs_cep[1].set_ylabel(r"$\frac{\phi}{\pi}$")
    axs_amp[-1].set_xlabel(r"$z$ ($\mu$m)")
    axs_time[-1].set_xlabel(r"$z$ ($\mu$m)")
    axs_cep[-1].set_xlabel(r"$z$ ($\mu$m)")
    fig_amp.set_tight_layout(True)
    fig_time.set_tight_layout(True)
    fig_cep.set_tight_layout(True)
    fig_cep.savefig(args.prefix+"/phases.pdf")
    fig_amp.savefig(args.prefix+"/amps.pdf")
    fig_time.savefig(args.prefix+"/sigs.pdf")
