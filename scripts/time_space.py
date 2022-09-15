import phases
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

N_COLS = 4

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Fit time a_pts data to a gaussian pulse envelope to perform CEP estimation.')
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
args = parser.parse_args()

with h5py.File(args.prefix+"/field_samples.h5", "r") as f:
    point_names = list(f['cluster_0'].keys())[1:]
    n_pts = len(point_names)
    n_x_slices = N_COLS
    n_z_slices = n_pts // N_COLS
    #figure out the number of time points by inspecting the first time series
    n_t_pts = len(f['cluster_0'][point_names[0]]['time'])
    n_f_pts = len(f['cluster_0'][point_names[0]]['frequency'])
    last_time = f['info']['time_bounds'][1]
    #read information about the source
    pulse_width = f['info']['sources'][0]['width']
    omega_0 = np.pi/(16*f['info']['sources'][0]['wavelen'])
    t_0 = (f['info']['sources'][0]['start_time'] + f['info']['sources'][0]['end_time'])/2
    #store locations for each point
    x_slices = [0.0 for i in range(N_COLS)]
    z_slices = [0.0 for i in range(n_z_slices)]
    sample_fields = [[np.zeros(n_t_pts) for j in range(n_x_slices)] for i in range(n_z_slices)]
    sample_fours = [[np.zeros(n_f_pts) for j in range(n_x_slices)] for i in range(n_z_slices)]
    #figure out the x locations by iterating over the first row
    for i in range(n_x_slices):
        x_slices[i] = f['cluster_0']['locations'][i][0]
    #iterate over each point
    for i in range(n_z_slices):
        #figure out the z locations by iterating over the first column
        z_slices[i] = f['cluster_0']['locations'][N_COLS*i][2]
        for j in range(n_x_slices):
            #iterate over times and take only the real parts
            sample_fields[i][j] = np.array(f['cluster_0'][point_names[i*N_COLS+j]]['time']['Re'])
            sample_fours[i][j] = np.array(f['cluster_0'][point_names[i*N_COLS+j]]['frequency']['Re'] + 1j*f['cluster_0'][point_names[i*N_COLS+j]]['frequency']['Re'])

#find the maximum values for each point we consider
max_t = 0.0
max_f = 0.0
for j in range(n_z_slices):
    for k in range(n_x_slices):
        tmp_max_t = max(np.abs(sample_fields[j][k]))
        tmp_max_f = max(np.abs(sample_fours[j][k]))
        if tmp_max_t > max_t:
            max_t = tmp_max_t
        if tmp_max_f > max_f:
            max_f = tmp_max_f
angle_scale = max_f/np.pi

#generate the time points
print("last_time = %f, %f" % (last_time, last_time))
time_pts = np.linspace(0, last_time, num=len(sample_fields[0][0]))

#get phase estimates
reses = [[phases.opt_pulse_env(time_pts, sample_fields[i][j]) for j in range(n_x_slices)] for i in range(n_z_slices)]

#make plots
fig_td, axs_td = plt.subplots(n_z_slices, n_x_slices)
#set labels for the plots
for i in range(n_z_slices):
    axs_td[i][0].set_ylabel(r"$z={0:.2g} \mu m$".format(z_slices[i]))
    for j in range(1, n_x_slices):
        axs_td[i][j].get_yaxis().set_visible(False)
for i in range(n_x_slices):
    axs_td[-1][i].set_xlabel(r"$x={0:.2g}$ fs".format(x_slices[i]))
    for j in range(n_x_slices-1):
        axs_td[j][i].get_xaxis().set_visible(False)

t_bounds = (-max_t*1.5, max_t*1.5)
for j, z in enumerate(z_slices):
    for k in range(n_x_slices):
        res_opt = reses[j][k]
        axs_td[j][k].plot(time_pts, sample_fields[j][k])
        axs_td[j][k].plot(time_pts, phases.gauss_env(res_opt[1].x, time_pts), linestyle=':', color='gray')
        axs_td[j][k].plot(time_pts, phases.gauss_series(res_opt[0].x, time_pts))
        phase = res_opt[0].x[4]
        if abs(phase) > 0.1:
            axs_td[j][k].annotate(r"$\phi/\pi = {0:.2g}$".format(phase/np.pi), (0.4, 0.1), xycoords='axes fraction')
        axs_td[j][k].set_ylim(t_bounds)

fig_td.set_tight_layout(True)
plt.ylabel(r"$E_x$ (meep units)")
plt.xlabel(r"$t$ (fs)")
fig_td.savefig(args.prefix+"/tdom_plot.pdf")
