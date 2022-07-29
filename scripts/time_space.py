import utils
import configparser
import numpy as np
import argparse
import h5py
from scipy.fft import rfft, rfftfreq

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIELD_LABEL = "ex-"

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Query sky-surveys for redshift data corresponding to a gravitational-wave detection.')
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
args = parser.parse_args()

#read information from the parameters file
config = configparser.ConfigParser()
config.read("params.conf")
devs = float(config["pulse_shape"]["pulse_cutoff"])
pulse_width = float(config["pulse_shape"]["pulse_width"])
post_pulse_runtime = float(config["pulse_shape"]["post_pulse_runtime"])
last_time = pulse_width*devs + post_pulse_runtime

#get a list of all h5 files with the electric field
ex_lst = utils.get_h5_list(FIELD_LABEL[:-1], args.prefix)

#initialize information classes
geom = utils.Geometry("params.conf")

#information about which points in space we're interested in looking at
z_slices = [geom.z_center-geom.junc_max_z/2, geom.z_center, geom.z_center+geom.junc_max_z/2]
n_slices = len(z_slices)
n_t_pts = len(ex_lst)

#store information about the fields as a function of time
center_fields = [np.zeros(n_t_pts) for z in z_slices]
mid_fields = [[] for z in range(n_slices)]

#iterate over the h5 files to get information in time
for i, fname in enumerate(ex_lst):
    fields,_ = geom.get_cross_section_fields(fname, z_slices)
    for j in range(n_slices):
        center_fields[j][i] = fields[j][fields.shape[1]//2]
        mid_fields[j].append(fields[j, :])

mid_fields = np.array(mid_fields)

fig_td, axs_td = plt.subplots(n_slices)
fig_fd, axs_fd = plt.subplots(n_slices)
fig_col, axs_col = plt.subplots(1, 2*n_slices)
for j in range(n_slices):
    four = rfft(center_fields[j])
    axs_td[j].plot(center_fields[j])
    axs_fd[j].plot(rfftfreq(center_fields[j].shape[0]), four)
    #draw the frequency color plots
    fours_mid = np.transpose(np.array([rfft(mid_fields[j, :, i]) for i in range(mid_fields.shape[-1])]))
    #convert into units of seconds and figure out the scale of the frequencies
    four_freqs = rfftfreq(mid_fields.shape[1], last_time/(n_t_pts*geom.um_scale))
    top_ind = fours_mid.shape[0]
    #use even columns for real, odd for imaginary
    axs_col[2*j].imshow(np.abs(fours_mid[:top_ind//2, :]))
    axs_col[2*j+1].imshow(np.angle(fours_mid[:top_ind//2, :]), cmap='twilight_shifted')
    #axs_col[j].set_xlim((0, geom.tot_len/geom.um_scale))
    #axs_col[j].set_ylim((four_freqs[0], four_freqs[-1]))
    #axs_col[j].axvline((geom.z_center+middle_w/2)/geom.um_scale, 0, 1, color='gray')
    #axs_col[j].axvline((geom.z_center-middle_w/2)/geom.um_scale, 0, 1, color='gray')
    axs_col[j].set_xlabel(r'$x$ (um)')
axs_col[0].set_ylabel(r'frequency ($1/\lambda$)')

fig_td.savefig(args.prefix+"/tdom_plot.pdf")
fig_fd.savefig(args.prefix+"/fdom_plot.pdf")
fig_col.savefig(args.prefix+"/fdom_plot_2d.pdf")
