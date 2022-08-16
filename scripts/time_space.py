import utils
import configparser
import numpy as np
import argparse
import h5py
from scipy.fft import rfft, rfftfreq

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

EPSILON = 0.1
FIELD_LABEL = "ex-"

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Query sky-surveys for redshift data corresponding to a gravitational-wave detection.')
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
parser.add_argument('--gap-width', type=float, help='The width of the junction gap', default=-1.0)
parser.add_argument('--gap-thick', type=float, help='The thickness of the junction', default=-1.0)
args = parser.parse_args()

#read information from the parameters file

#get a list of all h5 files with the electric field
ex_lst = utils.get_h5_list(FIELD_LABEL[:-1], args.prefix)
last_name =  ex_lst[-1].split('/')[-1]
per_ind = last_name.find('.')
#we expect there to be a decimal in the file name, so find the second instance
last_per_ind = last_name.find('.', per_ind+1)
if last_per_ind > 0:
    per_ind = last_per_ind
last_time = float(last_name[len(FIELD_LABEL):per_ind])

#initialize information classes
geom = utils.Geometry("params.conf", gap_width=args.gap_thick, gap_thick=args.gap_thick)

#information about which points in space we're interested in looking at. Look near the source, just below the interface, and at the center of the junction. For x points, look near the edge (far into the gold), just inside the gold, and at the center
z_slices = [geom.meep_len_to_um(2*geom.pml_thick), geom.meep_len_to_um(geom.t_junc + EPSILON), geom.meep_len_to_um(geom.z_center)]
x_slices = [geom.meep_len_to_um(2*geom.pml_thick), geom.meep_len_to_um(geom.l_junc - EPSILON), geom.meep_len_to_um(geom.z_center)]
x_inds = []
n_z_slices = len(z_slices)
n_x_slices = len(x_slices)
n_t_pts = len(ex_lst)

#store information about the fields as a function of time
sample_fields = [[np.zeros(n_t_pts) for i in range(n_x_slices)] for j in range(n_z_slices)]
mid_fields = [[] for z in range(n_z_slices)]

#iterate over the h5 files to get information in time
for i, fname in enumerate(ex_lst):
    fields,_ = geom.get_cross_section_fields(fname, z_slices)
    #figure out the x indices we want to read from based on the slices array. We check if it has already been initialized so that we don't have to repeat computations
    if len(x_inds) < len(x_slices):
        for x in x_slices:
            x_inds.append( int(geom.um_to_meep_len(x)*fields.shape[1]/geom.tot_len) )
    #add the time point to the sample field arrays
    for j in range(n_z_slices):
        for k, ii in enumerate(x_inds):
            sample_fields[j][k][i] = fields[j][ii]
            mid_fields[j].append(fields[j, :])

mid_fields = np.array(mid_fields)

fig_td, axs_td = plt.subplots(n_z_slices, n_x_slices)
fig_fd, axs_fd = plt.subplots(n_z_slices, n_x_slices)
fig_col, axs_col = plt.subplots(1, 2*n_z_slices)

#set labels for the plots
for j in range(n_z_slices):
    axs_td[j][0].set_ylabel(r"$z={} \mu m$".format(z_slices[j]))
    axs_fd[j][0].set_ylabel(r"$z={} \mu m$".format(z_slices[j]))
for k in range(n_x_slices):
    axs_td[0][k].set_xlabel(r"$x={} \mu m$".format(x_slices[k]))
    axs_fd[0][k].set_xlabel(r"$x={} \mu m$".format(x_slices[k]))

for j in range(n_z_slices):
    for k in range(n_x_slices):
        four = rfft(sample_fields[j][k])
        axs_td[j][k].plot(sample_fields[j][k])
        axs_fd[j][k].plot(rfftfreq(sample_fields[j][k].shape[0]), np.abs(four), color='black')
        axs_fd[j][k].plot(rfftfreq(sample_fields[j][k].shape[0]), np.angle(four), linestyle=':', color='red')
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
    #axs_col[j].axvline((geom.z_center+gap_width/2)/geom.um_scale, 0, 1, color='gray')
    #axs_col[j].axvline((geom.z_center-gap_width/2)/geom.um_scale, 0, 1, color='gray')
    axs_col[j].set_xlabel(r'$x$ (um)')
axs_col[0].set_ylabel(r'frequency ($1/\lambda$)')

fig_td.savefig(args.prefix+"/tdom_plot.pdf")
fig_fd.savefig(args.prefix+"/fdom_plot.pdf")
fig_col.savefig(args.prefix+"/fdom_plot_2d.pdf")
