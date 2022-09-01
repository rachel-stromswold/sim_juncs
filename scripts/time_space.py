import utils
import configparser
import numpy as np
import argparse
import h5py
import scipy.optimize as opt
from scipy.fft import rfft, rfftfreq

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

EPSILON = 0.1
FIELD_LABEL = "ex-"

DAT_ERR_SQ = 0.1

#parse arguments supplied via command line
parser = argparse.ArgumentParser(description='Query sky-surveys for redshift data corresponding to a gravitational-wave detection.')
parser.add_argument('--prefix', type=str, help='prefix to use when opening files', default='.')
parser.add_argument('--gap-width', type=float, help='The width of the junction gap', default=-1.0)
parser.add_argument('--gap-thick', type=float, help='The thickness of the junction', default=-1.0)
parser.add_argument('--pulse-center', type=float, help='center time of the Gaussian pulse', default=0.0)
parser.add_argument('--pulse-width', type=float, help='width the Gaussian pulse', default=1.0)
parser.add_argument('--pulse-frequency', type=float, help='frequency the Gaussian pulse', default=1.0)
parser.add_argument('--pulse-z', type=float, help='z location of the start of the Gaussian beam', default=1.0)
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

#convert args into more convenient shorthands
w_0 = args.pulse_frequency/geom.um_scale
sig = args.pulse_width / w_0
t_0 = args.pulse_center / w_0
w_0 *= 2*np.pi

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
for i in range(n_z_slices):
    axs_td[i][0].set_ylabel(r"$z={0:.2g} \mu m$".format(z_slices[i]))
    axs_fd[i][0].set_ylabel(r"$z={0:.2g} \mu m$".format(z_slices[i]))
    for j in range(1, n_x_slices):
        axs_td[i][j].get_yaxis().set_visible(False)
        axs_fd[i][j].get_yaxis().set_visible(False)
for i in range(n_x_slices):
    axs_td[-1][i].set_xlabel(r"$x={0:.2g} \mu m$".format(x_slices[i]))
    axs_fd[-1][i].set_xlabel(r"$x={0:.2g} \mu m$".format(x_slices[i]))
    for j in range(n_x_slices-1):
        axs_td[j][i].get_xaxis().set_visible(False)
        axs_fd[j][i].get_xaxis().set_visible(False)

#take the fourier transforms for each point we consider
sample_fours = [[rfft(sample_fields[j][k]) for k in range(n_x_slices)] for j in range(n_z_slices)]

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

#get the negative log likelihood of the params=(amplitude, phase) given the dataset d_comp and the vacuum field components v_comp
def phase_likelihood(params, d_comp, v_comp):
    terms = ( np.real(d_comp - params[0]*v_comp*np.exp(1j*params[1])) )**2
    #return -np.log( np.sum(np.exp(-terms/DAT_ERR_SQ)) )
    return np.sum(terms/DAT_ERR_SQ)

#make plots
time_pts = np.linspace(0, geom.meep_time_to_fs(last_time), num=len(sample_fields[0][0]))
freq_cutoff = time_pts.shape[0]//8
freq_pts = rfftfreq(time_pts.shape[0])
for j, z in enumerate(z_slices):
    for k in range(n_x_slices):
        four = rfft(sample_fields[j][k])
        axs_td[j][k].plot(time_pts, sample_fields[j][k])
        axs_td[j][k].set_ylim((-max_t, max_t))
        axs_fd[j][k].plot(freq_pts[:freq_cutoff], np.abs(sample_fours[j][k][:freq_cutoff]), color='black')
        axs_fd[j][k].plot(freq_pts[:freq_cutoff], np.angle(sample_fours[j][k][:freq_cutoff])*angle_scale, linestyle=':', color='red')
        axs_fd[j][k].set_ylim((-max_f, max_f))

        #figure out the posterior on the phase shift inside the material
        s_r = abs(z - args.pulse_z)
        #the fourier transform of the pulse, note that the amplitude scaling is a nuisance parameter so we only care that this is correct up to a constant scaling
        pulse_comps = np.exp(-0.5*(sig*(freq_pts-w_0))**2)*np.exp(1j*(w_0*t_0 - freq_pts*(t_0+s_r)))
        #start with a point with an amplitude 1.0 and no phase shift
        x0 = np.array([1., 0.])
        res = opt.minimize(phase_likelihood, x0, args=(sample_fours[j][k], pulse_comps))
        least_sq = res.x
        print("j={}, k={}: {}".format(j, k, least_sq))
        print(res)
        axs_td[j][k].text(0.01, 0.99, r"$\phi = {}$".format(least_sq[1]))
        axs_fd[j][k].text(0.01, 0.99, r"$\phi = {}$".format(least_sq[1]))
    #draw the frequency color plots
    fours_mid = np.transpose(np.array([rfft(mid_fields[j, :, i]) for i in range(mid_fields.shape[-1])]))
    #convert into units of seconds and figure out the scale of the frequencies
    four_freqs = rfftfreq(mid_fields.shape[1], last_time/(n_t_pts*geom.um_scale))
    top_ind = fours_mid.shape[0]
    #use even columns for real, odd for imaginary
    axs_col[2*j].imshow(np.abs(fours_mid[:top_ind//2, :]))
    axs_col[2*j+1].imshow(np.angle(fours_mid[:top_ind//2, :]), cmap='twilight_shifted')
    axs_col[j].set_xlabel(r'$x$ (um)')

axs_col[0].set_ylabel(r'frequency ($1/\lambda$)')

fig_td.set_tight_layout(True)
fig_fd.set_tight_layout(True)
fig_col.set_tight_layout(True)
plt.ylabel(r"$E_x$ (meep units)")
plt.xlabel(r"$t$ (fs)")
fig_td.savefig(args.prefix+"/tdom_plot.pdf")
plt.xlabel(r"$f$ (2*pi/fs)")
fig_fd.savefig(args.prefix+"/fdom_plot.pdf")
fig_col.savefig(args.prefix+"/fdom_plot_2d.pdf")
