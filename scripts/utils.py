import numpy as np
import configparser
import glob
import h5py

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PML_EXCLUDE_FACT = 2.0
FIELD_RANGE = (-0.5, 0.5)

LIGHT_SPEED = 1.0
PULSE_RES = 0.02
FREQ_RANGE_SCALE = 1.5
BEST_PHASE = -0.00204081632653 #taken from the maximum posterior phase obtained at the highest resolution

class dielectric_func:
    def __init__(self, suscept_str, eps_inf, um_scale=1.0):
        ent_lst = suscept_str.split(')')
        self.omegas = np.zeros(len(ent_lst)-1)
        self.gammas = np.zeros(len(ent_lst)-1)
        self.sigmas = np.zeros(len(ent_lst)-1)
        self.lorentz = np.ones(len(ent_lst)-1)
        self.eps_inf = eps_inf
        print("using susceptibilities: eps_inf={}".format(self.eps_inf))
        #(
        for i, st in enumerate(ent_lst):
            st = st[st.find("(")+1:]#)
            if len(st) > 0:
                st_comma_sep = st.split(',')
                self.omegas[i] = float(st_comma_sep[0])*(2*np.pi)/um_scale
                self.gammas[i] = float(st_comma_sep[1])*(2*np.pi)/um_scale
                self.sigmas[i] = float(st_comma_sep[2])*2*np.pi
                if len(st_comma_sep) > 3:
                    tok = st_comma_sep[3].strip().lower()
                    if tok == 'false' or tok == '0' or tok == 'lorentz':
                        self.lorentz[i] = 1
                    elif tok == 'true' or tok == '1' or tok == 'drude':
                        self.lorentz[i] = 0
                    else:
                        raise ValueError("Unknown susceptibility token: {}".format(tok))
                print("\t(omega, gamma, sigma)=({}, {}, {})".format(self.omegas[i], self.gammas[i], self.sigmas[i]))
        self.numers = self.sigmas*(self.omegas**2)

    def get_eps(self, omega):
        return np.sum( self.eps_inf + self.numers / (self.lorentz*self.omegas**2 - omega**2 - 1j*omega*self.gammas) )

    def get_eps_arr(self, omegas):
        #return np.array( [self.eps_inf + np.sum( self.sigmas*self.omegas**2 / (self.lorentz*self.omegas**2 - omega**2 - 1j*omega*self.gammas) ) for omega in omegas] )
        return np.array( [self.eps_inf + np.sum( self.numers / (self.lorentz*self.omegas**2 - omega**2 - 1j*omega*self.gammas) ) for omega in omegas] )

def get_h5_list(datname, directory):
    '''returns a sorted list of all h5 files matching the specified prefix, sorted alphabetically'''
    h5_lst = sorted(glob.glob(directory+"/"+datname+"*.h5"))
    prefix_len = len(datname)
    #remove all h5 files that aren't of the type we want (x componenet of field)
    for i, name in enumerate(h5_lst):
        final_name = name.split('/')[-1]
        if final_name[:prefix_len] != datname:
            h5_lst.pop(i)
    return h5_lst

def plot_h5_diel(directory, axs = None):
    eps_lst = get_h5_list('eps', directory)
    if len(eps_lst) > 0:
        if axs is None:
            fig, axs = plt.subplots(3)
        for i in range(3):
            vals,_ = get_fields_from_file(eps_lst[0], entry='eps', slice_ax=i)
            axs[i].imshow(vals, cmap='bwr')
        fig.savefig(directory + '/eps.pdf')
        plt.clf()

def get_closest_time(desired, prefix, field_label='ex-'):
    '''Finds the time closest to the desired time within the set of files specified by time_lst
    Returns: the filename of the closest time along with a floating point value representing this time
    '''
    ft_lft_digits = desired.find(".")
    ft_rht_digits = len(desired) - ft_lft_digits - 1#minus 1 because of the decimal point

    #it isn't guaranteed that the exact desired timestep will be available, so we use the closest available time
    fname_regex = prefix + "/" + field_label + "*" + desired + "*.h5"
    time_lst = sorted( glob.glob(fname_regex) )
    #make sure that the file name has the suffix ".h5"
    j = 0
    fname = ""
    while j < len(time_lst) and (len(fname) < 3 or fname[-3:] != '.h5'):
        fname = time_lst[j]
        j += 1
    #read the exact time at which the simulation is performed
    time_name = float(fname[fname.find(field_label)+len(field_label):-3])
    return fname, time_name

def get_point_fields(z):
    '''Iterates over all available hdf5 files and finds the field at the specified point z'''
    ret = np.zeros(len(h5_lst))
    z_ind = int(z/z_step + 0.5)
    #make sure that the specified z value is valid
    if z_ind > n_z_pts:
        raise ValueError("z={} is out of bounds".format(z))

    #if it is, read each h5 file
    for i, fname in enumerate(h5_lst):
        with h5py.File(fname, "r") as f:
            ret[i] = f["ex.r"][z_ind]
    return ret

def get_fields_from_file(fname, slice_ax=-1, slice_ind=-1, n_dims=3, max_z=1.0, entry='ex.r'):
    #load the simulated data from the hdf5 file and plot it
    with h5py.File(fname, "r") as f:
        field_list = np.array(f[entry])
        #check that the user wants to perform a slice
        if n_dims > 0 and (slice_ax >= 0 and slice_ax <= n_dims):
            #select the slice index to be halfway through by default
            if slice_ind < 0:
                slice_ind = field_list.shape[slice_ax]//2
            #return whichever slice the caller requested
            if slice_ax == 0:
                field_list = np.transpose(field_list[slice_ind, :, :])
            elif slice_ax == 1:
                field_list = np.transpose(field_list[:, slice_ind, :])
            elif slice_ax == 2:
                field_list = np.transpose(field_list[:, :, slice_ind])
        z_pts = np.linspace(0, max_z, num=field_list.shape[0])
    return field_list, z_pts

class Source:
    def __init__(self, config):
        #source info
        self.src_loc = float(config["pulse_shape"]["pulse_loc_z"])
        self.devs = float(config["pulse_shape"]["pulse_cutoff"])
        self.pulse_width = float(config["pulse_shape"]["pulse_width"])
        self.post_pulse_runtime = float(config["pulse_shape"]["post_pulse_runtime"])
        self.peak_freq = 2*np.pi*float(config["pulse_shape"]["frequency"])
        src_mon_dist = float(config["monitors"]["near_rad"])
        um_scale = float(config["simulation"]["um_scale"])
        #calculated values
        self.src_loc_l = self.src_loc - src_mon_dist
        self.src_loc_r = self.src_loc + src_mon_dist
        #for the analytic expression we decompose the gaussian source into its fourier components
        self.t_0 = self.devs*self.pulse_width/2
        self.E_0 = -1.0
        self.rt_2_i = np.sqrt(2)*self.t_0*1j
        self.min_freq = self.peak_freq - FREQ_RANGE_SCALE*self.devs*um_scale/self.pulse_width
        self.max_freq = self.peak_freq + FREQ_RANGE_SCALE*self.devs*um_scale/self.pulse_width
        self.pulse_freqs = np.arange(self.min_freq, self.max_freq, PULSE_RES)

        #create an array with the pulse components
        self.pulse_comps = self.E_0*self.pulse_width*np.exp(-0.5*(self.pulse_width*(self.pulse_freqs-self.peak_freq))**2) \
                            * np.exp(1j*BEST_PHASE)/np.sqrt(8*np.pi)

class Geometry:
    def __init__(self, fname):
        #read parameters from the params.conf file
        config = configparser.ConfigParser()
        config.read(fname)
        #geometry info
        self.pml_thick = float(config["simulation"]["pml_thickness"])
        self.n_dims = int(config["simulation"]["dimensions"])
        self.um_scale = float(config["simulation"]["um_scale"])
        self.length = float(config["simulation"]["length"])
        self.middle_w = float(config["junction"]["middle_w"])
        self.junc_max_z = float(config["junction"]["junc_max_z"])
        self.eps_1 = float(config["physical"]["eps_1"])
        try:
            self.eps_2 = float(config["physical"]["eps_2"])
        except KeyError:
            self.eps_2 = self.eps_1
        self.src = Source(config)

        #calculated values
        self.tot_len = self.length + 2*self.pml_thick
        self.z_center = self.tot_len / 2
        self.vol_loc_l = self.pml_thick/self.um_scale
        self.vol_loc_r = (self.tot_len - self.pml_thick)/self.um_scale
        self.l_junc = (self.z_center-self.middle_w/2)/self.um_scale
        self.r_junc = (self.z_center+self.middle_w/2)/self.um_scale
        self.b_junc = (self.z_center+self.junc_max_z/2)/self.um_scale
        self.t_junc = (self.z_center-self.junc_max_z/2)/self.um_scale
        self.lft_x = 0.5 - self.middle_w/(4*self.z_center)
        self.rht_x = 0.5 + self.middle_w/(4*self.z_center)
        self.top_y = 0.5 + self.junc_max_z/(4*self.z_center)
        self.bot_y = 0.5 - self.junc_max_z/(4*self.z_center)

        #calculate indices of refraction and reflection coefficients based on dielectric constants
        self.n_1 = np.sqrt(self.eps_1)
        self.n_2 = np.sqrt(self.eps_2)
        self.coeff_r = abs(self.n_1*self.eps_2 - self.n_2*self.eps_1)/(self.n_1*self.eps_2 + self.n_2*self.eps_1)
        self.coeff_t = 2*self.n_2/(self.n_1*self.eps_2 + self.n_2*self.eps_1)
        print("R={}, T={}, sum={}".format(self.coeff_r, self.coeff_t, self.coeff_r**2 + self.coeff_t**2))

        #if there is a frequency dependence in the dielectric constant, account for it
        try:
            suscep_str = config["physical"]["susceptibilities_2"]
        except KeyError:
            suscep_str = ""
        self.eps_2_func = dielectric_func(suscep_str, self.eps_2, um_scale=self.um_scale)
        self.eps_2s = self.eps_2_func.get_eps_arr(self.src.pulse_freqs*2*np.pi)
        self.n_2s = np.sqrt(self.eps_2s)

        #if the dielectric function has a frequency dependence, then so will the reflection and transmission coefficients and the refractive index
        self.k_2_rs = np.real(self.n_2s)/LIGHT_SPEED
        self.k_2_is = np.imag(self.n_2s)/LIGHT_SPEED
        self.coeff_rs = (self.n_2s*self.n_1 - self.eps_1) / (self.n_2s*self.n_1 + self.eps_1)
        self.coeff_ts = 2*self.n_2s*self.n_1 / (self.eps_2s + self.n_2s*self.eps_1)
        self.rt_cent = abs(self.z_center-self.src.src_loc)*self.n_1/LIGHT_SPEED

    #convert micrometers into the length units used by meep
    def um_to_meep_len(self, l):
        return l * self.um_scale

    #convert seconds into the time units used by meep
    def sec_to_meep_time(self, t):
        #speed of light in um/sec = 299792458000000
        return 299792458000000*t / self.um_scale

    #convert the length units used by meep into micrometers
    def meep_len_to_um(self, l):
        return l / self.um_scale

    #convert the time units used by meep into seconds
    def meep_time_to_sec(self, t):
        #speed of light in um/sec = 299792458000000
        return self.um_scale*t / 299792458000000

    #this is a more general version of get_field_x that uses the fourier decomposition, allowing for the inclusion of a dispersion relation
    def get_electric(r, t, c=LIGHT_SPEED, p=0):
        rt = abs(r-self.src.src_loc)*self.n_1/c
        #since we work in c=1 units k=omega
        #if we're to the left of the barrier, consider the source and reflected
        if r < self.z_center:
            rt_c = abs(r-self.z_center)*self.n_1/c
            return np.array([np.sum( self.src.pulse_comps*np.exp(1j*self.src.pulse_freqs*(t-t_0-rt)) \
                    - self.coeff_rs*self.src.pulse_comps*np.exp(1j*self.src.pulse_freqs*(t-t_0-rt_c-self.rt_cent)) ) \
                        *np.exp(1j*p)*PULSE_RES, 0, 0])
        else:
            decay = abs(r-self.z_center)*self.k_2_is*self.src.pulse_freqs
            rt_c = abs(r-self.z_center)*self.k_2_rs
            #decay = 0
            #rt_c = abs(r-z_center)*n_2/c
            return np.array([np.sum( self.coeff_ts*self.src.pulse_comps*np.exp(1j*self.src.pulse_freqs*(t-t_0-rt_c-self.rt_cent) - decay) )
                                        *np.exp(1j*p)*self.src.pulse_res, 0, 0])
        #return np.sum(pulse_comps*( np.exp(-1j*pulse_freqs*((r-src_loc)+t))))

    def get_magnetic(r, t, c=1.0, p=0):
        return np.roll(self.get_electric(r, t, c=c, p=p), 1)

    def light_posterior(c_pts, simul_pts, left_pml, right_pml, sigma_data=0.1):
        '''Perform a least squares fit of the speed of light to the observed data'''
        def get_sq_err(c):
            anyl_pts = np.array([np.real(self.get_electric(z, time, c=c))[0] for z in z_pts]) 
            return np.exp( np.sum(-0.5*(anyl_pts[left_pml+1:right_pml]-simul_pts[left_pml+1:right_pml])**2/sigma_data) )

        return np.array([get_sq_err(c) for c in c_pts])

    def phase_posterior(p_pts, simul_pts, left_pml, right_pml, sigma_data=0.1):
        '''Perform a least squares fit of the speed of light to the observed data'''
        def get_sq_err(p):
            anyl_pts = np.array([np.real(get_electric(z, time, p=p))[0] for z in z_pts]) 
            return np.exp( np.sum(-0.5*(anyl_pts[left_pml+1:right_pml]-simul_pts[left_pml+1:right_pml])**2/sigma_data) )

        return np.array([get_sq_err(p) for p in p_pts])

    def plot_h5_fields(self, fname, compare_anyl, time=-1, axs=None, er_axs=None):
        '''Generate a plot of the fields in the file specified by fname and place it on the axs specified
        fname: filename of .h5 file to plot
        compare_anyl: if True, then analytic solutions will be calculated and drawn
        axs: the axs on which the plot should be drawn. If this is set to None (default) then a new figure and set of axes is created and saved to fname.pdf
        er_axs: the axes on which the errors between simulated and analytic results should be drawn. If set to none, then no drawing is performed
        '''
        fig = None
        if axs is None:
            fig, axs = plt.subplots()
        if time < 0:
            compare_anyl=False

        field_list, z_pts = get_fields_from_file(fname, slice_ax=0, max_z=self.tot_len)
        #calculate the analytic data
        if compare_anyl:
            field_anyl = np.array([np.real(get_electric(z, time))[0] for z in z_pts])
        #figure out the index of the left and right absorbers
        left_pml = int( PML_EXCLUDE_FACT*(self.pml_thick) * field_list.shape[0] / (2*self.z_center) )
        right_pml = int( PML_EXCLUDE_FACT*(self.tot_len-self.pml_thick) * field_list.shape[0] / (2*self.z_center) )
        n_dims = len(field_list.shape)
        #plot the image in space
        if n_dims == 1:
            axs.plot(z_pts, field_list, linestyle=':', color='black')
            if compare_anyl:
                axs.plot(z_pts, field_anyl, color='black')
            axs.set_ylim(FIELD_RANGE)
            axs.set_xlim((z_pts[0], z_pts[-1]))
            #shade the pml region out since it is unphysical
            axs.fill_between(z_pts, FIELD_RANGE[0], FIELD_RANGE[1], where=z_pts<self.vol_loc_l, color='red', alpha=0.3)
            axs.fill_between(z_pts, FIELD_RANGE[0], FIELD_RANGE[1], where=z_pts>self.vol_loc_r, color='red', alpha=0.3)
            axs.axvline(self.z_center, color='gray')
            axs.set_ylabel(r"$t={}$".format(time))
            if i < len(field_times)-1:
                axs.tick_params('x', labelbottom=False)
            else:
                axs.set_xlabel(r"$z$")
        else:
            axs.imshow(field_list, vmin=FIELD_RANGE[0], vmax=FIELD_RANGE[1], extent=[0, 2*self.z_center/self.um_scale, 0, 2*self.z_center/self.um_scale], cmap='bwr')
            axs.axvline(self.l_junc, self.bot_y, self.top_y, color='gray')
            axs.axvline(self.r_junc, self.bot_y, self.top_y, color='gray')
            axs.axhline(self.t_junc, 0, self.lft_x, color='gray')
            axs.axhline(self.b_junc, 0, self.lft_x, color='gray')
            axs.axhline(self.t_junc, self.rht_x, 2*self.z_center, color='gray')
            axs.axhline(self.b_junc, self.rht_x, 2*self.z_center, color='gray')
            if time >= 0:
                axs.set_title("t={}".format(time))
            if n_dims == 2:
                axs.set_xlabel('y')
                axs.set_ylabel('x')
            else:
                axs.set_xlabel('y')
                axs.set_ylabel('z')

        #we can't calculate the errors if there's no analytic info provided
        if not compare_anyl:
            return fig, 0

        #find the integral square error between the analytic and simulated fields. Note that we have to exclude the absorbing boundary conditions as they are unphysical
        sq_err = (field_anyl-field_list)**2
        #plot the image in error space
        if er_axs is not None:
            if n_dims == 1:
                er_axs.plot(z_pts, sq_err, color='black')
                er_axs.set_ylim(sq_er_range)
                er_axs.set_xlim((z_pts[0], z_pts[-1]))
                er_axs.axvline(z_center, color='gray')
                er_axs.fill_between(z_pts, sq_er_range[0], sq_er_range[1], where=z_pts<self.vol_loc_l, color='red', alpha=0.3)
                er_axs.fill_between(z_pts, sq_er_range[0], sq_er_range[1], where=z_pts>self.vol_loc_r, color='red', alpha=0.3)
            else:
                tmp_im = er_axs.imshow(field_list, cmap='bwr')

        #integrate to find the mean squared error
        return fig, np.sum(sq_err[left_pml+1:right_pml])/(right_pml-left_pml-2)

    def get_cross_section_fields(self, fname, z_heights):
        field_list, z_pts = get_fields_from_file(fname, slice_ax=0, max_z=self.tot_len)
        ret_fields = []
        #iterate through each desired cross-section
        for i, z in enumerate(z_heights):
            ind = round(len(z_pts)*z/(z_pts[-1] - z_pts[0]))
            ret_fields.append(field_list[ind, :])
        return np.array(ret_fields), z_pts

    def plot_cross_section(self, fname, z_heights, axs=None, max_z=1.0):
        '''Helper function to plot an array of z cross sections onto the matplotlib axes specified by axs.
        h5_fields: name of the file to read from
        z_heights: heights for the desired cross sections
        axs: an array of axes to plot each z_height onto. If this is not None, then it must have the same dimensions as z_heights
        '''
        fig = None
        if axs is None:
            fig, axs = plt.subplots(len(z_heights))
        if len(axs) != len(z_heights):
            raise ValueError("z_heights={} and axs={} must have the same length".format(z_heights, axs));

        ret_fields,z_pts = self.get_cross_section_fields(fname, z_heights)
        #iterate through each desired cross-section
        for i, ax in enumerate(axs):
            ax.plot(z_pts, ret_fields[i])
            ax.set_ylim(FIELD_RANGE)
            ax.set_xlim((z_pts[0], z_pts[-1]))
            #shade the pml region out since it is unphysical
            ax.fill_between(z_pts, FIELD_RANGE[0], FIELD_RANGE[1], where=z_pts<self.vol_loc_l, color='red', alpha=0.3)
            ax.fill_between(z_pts, FIELD_RANGE[0], FIELD_RANGE[1], where=z_pts>self.vol_loc_r, color='red', alpha=0.3)
            ax.axvline(self.l_junc, color='gray')
            ax.axvline(self.r_junc, color='gray')

        return fig, axs
