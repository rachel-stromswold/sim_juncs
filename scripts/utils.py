import numpy as np
import configparser
import glob
import h5py

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PML_EXCLUDE_FACT = 2.0
FIELD_RANGE = (-1.0, 1.0)

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

class Geometry:
    def __init__(self, fname, gap_width=-1, gap_thick=-1):
        #read parameters from the params.conf file
        config = configparser.ConfigParser()
        config.read(fname)
        #geometry info
        self.pml_thick = float(config["simulation"]["pml_thickness"])
        self.n_dims = int(config["simulation"]["dimensions"])
        self.um_scale = float(config["simulation"]["um_scale"])
        self.length = float(config["simulation"]["length"])
        #calculate the total side length of the simulation
        self.tot_len = self.length + 2*self.pml_thick
        self.z_center = self.tot_len / 2
        self.gap_width = gap_width
        self.gap_thick = gap_thick
        if gap_width < 0:
            try:
                self.gap_width = float(config["junction"]["gap_width"])
            except:
                pass
        if gap_thick < 0:
            try:
                self.gap_thick = float(config["junction"]["gap_thick"])
            except:
                pass
        #if there is a frequency dependence in the dielectric constant, account for it
        try:
            suscep_str = config["physical"]["susceptibilities_2"]
        except KeyError:
            suscep_str = ""
        #figure out the ambient eps or set it to 1 (vacuum) if there isn't a provided value
        try:
            self.eps_1 = float(config["physical"]["ambient_eps"])
        except KeyError:
            self.eps_1 = 1
        #set the dielectric constant for the second material
        try:
            self.eps_2 = float(config["physical"]["eps_2"])
        except KeyError:
            self.eps_2 = self.eps_1

        #calculated values
        self.vol_loc_l = self.pml_thick/self.um_scale
        self.vol_loc_r = (self.tot_len - self.pml_thick)/self.um_scale
        self.l_junc = self.z_center - self.gap_width*self.um_scale/2
        self.r_junc = self.z_center + self.gap_width*self.um_scale/2
        self.b_junc = self.z_center + self.gap_thick*self.um_scale/2
        self.t_junc = self.z_center - self.gap_thick*self.um_scale/2
        self.lft_x = 0.5 - self.gap_width*self.um_scale/(4*self.z_center)
        self.rht_x = 0.5 + self.gap_width*self.um_scale/(4*self.z_center)
        self.top_y = 0.5 + self.gap_thick*self.um_scale/(4*self.z_center)
        self.bot_y = 0.5 - self.gap_thick*self.um_scale/(4*self.z_center)

    #convert micrometers into the length units used by meep
    def um_to_meep_len(self, l):
        return l * self.um_scale

    #convert seconds into the time units used by meep
    def sec_to_meep_time(self, t):
        #speed of light in um/sec = 299792458000000
        return 299792458000000*t / self.um_scale

    #convert electron-volts into the energy units used by meep
    def ev_to_meep_energy(self, ene):
        #The conversion is based on using hc=1 units so that energy and 1/length have the same units note hc=1.23984193 eV.um
        return self.um_scale*ene / 1.23984193

    #convert meep units for electric field into volts per meter
    def mks_e_field_to_meep_field(self, e_mag):
        #We already know how to express energy and charge in terms of meep units, length units are trivial. The electric field has units energy.charge^-1.length^-1. After some algebra we end up with a conversion factor of um_scale*3e14*hc*6.02e-19/sqrt(4pi)
        return 0.0595261143716197*e_mag/self.um_scale

    #convert coulombs to meep charge units
    def meep_charge_to_coulomb(self, chrg):
        #consider two charges of one meep charge unit separated by 1 um. The force (in meep unis is (q_m/um_scale)^2 where q_m is meep units charge. Note u_0=1 in meep. Note 3.5449... is sqrt(4pi)
        return chrg*self.um_scale*299792458000000/3.5449077018110318

    #convert the length units used by meep into micrometers
    def meep_len_to_um(self, l):
        return l / self.um_scale

    #convert the time units used by meep into seconds
    def meep_time_to_sec(self, t):
        #speed of light in um/sec = 299792458000000
        return self.um_scale*t / 299792458000000

    #convert electron-volts into the energy units used by meep
    def meep_energy_to_ev(self, ene):
        #The conversion is based on using hc=1 units so that energy and 1/length have the same units note hc=1.23984193 eV.um
        return 1.23984193*ene / self.um_scale

    #convert meep charge units to coulombs
    def meep_charge_to_coulomb(self, chrg):
        #consider two charges of one meep charge unit separated by 1 um. The force (in meep unis is (q_m/um_scale)^2 where q_m is meep units charge. Note u_0=1 in meep. Note 3.5449... is sqrt(4pi)
        return chrg*3.5449077018110318/(self.um_scale*299792458000000)

    #convert meep units for electric field into volts per meter
    def meep_field_to_mks_e_field(self, e_mag):
        #We already know how to express energy and charge in terms of meep units, length units are trivial. The electric field has units energy.charge^-1.length^-1. After some algebra we end up with a conversion factor of um_scale*3e14*hc*6.02e-19/sqrt(4pi)
        return 16.799349504942165*e_mag*self.um_scale

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
            um_max = self.meep_len_to_um(self.tot_len)
            axs.imshow(field_list, vmin=FIELD_RANGE[0], vmax=FIELD_RANGE[1], extent=[0, um_max, 0, um_max], cmap='bwr')
            axs.axvline(self.l_junc/self.um_scale, self.bot_y, self.top_y, color='gray')
            axs.axvline(self.r_junc/self.um_scale, self.bot_y, self.top_y, color='gray')
            axs.axhline(self.t_junc/self.um_scale, 0, self.lft_x, color='gray')
            axs.axhline(self.b_junc/self.um_scale, 0, self.lft_x, color='gray')
            axs.axhline(self.t_junc/self.um_scale, self.rht_x, 2*self.z_center, color='gray')
            axs.axhline(self.b_junc/self.um_scale, self.rht_x, 2*self.z_center, color='gray')
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
        field_list, z_pts = get_fields_from_file(fname, slice_ax=0, max_z=self.tot_len/self.um_scale)
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

        fields, z_pts = self.get_cross_section_fields(fname, z_heights)
        #iterate through each desired cross-section
        for i, (field, ax) in enumerate(zip(fields, axs)):
            ax.plot(z_pts, field)
            ax.set_ylim(FIELD_RANGE)
            ax.set_xlim((z_pts[0], z_pts[-1]))
            #shade the pml region out since it is unphysical
            ax.fill_between(z_pts, FIELD_RANGE[0], FIELD_RANGE[1], where=z_pts<self.vol_loc_l, color='red', alpha=0.3)
            ax.fill_between(z_pts, FIELD_RANGE[0], FIELD_RANGE[1], where=z_pts>self.vol_loc_r, color='red', alpha=0.3)
            ax.axvline(self.l_junc, color='gray')
            ax.axvline(self.r_junc, color='gray')

        return fig, axs, fields
