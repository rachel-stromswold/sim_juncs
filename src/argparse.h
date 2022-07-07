#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <cstring>
#include <errno.h>

#define DEFAULT_PML_THICK	0.3 //0.3 um
#define DEFAULT_LEN	8.0 //8 um
#define DEFAULT_PULSE_FREQ	1.43 //700 nm wavelength in eps_0=mu_0=c=1 units
#define DEFAULT_PULSE_START_T	0.0 //0.3 um
#define DEFAULT_PULSE_WIDTH	1.0
#define DEFAULT_RESOLUTION	20.0

#define BUF_SIZE 1024

typedef unsigned int _uint;

//information for parsing arguments
typedef struct {
    int n_dims = 1;
    const char* out_dir = ".";
    double pml_thickness = DEFAULT_PML_THICK;
    double um_scale = 1.0;
    double resolution = DEFAULT_RESOLUTION;
    long grid_num = -1;
    double len = DEFAULT_LEN;
    double ambient_eps = 1.0;
    char* geom_fname = NULL;
    double eps_2 = 2.0;
    double amp = 1.0;

    double source_z = 2.0;
    double freq = DEFAULT_PULSE_FREQ;
    double gauss_width = DEFAULT_PULSE_WIDTH;
    double gauss_start_time = DEFAULT_PULSE_START_T;
    double gauss_cutoff = 5.0;
    double post_source_t = 10.0;
    double src_mon_dist = 0.2;
} Settings;

inline void cleanup_settings(Settings* s) {
    if (s->geom_fname) free(s->geom_fname);
}

/**
  * Remove the whitespace surrounding a word
  */
inline char* trim_whitespace(char* str, size_t* len) {
    char* start = NULL;
    _uint last_non = 0;
    for (_uint i = 0; str[i] != 0; ++i) {
	if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n') {
	    if (start)
		last_non = i;
	    else
		start = str + i;
	}
    }
    str[last_non+1] = 0;
    if (len) *len = last_non+1;
    return start;
}

typedef struct {
    _uint n_susceptibilities = 0;
    double* eps_2_omega;
    double* eps_2_gamma;
    double* eps_2_sigma;
    bool* eps_2_use_denom;
} susceptibility_list;

inline void cleanup_susceptibility_list(susceptibility_list* s) {
    if (s->eps_2_omega) free(s->eps_2_omega);
    if (s->eps_2_gamma) free(s->eps_2_gamma);
    if (s->eps_2_sigma) free(s->eps_2_sigma);
    if (s->eps_2_use_denom) free(s->eps_2_use_denom);
}

/**
 * Add the list of susceptibilities to the Settings file s. The string should have the format (omega_0,gamma_0,sigma_0),(omega_1,gamma_1,sigma_1),...
 * returns: 0 on success or an error code
 * 		-1 null string
 * 		-2 invalid or empty string
 * 		-3 insufficient memory
 */
inline int parse_susceptibilities(susceptibility_list* s, char* const str) {
    //check that the Settings struct is valid and allocate memory
    if (!s) return -1;
    _uint buf_size = BUF_SIZE;
    s->eps_2_omega = (double*)malloc(buf_size*sizeof(double));
    s->eps_2_gamma = (double*)malloc(buf_size*sizeof(double));
    s->eps_2_sigma = (double*)malloc(buf_size*sizeof(double));
    s->eps_2_use_denom = (bool*)malloc(buf_size*sizeof(int));
    //TODO: get gud
    //if (!(s->eps_2_omega && s->eps_2_gamma && s->eps_2_sigma && s->eps_2_use_denom)) 

    double cur_omega = 0.0;
    double cur_gamma = 0.0;
    double cur_sigma = 0.0;
    int cur_use_denom = 0;

    //used for strtok_r
    char* save_str;
    char* tok;
    //find the first entry
    char* cur_entry = strchr(str, '(');
    char* end = strchr(str, ')');

    //only proceed if we have pointers to the start and end of the current entry
    while (cur_entry && end) {
	//resize buffers if necessary
	if (s->n_susceptibilities == buf_size) {
	    buf_size *= 2;
	    s->eps_2_omega = (double*)realloc(s->eps_2_omega, buf_size*sizeof(double));
	    s->eps_2_gamma = (double*)realloc(s->eps_2_gamma, buf_size*sizeof(double));
	    s->eps_2_sigma = (double*)realloc(s->eps_2_sigma, buf_size*sizeof(double));
	    s->eps_2_use_denom = (bool*)realloc(s->eps_2_use_denom, buf_size*sizeof(double));
	}
	cur_omega = 0.0;
	cur_gamma = 0.0;
	cur_use_denom = 1;
	//the open paren must occur before the end paren
	if (cur_entry > end) return -2;

	//null terminate the parenthesis and tokenize by commas
	end[0] = 0;
	//read the omega value
	tok = trim_whitespace( strtok_r(cur_entry+1, ",", &save_str), NULL );
	cur_omega = strtod(tok, NULL);
	if (errno) { errno = 0;return -2; }
	//read the gamma value
	tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	if (!tok) return -2;
	cur_gamma = strtod(tok, NULL);
	if (errno) { errno = 0;return -2; }
	//read the sigma value
	tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	if (!tok) return -2;
	cur_sigma = strtod(tok, NULL);
	if (errno) { errno = 0;return -2; }
	//read the (optional) use denom flag
	tok = strtok_r(NULL, ",", &save_str);
	if (tok) {
	    tok = trim_whitespace(tok, NULL);
	    cur_use_denom = strtol(tok, NULL, 10);
	    if (errno) {
		errno = 0;
		cur_use_denom=1;
		if (strcmp(tok, "drude") == 0 || strcmp(tok, "true") == 0) cur_use_denom = 0;
	    }
	    if (strcmp(tok, "lorentz") == 0 || strcmp(tok, "false") == 0) cur_use_denom = 1;
	}

	//save the information
	s->eps_2_omega[s->n_susceptibilities] = cur_omega;
	s->eps_2_gamma[s->n_susceptibilities] = cur_gamma;
	s->eps_2_sigma[s->n_susceptibilities] = cur_sigma;
	s->eps_2_use_denom[s->n_susceptibilities] = cur_use_denom;
	++s->n_susceptibilities;

	//advance to the next entry
	if (end[1] == 0) break;
	cur_entry = strchr(end+1, '(');
	end = strchr(cur_entry, ')');
    }

    //realloc arrays to exactly fit memory
    buf_size = s->n_susceptibilities;
    s->eps_2_omega = (double*)realloc(s->eps_2_omega, buf_size*sizeof(double));
    s->eps_2_gamma = (double*)realloc(s->eps_2_gamma, buf_size*sizeof(double));
    s->eps_2_sigma = (double*)realloc(s->eps_2_sigma, buf_size*sizeof(double));
    s->eps_2_use_denom = (bool*)realloc(s->eps_2_use_denom, buf_size*sizeof(int));

    return 0;
}

/**
 * Parse a single token value pair and save it to the settings struct
 */
inline void handle_pair(Settings* s, char* const tok, _uint toklen, char* const val, _uint vallen) {
    //remove trailing whitespace from value and token strings
    while (toklen > 0 && (tok[toklen-1] == ' ' || tok[toklen-1] == '\t')) --toklen;
    while (vallen > 0 && (val[vallen-1] == ' ' || val[vallen-1] == '\t')) --vallen;
    tok[toklen] = 0;val[vallen] = 0;

    //reset errors
    errno = 0;
    if (strcmp(tok, "pml_thickness") == 0) {
	s->pml_thickness = strtod(val, NULL);
    } else if (strcmp(tok, "dimensions") == 0) {
	s->n_dims = strtol(val, NULL, 10);
    } else if (strcmp(tok, "um_scale") == 0) {
	s->um_scale = strtod(val, NULL);
    } else if (strcmp(tok, "length") == 0) {
	s->len = strtod(val, NULL);
    } else if (strcmp(tok, "frequency") == 0) {
	s->freq = strtod(val, NULL);
    } else if (strcmp(tok, "post_pulse_runtime") == 0) {
	s->post_source_t = strtod(val, NULL);
    } else if (strcmp(tok, "pulse_ramp_start") == 0) {
	s->gauss_start_time = strtod(val, NULL);
    } else if (strcmp(tok, "pulse_amplitude") == 0) {
	s->amp = strtod(val, NULL);
    } else if (strcmp(tok, "pulse_width") == 0) {
	s->gauss_width = strtod(val, NULL);
    } else if (strcmp(tok, "pulse_cutoff") == 0) {
	s->gauss_cutoff = strtod(val, NULL);
    } else if (strcmp(tok, "pulse_loc_z") == 0) {
	s->source_z = strtod(val, NULL);
    } else if (strcmp(tok, "ambient_eps") == 0) {
	s->ambient_eps = strtod(val, NULL);
    } else if (strcmp(tok, "geom_fname") == 0) {
	s->geom_fname = trim_whitespace(strdup(val), NULL);
    } else if (strcmp(tok, "eps_2") == 0) {
	s->eps_2 = strtod(val, NULL);
    } else if (strcmp(tok, "near_rad") == 0) {
	s->src_mon_dist = strtod(val, NULL);
    }
}

/**
 * Parse the .ini style configuration file specified by fname
 */
inline int parse_conf_file(Settings* s, char* fname) {
    char buf[BUF_SIZE];
    _uint toklen = 0;
    char* tok = NULL;
    _uint val_s = 0;
    char* val = NULL;

    FILE* fp = fopen(fname, "r");
    //iterate until we reach the end of the file
    while (fgets(buf, BUF_SIZE, fp)) {
	_uint i = 0;
	for (; i < BUF_SIZE && (buf[i] == ' ' || buf[i] == '\t'); ++i);//ignore leading whitespace
	//ignore everything contained in a [] comment block. These are only allowed to occur on the first non whitespace character
	if (buf[i] == '[') {
	    for (; i < BUF_SIZE && buf[i] != ']'; ++i);
	    ++i;
	}
	//whatever precedes an = sign is a token
	tok = buf + i;
	toklen = i;
	//read until we reach the end of the line
	for (; i < BUF_SIZE-1 && buf[i] && buf[i] != '#'; ++i) {
	    //only read into the token string if we don't have an active value string
	    if (!val) {
		if (buf[i] == '=') {
		    buf[i] = 0;//null terminate the token
		    val = buf+i+1;
		    val_s = i;
		    //save the length of the token if we start reading a value
		    toklen = i-toklen;
		}
	    } else if (buf[i] == '\n' || buf[i] == ';') {
		buf[i] = 0;
		handle_pair(s, tok, toklen, val, i-val_s);
		//start a new token
		tok = buf+i+1;
		toklen = i+1;
		val = NULL;
	    }
	}
    }
    fclose(fp);

    //correct everything which involves length units to have the proper scale
    if (s->um_scale != 1.0) {
        s->freq /= s->um_scale;
	s->len *= s->um_scale;
	s->pml_thickness *= s->um_scale;
    }
    return 0;
}

//parse arguments used by the test program and remove them from the list so that they may safely be passed to open_mp initialization, returns 0 if no errors were detected and an error code otherwise
inline int parse_args(Settings* a, int* argc, char ** argv) {
    double tmp_resolution = -1.0;

    int n_args = *argc;//dereferencing that pointer all the time is kinda annoying
    if (a) {
	for (_uint i = 0; i < n_args; ++i) {
	    _uint to_rm = 0;
	    //look for known argument names
	    if (strstr(argv[i], "--out-dir") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --out-dir <prefix to output hdf5 files to>");
		    return 0;
		} else {
		    a->out_dir = argv[i+1];
		    to_rm = 2; //we need to remove both this and the next argument from the list
		}
	    } else if (strstr(argv[i], "--grid-res") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --grid-res <grid points per unit length>");
		    return 0;
		} else {
		    //a->resolution = strtod(argv[i+1], NULL);
		    double res = strtod(argv[i+1], NULL);
            //ensure that we have an odd number of grid points
            a->grid_num = 2*int(0.5*(1 + (a->len)*res)) + 1;
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --grid-res");
			return errno;
		    }
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--grid-num") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --grid-res <grid points per unit length>");
		    return 0;
		} else {
		    a->grid_num = strtol(argv[i+1], NULL, 10);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --grid-res");
			return errno;
		    }
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--length") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --length <length of each axis of the simulation volume in arbitrary units>");
		    return 0;
		} else {
		    a->len = strtod(argv[i+1], NULL);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --length");
			return errno;
		    }
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--freq") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --freq <frequency of incident light source>");
		    return 0;
		} else {
		    a->freq = strtod(argv[i+1], NULL);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --freq");
			return errno;
		    }
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--eps1") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --eps1 <epsilon1>");
		    return 0;
		} else {
		    a->ambient_eps = strtod(argv[i+1], NULL);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --eps1");
			return errno;
		    }
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--eps2") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --eps2 <epsilon1>");
		    return 0;
		} else {
		    a->eps_2 = strtod(argv[i+1], NULL);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --eps2");
			return errno;
		    }
		    to_rm = 2;
		}
	    }

	    //if we need to remove any arguments, do that now
	    if (to_rm) {
		_uint j = i;
		for (; j+to_rm < n_args; ++j) argv[j] = argv[j+to_rm];
		for (; j < n_args; ++j) argv[j] = NULL;
		n_args -= to_rm;
		--i; //don't skip the next entry after removal
	    }
	}

	//if a number of grid points was specified, use that
	if (a->grid_num > 0) {
	    double total_len = a->len + 2*a->pml_thickness;
	    a->resolution = (double)(a->grid_num) / total_len;
	}
    }
    *argc = n_args;

    return 0;
}

/**Places a string with in str the format 000aa_bb0 where aa.bb is the current time t and the result is zero padded such that it has a constant length regardless of input. str is guaranteed to be null terminated even in the event of an error.
 * returns the number of characters written (excluding the null terminator) on success or an error code listed below:
 *	-1: there was not sufficient space allocated to store the string
 *	-2: the number of digits exceeded the number specified
 **/
inline size_t make_dec_str(char* str, size_t len, double t, _uint n_digits_a, _uint n_digits_b) {
    //checks to avoid segfaults
    if (!str) return -1;
    if (len < n_digits_a+n_digits_b+2) {
	if (len > 0) str[0] = 0;
	return -1;
    }
    long bef = (long)t;
    t -= (double)bef;//we only want the parts of the decimal less than 1

    //write the portion to the left of the decimal
    size_t i = 0;
    do {
	if (i == n_digits_a) { 
	    str[n_digits_a] = 0;
	    return -2;
	}
	str[n_digits_a-i-1] = '0' + (char)(bef % 10);
	bef /= 10;
	++i;
    } while (bef && i < len-1);

    //write the preceeding zeros
    for (size_t j = 0; j < n_digits_a-i; ++j) str[j] = '0';

    //write the "decimal" place
    str[n_digits_a] = '.';
    i = n_digits_a+1;
    while (i < n_digits_b+n_digits_a+1) {
	t *= 10;
	str[i] = '0' + (char)t;
	t -= floor(t);
	++i;
    }

    //null terminate
    str[i] = 0;
    return i;
}

#endif //ARGPARSE_H
