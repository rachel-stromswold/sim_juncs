#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <cstring>
#include <errno.h>
#include <math.h>

#define DEFAULT_PML_THICK	0.3 //0.3 um
#define DEFAULT_LEN	8.0 //8 um
#define DEFAULT_PULSE_FREQ	1.43 //700 nm wavelength in eps_0=mu_0=c=1 units
#define DEFAULT_PULSE_START_T	0.0 //0.3 um
#define DEFAULT_PULSE_WIDTH	1.0
#define DEFAULT_RESOLUTION	2.0

#define DEFAULT_SMOOTH_N	1
#define DEFAULT_SMOOTH_RAD	0.05

#define BUF_SIZE 1024
#define SMALL_BUF_SIZE 256

typedef unsigned int _uint;

//information for parsing arguments
typedef struct {
    //system stuff
    const char* out_dir = "/tmp";
    char* geom_fname_al = NULL;
    char* geom_fname = NULL;
    char* conf_fname = NULL;
    //problem geometry
    int n_dims = -3;
    double pml_thickness = -DEFAULT_PML_THICK;
    double len = -DEFAULT_LEN;
    double um_scale = 1.0;
    long grid_num = -1;
    double resolution = -DEFAULT_RESOLUTION;
    double courant = 0.5;//0.5 is the meep default
    //simulation parameters
    double ambient_eps = -1.0;
    _uint smooth_n = 0;
    double smooth_rad = DEFAULT_SMOOTH_RAD;
    //misc
    double post_source_t = 10.0;
    _uint field_dump_span = 20;
    unsigned char dump_raw = 0;
    _uint verbosity = 1;//same as meep
    char* user_opts = NULL;
} parse_settings;

/**
 * Saves a string with the desired numbers with zeros appended to the front such that the total length of the string is n_digits
 * str: c string to write to
 * buf_size: at most buf_size-1 bytes will be written to str, the string is gauranteed to be NULL terminated
 * num: the number to write
 * n_digits: the number of digits to format the string to
 * returns 0: on success, -1 on failure
 */
inline int write_number(char* str, size_t buf_size, int num, size_t n_digits) {
    if (buf_size < n_digits + 1) return -1;
    //null terminate
    str[n_digits] = 0;
    if (n_digits == 0) return 0;
    //handle negative numbers
    if (num < 0) {
	str[0] = '-';
	num *= -1;
	str += 1;
	buf_size -= 1;
	n_digits -= 1;
    }
    int k = n_digits-1;
    do {
	//this means that there weren't enough digits to write the whole string
	if (k < 0) return -1;
	str[k--] = (num % 10) + '0';
	num /= 10;
    } while (num > 0);
    for (; k >= 0; --k) str[k] = '0';
    return 0;
}

/*
 * handy utility function which checks if the pointer is valid and saves to it if it is.
 */
inline void set_ercode(int* sto, int er) {
    if (sto) *sto = er;
}

inline void cleanup_settings(parse_settings* s) {
    if (s->geom_fname_al) free(s->geom_fname_al);
    if (s->user_opts) free(s->user_opts);
}

/**
  * Remove the whitespace surrounding a word
  */
inline char* trim_whitespace(char* str, size_t* len) {
    if (!str) return NULL;
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

/**
 * Parse a single token value pair and save it to the settings struct
 */
inline void handle_pair(parse_settings* s, char* const tok, _uint toklen, char* const val, _uint vallen) {
    //remove trailing whitespace from value and token strings
    while (toklen > 0 && (tok[toklen-1] == ' ' || tok[toklen-1] == '\t')) --toklen;
    while (vallen > 0 && (val[vallen-1] == ' ' || val[vallen-1] == '\t')) --vallen;
    tok[toklen] = 0;val[vallen] = 0;

    //reset errors
    errno = 0;
    if (strcmp(tok, "pml_thickness") == 0 && s->pml_thickness < 0) {
        s->pml_thickness = strtod(val, NULL);
    } else if (strcmp(tok, "resolution") == 0 && s->resolution < 0) {
        s->resolution = strtod(val, NULL);
    } else if (strcmp(tok, "dimensions") == 0 && s->n_dims < 0) {
        s->n_dims = strtol(val, NULL, 10);
    } else if (strcmp(tok, "um_scale") == 0) {
        s->um_scale = strtod(val, NULL);
    } else if (strcmp(tok, "length") == 0 && s->len < 0) {
        s->len = strtod(val, NULL);
    } else if (strcmp(tok, "smooth_rad") == 0) {
        s->smooth_rad = strtod(val, NULL);
    } else if (strcmp(tok, "smooth_n") == 0) {
        s->smooth_n = strtol(val, NULL, 10);
    } else if (strcmp(tok, "courant") == 0) {
        s->courant = strtod(val, NULL);
    } else if (strcmp(tok, "ambient_eps") == 0 && s->ambient_eps < 0) {
        s->ambient_eps = strtod(val, NULL);
    } else if (strcmp(tok, "geom_fname") == 0 && !s->geom_fname_al) {
        //copy only the non whitespace portion
	if (s->geom_fname_al) free(s->geom_fname_al);
        s->geom_fname_al = strdup(val);
        s->geom_fname = trim_whitespace(s->geom_fname_al, NULL);    
    } else if (strcmp(tok, "post_source_t") == 0) {
        s->post_source_t = strtod(val, NULL);
    } else if (strcmp(tok, "field_dump_span") == 0) {
        s->field_dump_span = strtol(val, NULL, 10);
    } else if (strcmp(tok, "dump_raw") == 0) {
        if (strcmp(val, "true") == 0)
            s->dump_raw = 1;
        else if (strcmp(val, "false") == 0)
            s->dump_raw = 0;
        else
            s->dump_raw = strtol(val, NULL, 10);
    }
}

//helper function to automatically adjust parameters if necessary
inline void correct_defaults(parse_settings* s) {
    double total_len = s->len + 2*s->pml_thickness;

    //if a number of grid points was specified, use that
    if (s->grid_num < 0) {
	//ensure that we have an odd number of grid points
	s->grid_num = 2*int( 0.5*(1 + (total_len)*(s->resolution)) ) + 1;
    }
	
    //simple but hacky way to let users override the params.conf file if needed
    s->n_dims = abs(s->n_dims);
    s->pml_thickness = fabs(s->pml_thickness);
    s->len = fabs(s->len);
    s->resolution = fabs(s->resolution);
    //simulation parameters
    s->ambient_eps = fabs(s->ambient_eps);

    //adjust the resolution if we are to use a grid number for specification
    s->resolution = (double)(s->grid_num) / total_len;
}

/**
 * Parse the .ini style configuration file specified by fname
 */
inline int parse_conf_file(parse_settings* s, char* fname) {
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

    correct_defaults(s);
    return 0;
}

//parse arguments used by the test program and remove them from the list so that they may safely be passed to open_mp initialization, returns 0 if no errors were detected and an error code otherwise
inline int parse_args(parse_settings* a, int* argc, char ** argv) {
    _uint n_args = (_uint)(*argc);//dereferencing that pointer all the time is kinda annoying
    if (a) {
	for (_uint i = 0; i < n_args; ++i) {
	    _uint to_rm = 0;
	    //look for known argument names
	    if (strstr(argv[i], "--conf-file") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --conf-file file name");
		    return 0;
		} else {
		    a->conf_fname = argv[i+1];
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--geom-file") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --out-dir <prefix to output hdf5 files to>");
		    return 0;
		} else {
		    //deallocate memory if necessary
		    if (a->geom_fname_al) free(a->geom_fname_al);
		    //copy only the non whitespace portion
		    a->geom_fname_al = strdup(argv[i+1]);
		    a->geom_fname = trim_whitespace(a->geom_fname_al, NULL);
		    to_rm = 2; //we need to remove both this and the next argument from the list
		}
	    } else if (strstr(argv[i], "--out-dir") == argv[i]) {
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
		    a->resolution = strtod(argv[i+1], NULL);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --grid-res");
			return errno;
		    }
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--grid-num") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --grid-num total number of grid points to use");
		    return 0;
		} else {
		    a->grid_num = strtol(argv[i+1], NULL, 10);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid floating point supplied to --grid-num");
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
	    } else if (strstr(argv[i], "-v") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep -v <verbosity (integer)>");
		    return 0;
		} else {
		    a->verbosity = strtol(argv[i+1], NULL, 10);
		    //check for errors
		    if (errno != 0) {
			printf("Invalid value supplied to -v, should be integer");
			return errno;
		    }
		    to_rm = 2;
		}
	    } else if (strstr(argv[i], "--opts") == argv[i]) {
		if (i == n_args-1) {
		    printf("Usage: meep --opts \"<name1=value1;name2=value2>\"");
		    return 0;
		} else {
		    if (a->user_opts) free(a->user_opts);
		    a->user_opts = strdup(argv[i+1]);
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
    }
    *argc = n_args;

    return 0;
}

/**Places a string with in str the format 000aa_bb0 where aa.bb is the current time t and the result is zero padded such that it has a constant length regardless of input. str is guaranteed to be null terminated even in the event of an error.
 * returns the number of characters written (excluding the null terminator) on success or an error code listed below:
 *	-1: there was not sufficient space allocated to store the string
 *	-2: the number of digits exceeded the number specified
 **/
inline long make_dec_str(char* str, size_t len, double t, _uint n_digits_a, _uint n_digits_b, char dec_char='.') {
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
    str[n_digits_a] = dec_char;
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
