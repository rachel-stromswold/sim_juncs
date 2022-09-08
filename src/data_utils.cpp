#include "data_utils.hpp"

// ================================ complex struct ================================

complex conj(complex z) {
    complex ret = {z.re, -z.im};
    return ret;
}

complex c_exp(complex z) {
    if (z.re == 0) {
	complex ret = {(_ftype)cos(z.im), (_ftype)sin(z.im)};
	return ret;
    } else {
	_ftype mult = exp(z.re);
	complex ret = {mult*(_ftype)cos(z.im), mult*(_ftype)sin(z.im)};
	return ret;
    }
}

complex c_mult(complex z, complex w) {
    complex ret = { z.re*w.re - z.im*w.im, z.re*w.im + z.im*w.re };
    return ret;
}

complex c_multd(complex z, _ftype w) {
    complex ret = { z.re*w, z.im*w };
    return ret;
}

_ftype abs_sq(complex z) {
    return z.re*z.re + z.im*z.im;
}

_ftype abs(complex z) {
    return sqrt(abs_sq(z));
}

// ================================ data struct ================================

/**
 * Initialize an empty data array
 */
void make_data_arr(data_arr* dat, _ulong size) {
    if (dat) {
	dat->size = size;
	dat->buf_size = size;
	dat->buf = (complex*)malloc(sizeof(complex)*dat->buf_size);
	//check to make sure allocation was successful
	if (!dat->buf) {
	    dat->size = 0;
	    dat->buf_size = 0;
	}
    }
}

void resize(data_arr* dat, _ulong new_size) {
    if (dat && new_size > dat->buf_size) {
	complex* tmp = (complex*)realloc(dat->buf, sizeof(complex)*new_size);
	if (dat->buf) {
	    dat->buf_size = new_size;
	    dat->buf = tmp;
	}
    }
}

#ifdef __cplusplus
data_arr::data_arr() {
    size = 0;
    buf_size = 0;
    buf = NULL;
}

/**
 * Swap data in o with this
 */
void data_arr::swap(data_arr& o) {
    //swap buf
    complex* tmp_buf = buf;
    buf = o.buf;
    o.buf = tmp_buf;
    //swap size
    _ulong tmp_size = size;
    size = o.size;
    o.size = tmp_size;
    //swap buf_size
    tmp_size = buf_size;
    buf_size = o.buf_size;
    o.buf_size = tmp_size;
}

/**
 * Access the value at location i
 */
complex& data_arr::operator[](size_t i) {
    return buf[i];
}

/**
 * Initialize a data array from a c++ vector
 */
data_arr::data_arr(std::vector<complex> in_pts) {
    size = in_pts.size();
    buf_size = size;
    buf = (complex*)malloc(buf_size);
    //check to make sure allocation was successful
    if (!buf) {
	size = 0;
	buf_size = 0;
    } else {
	//copy data
	for (_ulong i = 0; i < size; ++i) {
	    buf[i] = in_pts[i];
	}
    }
}

/**
 * Initialize an empty data array
 */
data_arr::data_arr(_ulong size) {
    make_data_arr(this, size);
}

/**
 * copy
 */
data_arr::data_arr(const data_arr& o) {
    buf = (complex*)malloc(sizeof(complex)*o.buf_size);
    if (buf) {
	size = o.size;
	buf_size = o.buf_size;
	for (_ulong k = 0; k < buf_size; ++k) buf[k] = o.buf[k];
    } else {
	size = 0;
	buf_size = 0;
    }
}

/**
 * assignment
 */
data_arr& data_arr::operator=(data_arr o) {
    swap(o);
    return *this;
}

/**
 * move
 */
data_arr::data_arr(data_arr&& o) {
    size = o.size;
    buf_size = o.buf_size;
    buf = o.buf;
    o.size = 0;
    o.buf_size = 0;
    o.buf = NULL;
}

/**
 * destructor
 */
data_arr::~data_arr() {
    cleanup_data_arr(this);
}
#endif

/**
 * Deallocate memory for a data array
 */
void cleanup_data_arr(data_arr* dat) {
    if (dat) {
	dat->size = 0;
	dat->buf_size = 0;
	if (dat->buf) free(dat->buf);
	dat->buf = NULL;
    }
}

/**
 * append a point to the end of the data array
 */
void add_point(data_arr* dat, complex val) {
    if (dat->size >= dat->buf_size) {
	dat->buf_size *= 2;
	dat->buf = (complex*)realloc(dat->buf, sizeof(complex)*(dat->buf_size));
    }
    dat->buf[(dat->size)++] = val;
}

void d_add_point(data_arr* dat, _ftype val) {
    complex pt = {0, val};
    add_point(dat, pt);
}

/**
 * These functions perform piecewise arithmetic operations. All operations are performed in place with a modified
 * returns 0 on success, -1 on failure if anything is returned
 */
#ifdef __cplusplus
void pw_conj(data_arr& a) {
#else
void pw_conj(data_arr a) {
#endif
    for (_ulong k = 0; k < a.size; ++k) a.buf[k].im *= -1;
}
#ifdef __cplusplus
void pw_exp(data_arr& a) {
#else
void pw_exp(data_arr a) {
#endif
    for (_ulong k = 0; k < a.size; ++k) a.buf[k] = c_exp(a.buf[k]);
}
#ifdef __cplusplus
int pw_mult(data_arr& a, const data_arr b) {
#else
int pw_mult(data_arr a, const data_arr b) {
#endif
    if (b.size != a.size) return -1;
    _ftype tmp_re = 0;
    _ftype tmp_im = 0;
    for (_ulong k = 0; k < a.size; ++k) {
	tmp_re = a.buf[k].re*b.buf[k].re - a.buf[k].im*b.buf[k].im;
	tmp_im = a.buf[k].re*b.buf[k].im + a.buf[k].im*b.buf[k].re;
	a.buf[k].re = tmp_re;
	a.buf[k].re = tmp_im;
    }
    return 0;
}
#ifdef __cplusplus
void pw_mult_scale(data_arr& a, complex z) {
#else
void pw_mult_scale(data_arr a, complex z) {
#endif
    for (_ulong k = 0; k < a.size; ++k) a.buf[k] = c_mult(a.buf[k], z);
}
#ifdef __cplusplus
int pw_add(data_arr& a, const data_arr b) {
#else
int pw_add(data_arr a, const data_arr b) {
#endif
    if (b.size != a.size) return -1;
    for (_ulong k = 0; k < a.size; ++k) {
	a.buf[k].re = a.buf[k].re + b.buf[k].re;
	a.buf[k].im = a.buf[k].im + b.buf[k].im;
    }
    return 0;
}
#ifdef __cplusplus
void pw_add_scale(data_arr& a, complex z) {
#else
void pw_add_scale(data_arr a, complex z) {
#endif
    for (_ulong k = 0; k < a.size; ++k) {
	a.buf[k].re = a.buf[k].re + z.re;
	a.buf[k].im = a.buf[k].im + z.im;
    }
}
#ifdef __cplusplus
int pw_sub(data_arr& a, const data_arr b) {
#else
int pw_sub(data_arr a, const data_arr b) {
#endif
    if (b.size != a.size) return -1;
    for (_ulong k = 0; k < a.size; ++k) {
	a.buf[k].re = a.buf[k].re - b.buf[k].re;
	a.buf[k].im = a.buf[k].im - b.buf[k].im;
    }
    return 0;
}
#ifdef __cplusplus
void pw_sub_scale(data_arr& a, complex z) {
#else
void pw_sub_scale(data_arr a, complex z) {
#endif
    for (_ulong k = 0; k < a.size; ++k) {
	a.buf[k].re = a.buf[k].re - z.re;
	a.buf[k].im = a.buf[k].im - z.im;
    }
}
#ifdef __cplusplus
void pw_abs_sq(data_arr& a) {
#else
void pw_abs_sq(data_arr a) {
#endif
    for (_ulong k = 0; k < a.size; ++k) {
	a.buf[k].re = a.buf[k].re*a.buf[k].re + a.buf[k].im*a.buf[k].im;
	a.buf[k].im = 0;
    }
}
#ifdef __cplusplus
void pw_abs(data_arr& a) {
#else
void pw_abs(data_arr a) {
#endif
    for (_ulong k = 0; k < a.size; ++k) {
	a.buf[k].re = sqrt(a.buf[k].re*a.buf[k].re + a.buf[k].im*a.buf[k].im);
	a.buf[k].im = 0;
    }
}

/**
 * this is a helper function for fft which computes the left and right terms for the fourier transform of the data specified by dat that only takes elements of the form span*n + rem for integer n between 0 and (dat_arr.size/span - 1)
 * NOTE: this function does not perform NULL safety checks. It is intended to be called by fft() and not end users.
 * dat_arr: data array
 * span: the span to take when indexing the array
 * rem: the remainder to take when indexing the array. For each integer n in the interval [0, dat_arr.size/span), take the nth index to be span*n+rem
 * left: where to save the left portion of the fft
 * right: where to save the right portion
 */
#ifdef __cplusplus
void part_rfft(const dat_helper& help, _ulong span, _ulong rem) {
#else
void part_rfft(dat_helper help, _ulong span, _ulong rem) {
#endif
    //base case, there are only two elements in the partial array. We multiply by 3 since that's the smallest number greater than 2
    if (span*2 >= help.truncated_size) {
	//for each k we have a simple 
	for (_ulong k = 0; k < help.sto.size; ++k) {
	    //val = dat[rem] + e^{-i*(2*pi/N)*span*k}*dat[span+rem]
	    complex phase_span = {0, -help.tau_by_n*span*k};
	    complex val = c_mult(c_exp(phase_span), help.dat.buf[span+rem]);
	    val.re += help.dat.buf[rem].re;
	    val.im += help.dat.buf[rem].im;
	    //sto[k] += e^{-i*(2*pi/N)*rem*k}*val
	    complex phase_rem = {0, -help.tau_by_n*rem*k};
	    val = c_mult(c_exp(phase_rem), val);
	    help.sto.buf[k].re += val.re;
	    help.sto.buf[k].im += val.im;
	    //sto.buf[k] = help.dat.buf[rem] + c_exp(-I*help.tau_by_n*k)*help.dat.buf[span+rem];
	}
    } else {
	//break the transform into two half transforms, each with twice the span and twice the phase factor
	_ulong new_span = 2*span;

	//recursively take left and right branches in the remainder
	part_rfft(help, new_span, rem);
	part_rfft(help, new_span, rem+span);
	/*for (_ulong k = 0; k < tmp_sto.size; ++k) {
	    help.sto.buf[k] += ( sto_left.buf[k] + sto_right.buf[k]*c_exp(-I*phase_fact*span) )*c_exp(-I*phase_fact*span*k);
	}*/
    }
}

/**
 * this is a helper function for fft which computes the left and right terms for the fourier transform of the data specified by dat that only takes elements of the form span*n + rem for integer n between 0 and (dat_arr.size/span - 1).
 * NOTE: this function does not perform NULL safety checks. It is intended to be called by fft() and not end users.
 * dat_arr: data array
 * span: the span to take when indexing the array
 * rem: the remainder to take when indexing the array. For each integer n in the interval [0, dat_arr.size/span), take the nth index to be span*n+rem
 * left: where to save the left portion of the fft
 * right: where to save the right portion
 */
#ifdef __cplusplus
void part_fft(const dat_helper& help, _ulong span, _ulong rem) {
#else
void part_fft(dat_helper help, _ulong span, _ulong rem) {
#endif
    //base case, there are only two elements in the partial array. We multiply by 3 since that's the smallest number greater than 2
    if (span*2 >= help.truncated_size) {
	//for each k we have a simple 
	long half_size = help.sto.size / 2;
	_ulong sto_k = half_size;
	for (long k = -half_size; k < half_size; ++k) {
	    //val = dat[rem] + e^{-i*(2*pi/N)*span*k}*dat[span+rem]
	    complex phase_span = {0, -help.tau_by_n*span*k};
	    complex val = c_mult(c_exp(phase_span), help.dat.buf[span+rem]);
	    val.re += help.dat.buf[rem].re;
	    val.im += help.dat.buf[rem].im;
	    //sto[k] += e^{-i*(2*pi/N)*rem*k}*val
	    complex phase_rem = {0, -help.tau_by_n*rem*k};
	    val = c_mult(c_exp(phase_rem), val);
	    help.sto.buf[sto_k].re += val.re;
	    help.sto.buf[sto_k].im += val.im;
	    //increment sto_k and check if it needs to be wrapped around
	    ++sto_k;
	    if (sto_k == help.sto.size) {
		sto_k = 0;
	    }
	    //sto.buf[k] = help.dat.buf[rem] + c_exp(-I*help.tau_by_n*k)*help.dat.buf[span+rem];
	}
    } else {
	//break the transform into two half transforms, each with twice the span and twice the phase factor
	_ulong new_span = 2*span;

	//recursively take left and right branches in the remainder
	part_rfft(help, new_span, rem);
	part_rfft(help, new_span, rem+span);
	/*for (_ulong k = 0; k < tmp_sto.size; ++k) {
	    help.sto.buf[k] += ( sto_left.buf[k] + sto_right.buf[k]*c_exp(-I*phase_fact*span) )*c_exp(-I*phase_fact*span*k);
	}*/
    }
}

/**
 * Compute the discrete fourier transform for the data series specified by dat. If it is known in advance that all terms in dat are real, calling rfft will be faster. Negative frequency terms are "wrapped around" at the Nyquist frequency. The behavior matches scipy's implementation, for details consult https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fft.html
 */
data_arr fft(const data_arr dat) {
    //figure out how many recursive steps are needed.
    _ulong log_2_n = (_ulong)(log(dat.size) / log(2));
    _ulong truncated_size = 1 << log_2_n;

    //allocate space for the final result
    data_arr final_result;
    make_data_arr(&final_result, truncated_size);
    complex null_val = {0.0, 0.0};
    for (_ulong k = 0; k < final_result.size; ++k) final_result.buf[k] = null_val;

    //make a helper object that keeps track of information over recursive steps
    dat_helper help = {dat, truncated_size, 2*(_ftype)M_PI/dat.size, final_result};

    part_fft(help, 1, 0);
    return help.sto;
}

/**
 * Compute the inverse discrete fourier transform for the data series specified by dat. If it is known in advance that all terms in dat are real, calling rfft will be faster. Negative frequency terms are "wrapped around" at the Nyquist frequency. The behavior matches scipy's implementation, for details consult https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fft.html
 */
data_arr ifft(const data_arr dat) {
    //figure out how many recursive steps are needed.
    _ulong log_2_n = (_ulong)(log(dat.size) / log(2));
    _ulong truncated_size = 1 << log_2_n;

    //allocate space for the final result
    data_arr final_result;
    make_data_arr(&final_result, truncated_size);
    complex null_val = {0.0, 0.0};
    for (_ulong k = 0; k < final_result.size; ++k) final_result.buf[k] = null_val;

    //make a helper object that keeps track of information over recursive steps
    dat_helper help = {dat, truncated_size, -2*(_ftype)M_PI/dat.size, final_result};

    part_fft(help, 1, 0);
    return help.sto;
}

/**
 * Compute the discrete fourier transform for the data series specified by dat. This version assumes that all terms in the data series are real. Thus, negative frequency terms correspond to the conjugate of their corresponding positive frequency terms. To reduce redundancy we only store up to the Nyquist frequency
 */
data_arr rfft(const data_arr dat) {
    //figure out how many recursive steps are needed.
    _ulong log_2_n = (_ulong)(log(dat.size) / log(2));
    _ulong truncated_size = 1 << log_2_n;

    //allocate space for the final result
    data_arr final_result;
    make_data_arr(&final_result, truncated_size/2);
    complex null_val = {0.0, 0.0};
    for (_ulong k = 0; k < final_result.size; ++k) final_result.buf[k] = null_val;

    //make a helper object that keeps track of information over recursive steps
    dat_helper help = {dat, truncated_size, 2*(_ftype)M_PI/dat.size, final_result};

    part_rfft(help, 1, 0);
    return help.sto;
}

/**
 * Compute the inverse discrete fourier transform for the data series specified by dat. This version assumes that all terms in the data series are real. Thus, negative frequency terms correspond to the conjugate of their corresponding positive frequency terms. To reduce redundancy we only store up to the Nyquist frequency
 */
data_arr irfft(const data_arr dat) {
    //figure out how many recursive steps are needed.
    _ulong log_2_n = (_ulong)(log(dat.size) / log(2));
    _ulong truncated_size = 1 << log_2_n;

    //allocate space for the final result
    data_arr final_result;
    make_data_arr(&final_result, truncated_size/2);
    complex null_val = {0.0, 0.0};
    for (_ulong k = 0; k < final_result.size; ++k) final_result.buf[k] = null_val;

    //make a helper object that keeps track of information over recursive steps
    dat_helper help = {dat, truncated_size, 2*(_ftype)M_PI/dat.size, final_result};

    part_rfft(help, 1, 0);
    return help.sto;
}
