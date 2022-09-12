#ifndef DATA_UTILS_H
#define DATA_UTILS_H

#ifdef __cplusplus
#include <vector>
#endif

#include <stdlib.h>
#include <math.h>

#ifdef STO_PREC_32
typedef float _ftype;
#else
typedef double _ftype;
#endif

typedef unsigned long long _ulong;

// ================================ complex struct ================================
//helper functions for complex number manipulation

#ifdef __cplusplus
struct complex {
    _ftype re;
    _ftype im;

    complex& operator=(complex o) { re = o.re;im = o.im;return *this; }
    complex& operator=(_ftype o) { re = o;return *this; }
    complex operator+=(complex o) { re += o.re; im += o.im;return *this; }
    operator _ftype() const { return re; }
    complex operator+(complex o) { return { re+o.re, im+o.im }; }
    complex operator-(complex o) { return { re-o.re, im-o.im }; }
    complex operator*(complex o) { return { re*o.re - im*o.im, re*o.im + im*o.re }; }
    complex operator*(_ftype b) { return { re*b, im*b }; }
    complex operator/(complex o) { _ftype norm2 = o.re*o.re + o.im*o.im;return { (re*o.re + im*o.im)/norm2, (im*o.re - re*o.im)/norm2 }; }
    complex operator/(_ftype b) { return { re/b, im/b }; }
};
#else
typedef struct {
    _ftype re;
    _ftype im;
} complex;
#endif

//return e^val
complex conj(complex z);
complex c_exp(complex z);
complex c_mult(complex z, complex w);
complex c_multd(complex z, _ftype w);
_ftype abs_sq(complex z);
_ftype abs(complex z);

// ================================ data struct ================================
#ifdef __cplusplus
struct data_arr {
    data_arr();
    data_arr(std::vector<complex> in_pts);
    void swap(data_arr& o);
    complex& operator[](size_t i);
    data_arr(_ulong size);
    data_arr(const data_arr& o);
    data_arr& operator=(data_arr o);
    data_arr(data_arr&& o);
    ~data_arr();
#else
typedef struct {
#endif
    complex* buf;
    _ulong size;
    _ulong buf_size;

#ifdef __cplusplus
};
#else
} data_arr;
#endif

void allocate(data_arr* dat);
void make_data_arr(data_arr* dat, _ulong size);
void resize(data_arr* dat, _ulong new_size);
void cleanup_data_arr(data_arr* dat);
void add_point(data_arr* dat, complex val);
void d_add_point(data_arr* dat, _ftype val);

//piecewise arithmetic operations on data_arrs
#ifdef __cplusplus
void pw_conj(data_arr& a);
void pw_exp(data_arr& a, const data_arr b);
int pw_mult(data_arr& a, const data_arr b);
void pw_mult_scale(data_arr& a, complex z);
int pw_add(data_arr& a, const data_arr b);
void pw_add_scale(data_arr& a, complex z);
int pw_sub(data_arr& a, const data_arr b);
void pw_sub_scale(data_arr& a, complex z);
void pw_abs_sq(data_arr& a);
void pw_abs(data_arr& a);
#else
void pw_conj(data_arr a);
void pw_exp(data_arr a, const data_arr b);
int pw_mult(data_arr a, const data_arr b);
void pw_mult_scale(data_arr a, complex z);
int pw_add(data_arr a, const data_arr b);
void pw_add_scale(data_arr a, complex z);
int pw_sub(data_arr a, const data_arr b);
void pw_sub_scale(data_arr a, complex z);
void pw_abs_sq(data_arr a);
void pw_abs(data_arr a);
#endif

typedef struct {
    const data_arr dat;
    const _ulong truncated_size;
    const _ftype tau_by_n;
    data_arr sto;
} dat_helper;

/**
 * this is a helper function for fft which computes the left and right terms for the fourier transform of the data specified by dat that only takes elements of the form span*n + rem for integer n between 0 and (dat_arr.size/span - 1)
 */
#ifdef __cplusplus
void part_rfft(const dat_helper& help, _ulong span, _ulong rem);
void part_fft(const dat_helper& help, _ulong span, _ulong rem);
#else
void part_rfft(dat_helper help, _ulong span, _ulong rem);
void part_fft(dat_helper help, _ulong span, _ulong rem);
#endif
//fetch the discrete fourier transform of the data array given by arr
data_arr rfft(const data_arr dat);
data_arr irfft(const data_arr dat);
data_arr fft(const data_arr dat);
data_arr ifft(const data_arr dat);

#endif //DATA_UTILS_H
