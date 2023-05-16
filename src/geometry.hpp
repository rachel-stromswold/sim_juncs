#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <math.h>

#ifdef STO_PREC_32
typedef float _ftype;
#else
typedef double _ftype;
#endif

// ================================ matrix struct ================================

template <size_t M, size_t N>
struct matrix {
public:
    _ftype el[M*N];
    //Default constructor to zero matrix
    matrix(int val=0) {
	for (size_t i = 0; i < M*N; ++i) {
	    el[i] = 0;
	}
	if (val != 0) {
	    size_t min_m = M;
	    if (N < M) min_m = N;
	    for (size_t i = 0; i < min_m; ++i) {
		el[N*i+i] = val;
	    }
	}
    }
    /**
     * Set the row with index i to be specified by the matrix r
     */
    void set_row(size_t i, matrix<N,1> r) {
	if (i < M) {
	    for (size_t j = 0; j < N; ++j) {
		el[N*i+j] = r.el[j];
	    }
	}
    }
    /**
     * Set the row with index i to be specified by the matrix r
     */
    void set_col(size_t j, matrix<N,1> c) {
	if (j < N) {
	    for (size_t i = 0; i < N; ++i) {
		el[N*i+j] = c.el[i];
	    }
	}
    }
    _ftype get(size_t i, size_t j) {
	return el[N*i + j];
    }
    //copy constructor
    matrix(const matrix<M,N>& o) {
	for (size_t i = 0; i < M*N; ++i) el[i] = o.el[i];
    }
    //assignment operator
    matrix<M,N>& operator=(const matrix<M,N>& o) {
	for (size_t i = 0; i < M*N; ++i) {
	    el[i] = o.el[i];
	}
	return *this;
    }
    //scalar operations
    matrix<M,N> operator*(const _ftype s) const {
	matrix<M,N> ret;
	for (size_t i = 0; i < M*N; ++i) { ret.el[i] = s*el[i]; }
	return ret;
    }
    matrix<M,N> operator/(const _ftype s) const {
	matrix<M,N> ret;
	for (size_t i = 0; i < M*N; ++i) { ret.el[i] = el[i]/s; }
	return ret;
    }
    //add two matrices
    matrix<M,N> operator+(const matrix<M,N>& o) const {
	matrix<M,N> ret;
	for (size_t i = 0; i < M*N; ++i) { ret.el[i] = el[i] + o.el[i]; }
	return ret;
    }
    matrix<M,N> operator-(const matrix<M,N>& o) const {
	matrix<M,N> ret;
	for (size_t i = 0; i < M*N; ++i) { ret.el[i] = el[i] - o.el[i]; }
	return ret;
    }
    //multiply two matrices
    template <size_t L>
    matrix<M,L> operator*(const matrix<N,L> o) const {
	matrix<M,L> ret;
	for (size_t i = 0; i < M; ++i) {
	    for (size_t j = 0; j < L; ++j) {
		size_t ind = L*i+j;
		ret.el[ind] = 0;
		for (size_t k = 0; k < N; ++k) {
		    ret.el[ind] += el[N*i+k]*o.el[L*k+j];
		}
	    }
	}
	return ret;
    }
    //comparisons
    bool operator==(const matrix<M,N>& o) const {
	for (size_t i = 0; i < M*N; ++i) {
	    if (el[i] != o.el[i]) return false;
	}
	return true;
    }
    static matrix<M,N> id() {
	matrix<M,N> ret;
	//partial identity for non square matrices
	size_t min_ind = M;
	if (N < min_ind) min_ind = N;
	for (size_t i = 0; i < min_ind; ++i) {
	    ret.el[i*N+i] = 1.0;
	}
	return ret;
    }
    size_t to_str(char* str, size_t n) {
	//we subtract 2*M characters for the open and close brackets for each row, M-1 for the commas separating rows, and M*(N-1) for the commas in columns
	//(n/MN) - 2M - (M-1) - M(N-1) - 2 = (n/MN) - 2M - MN -1
	size_t c_per_el = (n / (M*N)) - 2*M - M*N - 1;
	if (str && c_per_el > 0 && n > 2) {
	    str[0] = '[';str[1] = '[';//]]
	    size_t off = 2;
	    for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
		    //make sure we don't go past the end
		    if (off >= n-2) {
			str[n-1] = 0;
			return off;
		    }
		    //write the row with a terminating ',' or not depending on whether we're at the end
		    if (j < N-1)
			off += snprintf(str+off, c_per_el, "%f,", el[N*i+j]);
		    else
			off += snprintf(str+off, c_per_el, "%f", el[N*i+j]);
		}
		if (i < M-1) str[off++] = ',';
	    }
	    return off;
	}
	if (str && n > 0) str[0] = 0;
	return 0;
    }
};
template <size_t M,size_t N> inline matrix<M,N> operator*(_ftype s, const matrix<M,N>& m) { return m*s; }

template <size_t N>
struct rvector : public matrix<N,1> {
    rvector<N>() { for(size_t i = 0; i < N; ++i) this->el[i] = 0; }
    /**
     * freely enable conversions between row and column rvectors
     */
    rvector<N>(const matrix<1,N>& o) {
	for (size_t i = 0; i < N; ++i) this->el[i] = o.el[i];
    }
    rvector<N>& operator=(const matrix<1,N>& o) {
	for (size_t i = 0; i < N; ++i) this->el[i] = o.el[i];
	return this;
    }
    /**
     * create an NxN matrix with elements from the rvector along the diagonal
     */
    matrix<N,N> make_diagonal() const {
	matrix<N,N> ret;
	for (size_t i = 0; i < N; ++i) {
	    ret.el[i*N + i] = this->el[i];
	}
	return ret;
    }
    //return the square of the normal of the rvector
    _ftype normsq() {
	_ftype ret = 0.0;
	for (size_t i = 0; i < N; ++i) ret += (this->el[i])*(this->el[i]);
	return ret;
    }
    //return the normal of the rvector
    _ftype norm() {
	return sqrt(normsq());
    }
    //return a rvector which points along the same direction with unit length
    rvector<N> normalize() {
	rvector<N> ret;
	_ftype n = norm();
	for (size_t i = 0; i < N; ++i) ret.el[i] = this->el[i]/n;
	return ret;
    }
};
template <size_t N> inline rvector<N> operator*(_ftype s, const rvector<N>& m) { return m*s; }

struct vec3 : public rvector<3> {
    vec3() { this->el[0] = 0;this->el[1] = 0;this->el[2] = 0; }
    vec3(const matrix<3,1>& o) { el[0] = o.el[0];el[1] = o.el[1];el[2] = o.el[2]; }
    vec3& operator=(const matrix<3,1>& o) { el[0] = o.el[0];el[1] = o.el[1];el[2] = o.el[2];return *this; }
    //construct based on coordinates
    vec3(_ftype x, _ftype y, _ftype z) {
	el[0] = x;el[1] = y;el[2] = z;
    }
    _ftype dot(const vec3& b) const {
	return el[0]*b.el[0] + el[1]*b.el[1] + el[2]*b.el[2];
    }
    vec3 cross(const vec3& b) const {
	vec3 ret;
	ret.el[0] = el[1]*b.el[2] - el[2]*b.el[1];
	ret.el[1] = el[2]*b.el[0] - el[0]*b.el[2];
	ret.el[2] = el[0]*b.el[1] - el[1]*b.el[0];
	return ret;
    }
    _ftype& x() { return el[0]; }
    _ftype& y() { return el[1]; }
    _ftype& z() { return el[2]; }
};
inline vec3 operator*(_ftype s, const vec3& m) { return m*s; }
//make a rotation matrix about the specified axis with the specified rotation
inline matrix<3,3> make_rotation(_ftype theta, vec3 a) {
    a = a.normalize();
    matrix<3,3> ret;
    _ftype ct = cos(theta);
    _ftype st = sin(theta);
    _ftype ctc = 1.0-ct;
    ret.el[0] = ct+a.el[0]*a.el[0]*ctc;	ret.el[1] = a.el[0]*a.el[1]*ctc-a.el[2]*st;	ret.el[2] = a.el[0]*a.el[2]*ctc + a.el[1]*st;
    ret.el[3] = ret.el[1]+2*a.el[2]*st;	ret.el[4] = ct+a.el[1]*a.el[1]*ctc;		ret.el[5] = a.el[1]*a.el[2]*ctc - a.el[0]*st;
    ret.el[6] = ret.el[2]-2*a.el[1]*st;	ret.el[7] = ret.el[5]+2*a.el[0]*st;		ret.el[8] = ct+a.el[2]*a.el[2]*ctc;
    return ret;
}
//TODO: finish
inline matrix<3,3> invert(matrix<3,3> a) {
    _ftype inv_det = a.get(0,0)*(a.get(1,1)*a.get(2,2) - a.get(1,2)*a.get(2,1))
		    + a.get(0,1)*(a.get(1,2)*a.get(2,0) - a.get(1,0)*a.get(2,2))
		    + a.get(0,2)*(a.get(1,0)*a.get(2,1) - a.get(1,1)*a.get(2,0));
    matrix<3,3> ret;
    return ret;
}

struct quaternion : public rvector<4> {
    quaternion() : rvector<4>() {}
    quaternion operator*(const quaternion& o) {
	quaternion ret;
	ret.el[0] = el[0]*o.el[0] - el[1]*o.el[1] - el[2]*o.el[2] - el[3]*o.el[3];
	ret.el[1] = el[2]*o.el[3] - el[3]*o.el[2];
	ret.el[2] = el[3]*o.el[1] - el[1]*o.el[3];
	ret.el[3] = el[1]*o.el[2] - el[2]*o.el[1];
	return ret;
    }
    quaternion operator*(const vec3& o) {
	quaternion ret;
	ret.el[0] = -(el[1]*o.el[0] + el[2]*o.el[1] + el[3]*o.el[2]);
	ret.el[1] = el[0]*o.el[1] + el[2]*o.el[2] - el[3]*o.el[1];
	ret.el[2] = el[0]*o.el[2] + el[3]*o.el[0] - el[1]*o.el[2];
	ret.el[3] = el[0]*o.el[3] + el[1]*o.el[1] - el[2]*o.el[0];
	return ret;
    }
};

#endif
