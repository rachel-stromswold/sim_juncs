#ifndef MAIN_TEST_H
#define MAIN_TEST_H

#include <cmath>
#include <meep.hpp>
#include "argparse.h"
#include "data_utils.hpp"
#include "cgs.hpp"
#include "disp.hpp"

#define POW_N 4

#define EPSILON 0.01
#define REP_TEST_N 3
#define VAL_TEST_S  6
#define SMALL_TEST_N 10
#define TEST_N	10000

//from Numerical recipes in C
#define LCG_MULT 1664525
#define LCG_INC 1013904223
#define LCG_MOD 0xffffffff

#define TEST_SEED 31415926

typedef unsigned int _uint;
typedef unsigned char _uint8;

inline _uint lcg(_uint state) {
    return (state*LCG_MULT + LCG_INC) & LCG_MOD;
}

inline double floatize(_uint state) {
    return (double)state / LCG_MOD;
}

#endif
