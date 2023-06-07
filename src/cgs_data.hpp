#ifndef CGS_DATA_HPP
#define CGS_DATA_HPP

#include <limits>
#include "cgs_read.hpp"

#define DEF_IM_RES	255
#define DEF_TEST_N	50000
#define WALK_STEP_STC	0.05
//This describes the fraction of the camera look depth that is used as a step size. Smaller values make more accurate images but take longer.
#define WALK_STEP	0.01

enum basis_comp_vectors { C_EX, C_EY, C_EZ, C_HX, C_HY, C_HZ };

basis_comp_vectors read_component_string(const char* str);
value cgs_gen_gaussian_source(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_continuous_source(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_box(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_plane(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_sphere(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_cylinder(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_composite(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_union(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_intersect(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_complement(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_difference(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_snapshot(context& c, cgs_func f, parse_ercode& er);
value cgs_gen_monitor(context& c, cgs_func f, parse_ercode& er);

void setup_geometry_context(context& con);

#endif
