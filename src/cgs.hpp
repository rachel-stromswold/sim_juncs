#define DEBUG_INFO 1

#ifndef CGS_H
#define CGS_H

#include "cgs_read.hpp"
#include <unordered_map>
#include <vector>
#include <cstring>
#include <cstdint>
#include <math.h>

#define DEF_TEST_N	50000
#define WALK_STEP_STC	0.05
#define WALK_STEP	0.01
#define DEF_IM_RES	255
#define IM_DEPTH	255
#define IM_LN_DEPTH	8

#define SIDE_END	2
#define SIDE_UNDEF	3

#define FLAG_NORMAL	0
#define FLAG_OPER	1
#define FLAG_COMP	2
#define FLAG_AND	3

typedef enum { CGS_UNION, CGS_INTERSECT, CGS_DIFFERENCE, CGS_CMB_NOOP } combine_type;
//note that ROOTS are a special type of COMPOSITES
typedef enum { CGS_UNDEF, CGS_ROOT, CGS_DATA, CGS_COMPOSITE, CGS_SPHERE, CGS_BOX, CGS_PLANE, CGS_CYLINDER } object_type;

/*
 * A virtual class which describes a simple shape (such as cubes, cylinders spheres etc)
 */
class object {
private:
    friend class Geometricobject;
protected:
    mat3x3 trans_mat;
    int invert;

public:
    object(int p_invert=0);
    object(const quaternion& orientation, int p_invert=0);
    virtual int in(const vec3& r) = 0;
    void set_inversion(int p_invert);
    void rescale(const vec3& components);
    void set_trans_mat(const mat3x3& new_mat) { trans_mat = new_mat; }
};

class sphere : public object {
private:
    vec3 center;
    double rad;

public:
    sphere(vec3& p_center, double p_rad, int p_invert=0);
    sphere(vec3& p_center, vec3& p_rad, int p_invert=0);
    int in(const vec3& r);

    vec3 get_center() const { return center; }
    double get_rad() const { return rad; }
};

class box : public object {
private:
    vec3 center;
    vec3 offset;

public:
    box(vec3& p_corner_1, vec3& p_corner_2, int p_invert=0);
    box(vec3& p_corner_1, vec3& p_corner_2, quaternion p_orientation, int p_invert=0);
    int in(const vec3& r);

    vec3 get_center() const { return center; }
    vec3 get_offset() const { return offset; }
};

class plane : public object {
private:
    vec3 normal;
    double offset;

public:
    plane(vec3& normal, double offset, int p_invert=0);
    plane(vec3& point_1, vec3& point_2, vec3& point_3, int p_invert=0);
    int in(const vec3& r);

    vec3 get_normal() const { return normal; }
    double get_offset() const { return offset; }
};

class cylinder : public object {
private:
    vec3 center;
    double height;
    double r1_sq;
    double r1_sq_x_h;
    double r2_sq;

public:
    cylinder(vec3& p_center, double p_height, double p_r1, double p_r2, int p_invert=0);
    int in(const vec3& r);

    vec3 get_center() const { return center; }
    double get_height() const { return height; }
    double get_r1() const { return sqrt(r1_sq); }
    double get_r2() const { return sqrt(r2_sq); }
};

class scene;
/*
 * A composite object is a node in a binary tree that represents one or more primitives combined via unions and intersections
 */
class composite_object : public object {
    friend class scene;
protected:
    object* children[2];
    object_type child_types[2];
    combine_type cmb;
    bool is_leaf;

    //metadata allows users to attach miscelaneous information to an object such as a name or material properties
    std::unordered_map<std::string, value> metadata;

    //since we need to perform casts to use the appropriate object type, it's helpful to add a layer of abstraction
    int call_child_in(_uint side, const vec3& r);
    //this is a helper function which returns a pointer to a copy of the object pointed to by the child on the specified side
    object* copy_child(_uint side) const;

public:
    composite_object(combine_type p_cmb = CGS_UNION);
    composite_object(combine_type p_cmb, const cgs_func& spec, int p_invert);
    composite_object(const composite_object& o);
    composite_object(composite_object&& o);
    composite_object& operator=(composite_object&& o);
    ~composite_object();

    //add a child to the composite object by parsing the string description
    void add_child(_uint side, object* o, object_type p_type);
    int in(const vec3& r);
    const object* get_child_l() { return children[0]; }
    const object* get_child_r() { return children[1]; }
    object_type get_child_type_l() const { return child_types[0]; }
    object_type get_child_type_r() const { return child_types[1]; }
    combine_type get_combine_type() const { return cmb; }
    int has_metadata(std::string key) const { return metadata.count(key); }
    value fetch_metadata(std::string key) { return metadata[key]; }
};

/**
 * A helper class for scene which maintains two parallel stacks used when constructing binary trees. The first specifies the integer code for the side occupied and the second stores pointers to objects. Each index specified in the side array is a "tribit" with 0 specifying the left side, 1 the right and 2 the end end of the stack
 * WARNING: The stack only handles pointers to objects. The caller is responsible for managing the lifetime of each object placed on the stack.
 */
struct side_obj_pair {
public:
    _uint side;
    composite_object* obj;
    side_obj_pair(_uint p_side, composite_object* p_obj) { side = p_side;obj=p_obj; }
    side_obj_pair(const side_obj_pair& o) { side = o.side;obj = o.obj; }
    side_obj_pair(side_obj_pair&& o) { side = o.side;obj = o.obj;o.side = 0;o.obj = NULL; }
    side_obj_pair& operator=(const side_obj_pair& o) { side = o.side;obj = o.obj;return *this; }
};
class object_stack : public stack<side_obj_pair> {
public:
    //parse_ercode push(_uint side, composite_object* obj);
    parse_ercode emplace_obj(object* obj, object_type p_type);
    parse_ercode pop(_uint* side, composite_object** obj);
    _uint look_side();
    composite_object* look_obj();
    composite_object* get_root();
};
class scene {
private:
    //std::vector<object*> objects;
    std::vector<composite_object*> roots;
    std::vector<composite_object*> data_objs;

    //let users define constants
    context named_items;
    parse_ercode fail_exit(parse_ercode er, FILE* fp);
    void save_imbuf(const char* out_fname, _uint8* z_buf, _uint8* c_buf, size_t res_x, size_t res_y);

public:
    scene() {}
    scene(const char* p_fname, parse_ercode* ercode = NULL);
    scene(const char* p_fname, context con, parse_ercode* ercode = NULL);
    scene(const scene& o);
    scene(scene&& o);
    scene& operator=(scene& o);
    ~scene();

    std::vector<composite_object*> get_roots() { return roots; }
    std::vector<composite_object*> get_data() { return data_objs; }
    void read();
 
    parse_ercode read_file(const char* p_fname);
    parse_ercode make_object(const cgs_func& f, object** ptr, object_type* type, int p_invert) const;
    parse_ercode make_transformation(const cgs_func& f, mat3x3& res) const;
    context& get_context() { return named_items; }
    //void cleanup_func(cgs_func& f);
    void draw_stochastic(const char* out_fname, vec3 cam_pos, vec3 cam_look, vec3 cam_up, size_t res=DEF_IM_RES, size_t n_samples=DEF_TEST_N);
    void draw(const char* out_fname, vec3 cam_pos, vec3 cam_look, vec3 cam_up, rvector<2> scale, size_t res_x=DEF_IM_RES, size_t res_y=DEF_IM_RES, size_t n_samples=DEF_TEST_N);
    void draw(const char* out_fname, vec3 cam_pos);
};

#endif //CGS_H
