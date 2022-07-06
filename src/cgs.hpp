#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

#include <cstring>

#define BUF_SIZE 1024
#define ARGS_BUF_SIZE 256

#define SIDE_END	2
#define SIDE_UNDEF	3

typedef unsigned int _uint;
typedef unsigned char _uint8;

typedef enum { E_SUCCESS, E_LACK_TOKENS, E_BAD_TOKEN, E_BAD_SYNTAX, E_NOMEM, E_EMPTY_STACK, E_NOT_BINARY } parse_ercode;

typedef enum { CGS_UNION, CGS_INTERSECT, CGS_DIFFERENCE } combine_type;
typedef enum { CGS_UNDEF, CGS_COMPOSITE, CGS_SPHERE, CGS_BOX, CGS_CYLINDER } object_type;

typedef Eigen::Vector3d evec3;
typedef Eigen::Vector4d evec4;

class CompositeObject;

typedef struct {
    char* name;
    char* args[ARGS_BUF_SIZE];
    _uint n_args = 0;
} cgs_func;

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

/*
 * A virtual class which describes a simple shape (such as cubes, cylinders spheres etc)
 */
class Object {
private:
    friend class GeometricObject;
protected:
    Eigen::Matrix3d trans_mat;
    int invert;

public:
    Object(int p_invert=0);
    Object(Eigen::Quaterniond& orientation, int p_invert=0);
    virtual int in(evec3& r) = 0;
    void set_inversion(int p_invert);
};

class Sphere : public Object {
private:
    evec3 center;
    double rad;

public:
    Sphere(evec3& p_center, double p_rad, int p_invert=0);
    int in(evec3& r);
};

class Box : public Object {
private:
    evec3 center;
    evec3 offset;

public:
    Box(evec3& p_corner, evec3& p_offset, int p_invert=0);
    Box(evec3& p_corner, evec3& p_offset, Eigen::Quaterniond p_orientation, int p_invert=0);
    int in(evec3& r);
};

class Cylinder : public Object {
private:
    evec3 center;
    double height;
    double r1_sq;
    double r1_sq_div_h;
    double r2_sq;

public:
    Cylinder(evec3& p_center, double p_height, double p_r1, double p_r2, int p_invert=0);
    int in(evec3& r);
};

class Scene;

/*
 * A composite object is a node in a binary tree that represents one or more primitives combined via unions and intersections
 */
class CompositeObject : public Object {
    friend class Scene;
protected:
    Object* children[2];
    object_type child_types[2];
    combine_type cmb;
    bool is_leaf;

    //metadata allows users to attach miscelaneous information to an object such as a name or material properties
    std::unordered_map<std::string, std::string> metadata;

    //since we need to perform casts to use the appropriate object type, it's helpful to add a layer of abstraction
    int call_child_in(_uint side, evec3& r);

public:
    CompositeObject(combine_type p_cmb = CGS_UNION);
    CompositeObject(combine_type p_cmb, const cgs_func& spec);
    // TODO: implement
    ~CompositeObject();
    /*CompositeObject(const CompositeObject& o);
    CompositeObject(CompositeObject&& o);*/

    //add a child to the composite object by parsing the string description
    void add_child(_uint side, Object* o, object_type p_type);
    int in(evec3& r);
    Object* get_child_l() { return children[0]; }
    Object* get_child_r() { return children[1]; }
    object_type get_child_type_l() { return child_types[0]; }
    object_type get_child_type_r() { return child_types[1]; }
    combine_type get_combine_type() { return cmb; }
};

parse_ercode parse_vector(char* str, evec3& sto);

parse_ercode make_object(const cgs_func& f, Object** ptr, object_type* type);

parse_ercode parse_func(char* token, size_t open_par_ind, cgs_func& f);

/**
 * A helper class (and enumeration) that can track the context and scope of a curly brace block while parsing a file
 */
typedef enum {BLK_UNDEF, BLK_MISC, BLK_INVERT} block_type;
class BlockStack {
private:
    size_t stack_ptr;
    size_t buf_len;
    block_type* blk_stack;
    parse_ercode grow(size_t new_size);
    void safe_alloc();
public:
    BlockStack();
    ~BlockStack();
    BlockStack(const BlockStack& o);
    BlockStack(BlockStack&& o);

    parse_ercode push(block_type b);
    block_type pop();
};

/**
 * A helper class for scene which maintains two parallel stacks used when constructing binary trees. The first specifies the integer code for the side occupied and the second stores pointers to objects. Each index specified in the side array is a "tribit" with 0 specifying the left side, 1 the right and 2 the end end of the stack
 * WARNING: The stack only handles pointers to objects. The caller is responsible for managing the lifetime of each object placed on the stack.
 */
class ObjectStack {
private:
    size_t stack_ptr;
    size_t buf_len;
    _uint* side_stack;
    CompositeObject** obj_stack;
    parse_ercode grow(size_t new_size);
    void safe_alloc();

public:
    ObjectStack();
    ~ObjectStack();
    ObjectStack(const ObjectStack& o);
    ObjectStack(ObjectStack&& o);

    parse_ercode push(_uint side, CompositeObject* obj);
    parse_ercode emplace_obj(Object* obj, object_type p_type);
    parse_ercode pop(_uint* side, CompositeObject** obj);
    _uint look_side();
    CompositeObject* look_obj();
    CompositeObject* get_root();
};

class Scene {
private:
    std::vector<Object*> objects;

public:
    Scene(const char* p_fname);
    ~Scene();

    CompositeObject* get_root();
    void read();
};
