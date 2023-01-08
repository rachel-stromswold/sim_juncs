#define DEBUG_INFO 1

#ifndef CGS_READ_H
#define CGS_READ_H

#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>
#include <math.h>
#include "geometry.hpp"

//hints for dynamic buffer sizes
#define BUF_SIZE 	1024
#define LINE_SIZE 	512
#define DEF_STACK_SIZE	16
#define ARGS_BUF_SIZE 	256
#define FUNC_BUF_SIZE 	256

//keywords
#define KEY_CLASS_LEN	5
#define KEY_FOR_LEN	3
#define KEY_DEF_LEN	3
#define KEY_IN_LEN	2
#define KEY_IF_LEN	2

typedef matrix<3,3> mat3x3;

typedef unsigned int _uint;
typedef unsigned char _uint8;

typedef enum { E_SUCCESS, E_NOFILE, E_LACK_TOKENS, E_BAD_TOKEN, E_BAD_SYNTAX, E_BAD_VALUE, E_BAD_TYPE, E_NOMEM, E_EMPTY_STACK, E_NOT_BINARY, E_NAN, E_NOT_DEFINED, E_OUT_OF_RANGE } parse_ercode;
typedef enum {VAL_UNDEF, VAL_STR, VAL_NUM, VAL_LIST, VAL_3VEC, VAL_MAT, VAL_FUNC, VAL_INST} valtype;
typedef enum {BLK_UNDEF, BLK_MISC, BLK_INVERT, BLK_TRANSFORM, BLK_DATA, BLK_ROOT, BLK_COMPOSITE, BLK_FUNC_DEC, BLK_LITERAL, BLK_COMMENT, BLK_SQUARE, BLK_QUOTE, BLK_QUOTE_SING, BLK_PAREN, BLK_CURLY} block_type;

/** ======================================================== utility functions ======================================================== **/

bool is_char_sep(char c);
bool is_token(const char* str, size_t i, size_t len);
char* find_token_before(char* str, size_t i);
char* strchr_block(char* str, char c);
char* token_block(char* str, const char* comp);
int read_cgs_line(char** bufptr, size_t* n, FILE* fp, size_t* lineno);
char* CGS_trim_whitespace(char* str, size_t* len);

/** ============================ line_buffer ============================ **/

class line_buffer {
private:
    char** lines;
    size_t* line_sizes;
    size_t n_lines;
public:
    line_buffer() { lines = NULL; line_sizes = NULL; n_lines = 0; }
    line_buffer(char* p_fname);
    line_buffer(const char** p_lines, size_t pn_lines);
    line_buffer(const line_buffer& o);
    line_buffer(line_buffer&& o);
    ~line_buffer();
    line_buffer get_enclosed(size_t start_line, long* end_line, char start_delim, char end_delim, size_t line_offset = 0, bool include_delims=false);
    char* get_line(size_t i) { return (i < n_lines) ? lines[i] : NULL; }
    size_t get_n_lines() { return n_lines; }
    char* flatten(char sep_char = 0);
};

/** ============================ stack ============================ **/

/**
 * A helper class (and enumeration) that can track the context and scope of a curly brace block while parsing a file
 */
template <typename T>
class stack {
protected:
    size_t stack_ptr;
    size_t buf_len;
    T* buf;
    parse_ercode grow(size_t new_size) {
	size_t old_size = buf_len;
	buf_len = new_size;

	T* tmp_buf = (T*)realloc(buf, sizeof(T)*buf_len);
	if (tmp_buf) {
	    buf = tmp_buf;
	    char* clear_buf = (char*)(tmp_buf+old_size);
	    size_t clear_len = sizeof(T)*(buf_len-old_size);
	    for (size_t i = 0; i < clear_len; ++i) clear_buf[i] = 0;
	} else {
	    buf_len = old_size;
	    return E_NOMEM;
	}

	return E_SUCCESS;
    }

public:
    stack() {
	stack_ptr = 0;
	buf_len = DEF_STACK_SIZE;
	buf = (T*)calloc(buf_len, sizeof(T));
	//check that allocation was successful
	if (!buf) {
	    buf = NULL;
	    buf_len = 0;
	    stack_ptr = 0;
	}
    }
    ~stack<T>() {
	if (buf) {
	    for (_uint i = 0; i < stack_ptr; ++i) buf[i].~T();
	    free(buf);
	    buf = NULL;
	}
	stack_ptr = 0;
	buf_len = 0;
    }
    //swap
    void swap(stack<T>& o) {
	size_t tmp = buf_len;
	buf_len = o.buf_len;
	o.buf_len = tmp;
	tmp = stack_ptr;
	stack_ptr = o.stack_ptr;
	o.stack_ptr = tmp;
	T* tmp_buf = buf;
	buf = o.buf;
	o.buf = tmp_buf;
    }
    //copy constructor, assumes that the '=' operator does something sane
    stack(const stack<T>& o) {
	stack_ptr = o.stack_ptr;
	//to save memory we'll only allocate however many entries the old object had
	buf_len = stack_ptr;
	//make sure that we set aside some space for the buffer
	if (buf_len == 0) buf_len = DEF_STACK_SIZE;
	buf = (T*)calloc(buf_len, sizeof(T));
	if (!buf) {
	    buf = NULL;
	    buf_len = 0;
	    stack_ptr = 0;
	}

	//copy the entries from the old object stack
	for (_uint i = 0; i < stack_ptr; ++i) buf[i] = o.buf[i];
    }
    //move constructor
    stack(stack<T>&& o) {
	stack_ptr = o.stack_ptr;
	buf_len = o.buf_len;
	buf = o.buf;
	//invalidate the object that just got moved
	o.stack_ptr = 0;
	o.buf_len = 0;
	o.buf = NULL;
    }
    //assignment operator
    stack<T>& operator=(stack<T> o) {
	swap(o);
	return *this;
    }
    /**
     * Push a side-object pair onto the stack.
     * returns: An error code if one occurred or E_SUCCESS if the operation was successful. It is possible for this function to fail if there wasn't sufficient space to push the object onto the stack.
     */
    parse_ercode push(T b) {
	//we might need to allocate more memory
	parse_ercode res = E_SUCCESS;
	if (stack_ptr == buf_len)
	    res = grow(2*buf_len);
	if (res == E_SUCCESS)
	    buf[stack_ptr++] = b;

	return res;
    }
    /**
     * Pop a side-object pair into the values pointed to by side and obj respectively. This function is guaranteed to succeed, but if the stack is empty then SIDE_END will be assigned to side and NULL will be assigned to obj. The caller should check that the returned values are valid.
     */
    parse_ercode pop(T* ptr) {
	//check if there is anything that can be popped
	if (stack_ptr == 0) return E_EMPTY_STACK;
	//assign the return values
	if (ptr) *ptr = buf[stack_ptr-1];
	//decrement the stack pointer
	--stack_ptr;

	return E_SUCCESS;
    }
    /**
     * Reset the stack to be empty
     */
    void reset() { stack_ptr = 0; }
    /**
     * Check wheter an entry matching key is already somewhere in the stack
     */
    bool has(const T& key) {
	for (_uint i = 0; i < stack_ptr; ++i) {
	    if (buf[i] == key) return true;
	}
	return false;
    }
    /**
     * Check wheter the stack is empty
     */
    bool is_empty() {
	return stack_ptr == 0;
    }
    /**
     * Get the number of entries on the stack
     */
    size_t size() {
	return stack_ptr;
    }
    /**
     * Returns a copy of the object at the top of the stack. If the stack is empty then the object casted from 0 is returned.
     */
    T peek(size_t ind = 1) {
	if (stack_ptr >= ind && ind > 0)
	    return buf[stack_ptr - ind];
	return (T)0;
    }
};

class CompositeObject;
/**
 * This acts similar to getline, but stops at a semicolon, newline (unless preceeded by a \), {, or }.
 * bufptr: a pointer to which the buffer is saved. If bufptr is NULL than a new buffer is allocated through malloc()
 * n: a pointer to a size_t with the number of characters in the buffer pointed to by bufptr. The call will return do nothing if n is null but *bufptr is not.
 * fp: file pointer to read from
 * linecount: a pointer to an integer specifying the number of new line characters read.
 * Returns: 0 if the end of the file was reached, 1 otherwise
 */
int read_cgs_line(char** bufptr, size_t* n, FILE* fp, size_t* line);
//
class value;
class context;
class user_func;
union V {
    char* s;
    double x;
    value* l;
    vec3* v;
    context* c;
    user_func* f;
    mat3x3* m;
};
struct value {
    valtype type;
    V val;
    size_t n_els; //only applicable for string and list types

    value() { val.x = 0;type=VAL_UNDEF;n_els=0; }
    bool operator==(std::string str);
    bool operator!=(std::string str);
    valtype get_type() { return type; }
    size_t size() { return n_els; }
    V& get_val() { return val; }
    char* to_c_str();
    double to_float();
    value cast_to(valtype type, parse_ercode& er) const;
};
void cleanup_val(value* o);
value copy_val(const value o);
void swap_val(value* a, value* b);

struct cgs_func {
    char* name;
    value args[ARGS_BUF_SIZE];
    char* arg_names[ARGS_BUF_SIZE];
    size_t n_args;
    cgs_func() { name = NULL;n_args = 0; }
};
value make_val_undef();
value make_val_num(double x);
value make_val_str(const char* s);
value make_val_std_str(std::string s);
value make_val_list(const value* vs, size_t n_vs);
value make_val_mat(mat3x3 m);
value make_val_vec3(vec3 vec);
value make_val_func(const char* name, size_t n_args, value (*p_exec)(context&, cgs_func, parse_ercode&));
cgs_func parse_func_decl(char* str);
cgs_func copy_func(const cgs_func o);
void cleanup_func(cgs_func* o);
void swap(cgs_func* a, cgs_func* b);
value lookup_named(const cgs_func f, const char* name);

//collections (like python dicts) are simply aliases for functions under the hood
typedef context collection;

/**
 * A class which stores a labeled value.
 */
class name_val_pair {
private:
    char* name;
    value val;
public:
    void swap(name_val_pair& o);
    void copy(const name_val_pair& o);
    name_val_pair();
    name_val_pair(int a);
    name_val_pair(const char* p_name, value p_val);
    name_val_pair(const name_val_pair& o);
    name_val_pair(name_val_pair&& o);
    ~name_val_pair();
    name_val_pair& operator=(name_val_pair& o);
    name_val_pair& operator=(const name_val_pair& o);
    bool name_matches(const char* str) const;
    value& get_val();
};
struct type_ind_pair {
public:
    block_type t;
    size_t i;
    type_ind_pair() { t = BLK_UNDEF;i = 0; }
    type_ind_pair(size_t ii) { t = BLK_UNDEF;i = ii; }
    type_ind_pair(block_type tt, size_t ii) { t = tt;i = ii; }
};

/** ============================ context ============================ **/
class context : public stack<name_val_pair> {
private:
    context* parent;
    value do_op(char* tok, size_t ind, parse_ercode& er);
public:
    context() : stack<name_val_pair>() { parent = NULL; }
    context(context* p_parent) : stack<name_val_pair>() { parent = p_parent; }
    //parse_ercode push(_uint side, CompositeObject* obj);
    void emplace(const char* p_name, value p_val) { name_val_pair inst(p_name, p_val);push(inst); }
    value lookup(const char* name) const;
    parse_ercode pop_n(size_t n);
    value parse_value(char* tok, parse_ercode& er);
    cgs_func parse_func(char* token, long open_par_ind, parse_ercode& f, char** end, int name_only=0);
    value parse_list(char* str, parse_ercode& sto);
    void swap(stack<name_val_pair>& o) { stack<name_val_pair>::swap(o); }
    parse_ercode set_value(const char* name, value new_val);
    parse_ercode read_from_lines(line_buffer b);

    void register_func(cgs_func sig, value (*p_exec)(context&, cgs_func, parse_ercode&));
};
/**
 * A class for functions defined by the user along with the implementation code
 */
class user_func {
private:
    cgs_func call_sig;
    line_buffer code_lines;
    value (*exec)(context&, cgs_func, parse_ercode&);
public:
    //read the function with contents stored in the file pointer fp at the current file position
    user_func(cgs_func sig, line_buffer p_buf);
    user_func(value (*p_exec)(context&, cgs_func, parse_ercode&));
    ~user_func();
    user_func(const user_func& o);
    user_func(user_func&& o);
    line_buffer& get_buffer() { return code_lines; }
    value eval(context& c, cgs_func call, parse_ercode& er);
};

#endif //CGS_READ_H
