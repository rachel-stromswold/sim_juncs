#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <meep.hpp>
#include "cgs.hpp"
#include "disp.hpp"

TEST_CASE("Test function parsing") {
    char buf[BUF_SIZE];

    const char* test_func_1 = "f()";
    const char* test_func_2 = "f(a, b, c)";
    const char* test_func_3 = "foo([0,1,2,3],a,banana)";
    const char* test_func_4 = "foo(a, Box(0,1,2,3), banana)";
    const char* test_func_5 = "foo ( a , b , c )";
    const char* bad_test_func_1 = "foo ( a , b , c";
    const char* bad_test_func_2 = "foo ( a , b , c, )";
    const char* bad_test_func_3 = "foo ( a ,, c )";

    Scene sc;
    cgs_func cur_func;
    //check string 1
    strncpy(buf, test_func_1, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, 1, cur_func, NULL);
    CHECK(cur_func.n_args == 0);
    CHECK(strcmp(cur_func.name, "f") == 0);
    //check string 2
    strncpy(buf, test_func_2, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, 1, cur_func, NULL);
    CHECK(cur_func.n_args == 3);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "f") == 0);
    INFO("func arg=", cur_func.args[0]);
    CHECK(strcmp(cur_func.args[0], "a") == 0);
    INFO("func arg=", cur_func.args[1]);
    CHECK(strcmp(cur_func.args[1], "b") == 0);
    INFO("func arg=", cur_func.args[2]);
    CHECK(strcmp(cur_func.args[2], "c") == 0);
    //check string 3
    strncpy(buf, test_func_3, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, 3, cur_func, NULL);
    CHECK(cur_func.n_args == 3);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "foo") == 0);
    INFO("func arg=", cur_func.args[0]);
    CHECK(strcmp(cur_func.args[0], "[0,1,2,3]") == 0);
    INFO("func arg=", cur_func.args[1]);
    CHECK(strcmp(cur_func.args[1], "a") == 0);
    INFO("func arg=", cur_func.args[1]);
    CHECK(strcmp(cur_func.args[2], "banana") == 0);
    //check string 4
    strncpy(buf, test_func_4, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, 3, cur_func, NULL);
    CHECK(cur_func.n_args == 3);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "foo") == 0);
    INFO("func arg=", cur_func.args[0]);
    CHECK(strcmp(cur_func.args[0], "a") == 0);
    INFO("func arg=", cur_func.args[1]);
    CHECK(strcmp(cur_func.args[1], "Box(0,1,2,3)") == 0);
    INFO("func arg=", cur_func.args[2]);
    CHECK(strcmp(cur_func.args[2], "banana") == 0);
    //check string 5
    strncpy(buf, test_func_5, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, 4, cur_func, NULL);
    CHECK(cur_func.n_args == 3);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "foo") == 0);
    INFO("func arg=", cur_func.args[0]);
    CHECK(strcmp(cur_func.args[0], "a") == 0);
    INFO("func arg=", cur_func.args[1]);
    CHECK(strcmp(cur_func.args[1], "b") == 0);
    INFO("func arg=", cur_func.args[2]);
    CHECK(strcmp(cur_func.args[2], "c") == 0);
    //check bad string 1
    strncpy(buf, bad_test_func_1, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    parse_ercode er = sc.parse_func(buf, 4, cur_func, NULL);
    CHECK(er == E_BAD_SYNTAX);
    //check bad string 2
    strncpy(buf, bad_test_func_2, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    er = sc.parse_func(buf, 4, cur_func, NULL);
    CHECK(er == E_LACK_TOKENS);
    //check bad string 3
    strncpy(buf, bad_test_func_3, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    er = sc.parse_func(buf, 4, cur_func, NULL);
    CHECK(er == E_LACK_TOKENS);
}

TEST_CASE("Test Object Trees") {
    //declare variables
    char buf[BUF_SIZE];
    cgs_func cur_func;
    ObjectStack test_stack;
    Object* cur_obj;object_type cur_type;
    parse_ercode er;

    //setup a bunch of strings describing objects
    const char* root_obj_str = "Composite(eps = 3.5)";
    const char* l_str = "union()";
    const char* ll_str = "Box([0,0,0], [1,1,1])";
    const char* lr_str = "Sphere([2,0,0], 1)";
    const char* r_str = "intersect()";
    const char* rl_str = "Box([0,0,0], [1,1,1])";
    const char* rr_str = "Cylinder([2,0,0], 1, 1)";

    Scene sc;
    //Insert root object
    strncpy(buf, root_obj_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, (size_t)(strchr(buf, '(')-buf), cur_func, NULL);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    //Insert object 1
    strncpy(buf, l_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, (size_t)(strchr(buf, '(')-buf), cur_func, NULL);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    //Insert left union object
    strncpy(buf, ll_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, (size_t)(strchr(buf, '(')-buf), cur_func, NULL);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    //Insert right union object
    strncpy(buf, lr_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, (size_t)(strchr(buf, '(')-buf), cur_func, NULL);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    //Insert object 2
    strncpy(buf, r_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, (size_t)(strchr(buf, '(')-buf), cur_func, NULL);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    //Insert left union object
    strncpy(buf, rl_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, (size_t)(strchr(buf, '(')-buf), cur_func, NULL);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    //Insert right union object
    strncpy(buf, rr_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, (size_t)(strchr(buf, '(')-buf), cur_func, NULL);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);

    //get all composite objects in the tree
    CompositeObject* root = test_stack.get_root();
    CHECK(root != NULL);
    CompositeObject* comp_l = (CompositeObject*)(root->get_child_l());
    CompositeObject* comp_r = (CompositeObject*)(root->get_child_r());
    //check that all are not NULL and that types are correct
    CHECK(root->get_child_type_l() == CGS_COMPOSITE);
    CHECK(root->get_child_type_r() == CGS_COMPOSITE);
    CHECK(comp_l != NULL);
    CHECK(comp_r != NULL);

    //now test that the composite object has the right structure
    CHECK(comp_l->get_combine_type() == CGS_UNION);
    CHECK(comp_l->get_child_l() != NULL);
    CHECK(comp_l->get_child_type_l() == CGS_BOX);
    CHECK(comp_l->get_child_r() != NULL);
    CHECK(comp_l->get_child_type_r() == CGS_SPHERE);
    //check the right branch
    CHECK(comp_r->get_combine_type() == CGS_INTERSECT);
    CHECK(comp_r->get_child_l() != NULL);
    CHECK(comp_r->get_child_type_l() == CGS_BOX);
    CHECK(comp_r->get_child_r() != NULL);
    CHECK(comp_r->get_child_type_r() == CGS_CYLINDER);
    delete root;
}

TEST_CASE("Test File Parsing") {
    Scene s("test.mt");
    std::vector<CompositeObject*> roots_vec = s.get_roots();
    CHECK(roots_vec.size() > 0);
    CompositeObject* root = roots_vec[0];

    CHECK(root != NULL);
    CompositeObject* comp_l = (CompositeObject*)(root->get_child_l());
    CompositeObject* comp_r = (CompositeObject*)(root->get_child_r());
    CompositeObject* comp_rl = (CompositeObject*)(comp_r->get_child_l());
    //check that all are not NULL and that types are correct
    CHECK(root->get_child_type_l() == CGS_COMPOSITE);
    CHECK(root->get_child_type_r() == CGS_COMPOSITE);
    CHECK(comp_r->get_child_type_l() == CGS_COMPOSITE);
    CHECK(comp_l != NULL);
    CHECK(comp_r != NULL);
    CHECK(comp_rl != NULL);

    //now test that the composite object has the right structure
    CHECK(comp_l->get_combine_type() == CGS_UNION);
    CHECK(comp_l->get_child_l() != NULL);
    CHECK(comp_l->get_child_type_l() == CGS_BOX);
    CHECK(comp_l->get_child_r() != NULL);
    CHECK(comp_l->get_child_type_r() == CGS_SPHERE);
    //check the right branch
    CHECK(comp_r->get_combine_type() == CGS_INTERSECT);
    CHECK(comp_r->get_child_l() != NULL);
    CHECK(comp_rl->get_combine_type() == CGS_DIFFERENCE);
    CHECK(comp_rl->get_child_l() != NULL);
    CHECK(comp_rl->get_child_type_l() == CGS_BOX);
    CHECK(comp_rl->get_child_r() != NULL);
    CHECK(comp_rl->get_child_type_r() == CGS_BOX);
    CHECK(comp_r->get_child_r() != NULL);
    CHECK(comp_r->get_child_type_r() == CGS_CYLINDER);
}

TEST_CASE("Test Geometric Inclusion") {
    Scene s("test.mt");
    CompositeObject* root = s.get_roots()[0];

    CHECK(root->in(Eigen::Vector3d(0.5,0.5,0.5)) == 1);
    CHECK(root->in(Eigen::Vector3d(0.5,0.1,2)) == 1);
    CHECK(root->in(Eigen::Vector3d(2.5,0.5,0.1)) == 1);
    CHECK(root->in(Eigen::Vector3d(3.1,1.1,1.1)) == 0);
    CHECK(root->in(Eigen::Vector3d(1.5,0.5,0.5)) == 0);
    CHECK(root->in(Eigen::Vector3d(1.5,0.1,4.5)) == 0);
}

TEST_CASE("Test dispersion material volumentric inclusion") {
    //load settings from the configuration file
    Settings args;
   std::string name = "params.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
    free(name_dup);
    CHECK(ret == 0);

    //create the geometry object
    Scene s("test.mt");
    CompositeObject* root = s.get_roots()[0];
    cgs_material_function mat_func(root, 2, 1);

    meep::vec test_loc_1(0.5,0.5,0.5);
    meep::vec test_loc_2(0.5,0.1,2);
    meep::vec test_loc_3(2.5,0.5,0.1);
    meep::vec test_loc_4(3.1,1.1,1.1);
    meep::vec test_loc_5(1.5,0.5,0.5);
    meep::vec test_loc_6(1.5,0.1,4.5);
    CHECK(mat_func.in_bound(test_loc_1) == 1.0);
    CHECK(mat_func.in_bound(test_loc_2) == 1.0);
    CHECK(mat_func.in_bound(test_loc_3) == 1.0);
    CHECK(mat_func.in_bound(test_loc_4) == 0.0);
    CHECK(mat_func.in_bound(test_loc_5) == 0.0);
    CHECK(mat_func.in_bound(test_loc_6) == 0.0);
}
