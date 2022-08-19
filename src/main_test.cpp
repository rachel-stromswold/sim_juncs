#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>
#include "main_test.hpp"

TEST_CASE("Test Fourier transforms") {
    data_arr dat;
    make_data_arr(&dat, 1 << POW_N);
    for (size_t k = 0; k < dat.size; ++k) {
	dat.buf[k].re = sin(4*M_PI*k / dat.size);
	dat.buf[k].im = 0.0; 
    }

    //fourier transform the data
    data_arr fft_dat = rfft(dat);
    CHECK(fft_dat.size == dat.size/2);
    CHECK(fft_dat[2].re == doctest::Approx(0.0));
    CHECK(fft_dat[2].im == doctest::Approx(-8.0));

    //print out the fourier transform for debugging purposes
    printf("Fourier transform of dat[k]=sin(2*pi*2*k/size): \n");
    for (size_t k = 0; k < fft_dat.size; ++k) printf("%f+%fi ", fft_dat[k].re, fft_dat[k].im);
    printf("\n");

    //Only the magnitude is of interest to us for testing
    pw_abs(fft_dat);

    //since we initialize the raw data with dat[k]=sin(2*pi*2*k/size), we expect the second Fourier component to be very large compared to everything else.
    CHECK(fft_dat.buf[2] > fft_dat.buf[0]);
    CHECK(fft_dat.buf[2] > fft_dat.buf[1]);
    for (size_t k = 3; k < fft_dat.size; ++k) {
	CHECK(fft_dat.buf[2] > fft_dat.buf[k]);
    }

    //perform initial cleanup
    cleanup_data_arr(&fft_dat);

    //check that the fourier transform of a gaussian is a gaussian
    const double gauss_sigma_sq = 2;
    const int k_0 = 8;//use ints to make k-k_0 and k_0-k both valid. Totally not sketch :p
    for (int k = 0; k < dat.size; ++k) {
	//fill up each index with a pseudo random number between 0 and 1
	dat.buf[k].re = exp( (k-k_0)*(k_0-k)/(2*gauss_sigma_sq) );
	dat.buf[k].im = 0.0; 
    }
    //print out the result
    printf("Input Gaussian: \n");
    for (size_t k = 0; k < dat.size; ++k) printf("%f+%fi ", dat[k].re, dat[k].im);
    printf("\n");

    //take the fourier transform and multiply by a normalization
    fft_dat = rfft(dat);
    complex scale = {1/sqrt(2*M_PI*gauss_sigma_sq), 0};
    pw_mult_scale(fft_dat, scale);
    //print out the result
    printf("Fourier transform Gaussian: \n");
    for (size_t k = 0; k < fft_dat.size; ++k) printf("%f+%fi ", fft_dat[k].re, fft_dat[k].im);
    printf("\n");

    //check that the inverse fft of the fft is approximately the original
    double omega_scale = 2*M_PI/dat.size;
    double k_0_by_sigma = k_0/sqrt(gauss_sigma_sq);
    for (size_t k = 0; k < fft_dat.size; ++k) {
	double omega = omega_scale*k;
	//the phase should be -\frac{\sigma^2\omega^2}{2} - i\omega k_0
	complex expected_phase = {-omega*omega*gauss_sigma_sq/2, -omega*k_0};
	complex expected = c_exp(expected_phase);
	CHECK(fft_dat[k].re == doctest::Approx(expected.re));
	CHECK(fft_dat[k].im == doctest::Approx(expected.im));
    }

    cleanup_data_arr(&dat);
}

TEST_CASE("Test function parsing") {
    char buf[BUF_SIZE];

    const char* test_func_1 = "f()";
    const char* test_func_2 = "f(a, b, c)";
    const char* test_func_3 = "foo([0,1,2,3],a,banana)";
    const char* test_func_4 = "foo(a, Box(0,1,2,3), banana)";
    const char* test_func_5 = "foo ( a , b , c )";
    const char* test_func_6 = "f(eps = 3.5)";
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
    //check string 6
    strncpy(buf, test_func_6, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    sc.parse_func(buf, 1, cur_func, NULL);
    CHECK(cur_func.n_args == 1);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "f") == 0);
    INFO("func arg=", cur_func.args[0]);
    CHECK(strcmp(cur_func.args[0], "eps = 3.5") == 0);
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
    parse_ercode er;
    Scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    //check that metadata works
    std::vector<CompositeObject*> data_vec = s.get_data();
    CHECK(data_vec.size() > 0);
    CHECK(data_vec[0]->has_metadata("name"));
    CHECK(data_vec[0]->has_metadata("entry"));
    std::string name_str = data_vec[0]->fetch_metadata("name");
    std::string entry_str = data_vec[0]->fetch_metadata("entry");
    CHECK(name_str == "foo");
    CHECK(entry_str == "bar,(arr),[blah]");

    //test geometric information
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
    parse_ercode er;
    Scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    CompositeObject* root = s.get_roots()[0];

    CHECK(root->in(Eigen::Vector3d(4.5,4.5,4.5)) == 1);
    CHECK(root->in(Eigen::Vector3d(4.5,4.1,6)) == 1);
    CHECK(root->in(Eigen::Vector3d(6.5,4.5,4.1)) == 1);
    CHECK(root->in(Eigen::Vector3d(7.1,5.1,5.1)) == 0);
    CHECK(root->in(Eigen::Vector3d(5.5,4.5,4.5)) == 0);
    CHECK(root->in(Eigen::Vector3d(5.5,4.1,8.5)) == 0);
}

TEST_CASE("Test dispersion material volumentric inclusion") {
    parse_ercode er;
    //load settings from the configuration file
    Settings args;
    std::string name = "tests/test.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
    free(name_dup);
    CHECK(ret == 0);

    //create the geometry object
    Scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    CompositeObject* root = s.get_roots()[0];
    cgs_material_function mat_func(root);

    //check locations
    meep::vec test_loc_1(4.5,4.5,4.5);
    meep::vec test_loc_2(4.5,4.1,6);
    meep::vec test_loc_3(6.5,4.5,4.1);
    meep::vec test_loc_4(7.1,5.1,5.1);
    meep::vec test_loc_5(5.5,4.5,4.5);
    meep::vec test_loc_6(5.5,4.1,8.5);
    CHECK(mat_func.in_bound(test_loc_1) == 3.5);
    CHECK(mat_func.in_bound(test_loc_2) == 3.5);
    CHECK(mat_func.in_bound(test_loc_3) == 3.5);
    CHECK(mat_func.in_bound(test_loc_4) == 1.0);
    CHECK(mat_func.in_bound(test_loc_5) == 1.0);
    CHECK(mat_func.in_bound(test_loc_6) == 1.0);
}

TEST_CASE("Test reading of configuration files") {
    std::string name = "tests/test.conf";
    char* name_dup = strdup(name.c_str());

    SUBCASE("Reading just a config file works") {
	Settings args;
	int ret = parse_conf_file(&args, name_dup);

	CHECK(args.n_dims == 3);
	CHECK(args.pml_thickness == 2.0);
	CHECK(args.len == 16.0);
	CHECK(args.um_scale == 1.0);
	//NOTE that resolutions are rounded up to accomodate the nearest odd integer number of grid points. This value should NOT be the same as that given in the .conf file.
	CHECK(args.resolution == 1.55);
	CHECK(args.courant == 0.3);
	CHECK(args.smooth_n == 1);
	CHECK(args.smooth_rad == 0.25);
	CHECK(strcmp(args.monitor_locs, "(1.0,1.0,1.0)") == 0);
	CHECK(args.post_source_t == 1.0);
	CHECK(args.ambient_eps == 1.0);
	CHECK(strcmp(args.geom_fname, "tests/test.geom") == 0);
    }

    SUBCASE("Command line arguments override defaults") {
	Settings args;
	//parse commandline arguments
	std::string sim_argv_cpp[] = { "./test", "--conf-file", "blah.conf", "--geom-file", "blah.geom", "--out-dir", "/blah", "--grid-res", "3.0", "--length", "9.0", "--eps1", "2.0" };
	//gcc doesn't like using string literals as c-strings, ugh
	int n_args = sizeof(sim_argv_cpp)/sizeof(std::string);
	char** sim_argv_c = (char**)malloc(sizeof(char*)*n_args);
	for (_uint i = 0; i < n_args; ++i) sim_argv_c[i] = strdup(sim_argv_cpp[i].c_str());

	//finally we can parse the command line arguments
	int ret = parse_args(&args, &n_args, sim_argv_c);

	//this is used when calling parse_args, so it should be checked before everything else
	CHECK(args.conf_fname != NULL);
	CHECK(strcmp(args.conf_fname, "blah.conf") == 0);

	//read the config file
	ret = parse_conf_file(&args, name_dup);
	CHECK(args.geom_fname != NULL);
	CHECK(strcmp(args.geom_fname, "blah.geom") == 0);
	CHECK(args.out_dir != NULL);
	CHECK(strcmp(args.out_dir, "/blah") == 0);
	CHECK(args.len == 9.0);
	//NOTE that resolutions are rounded up to accomodate the nearest odd integer number of grid points. This value should NOT be the same as that given in the command line arguments.
	CHECK(args.resolution == doctest::Approx(3.1538));
	CHECK(args.ambient_eps == 2.0);
	CHECK(args.n_dims == 3);
	CHECK(args.pml_thickness == 2.0);
	CHECK(args.um_scale == 1.0);
	//NOTE that resolutions are rounded up to accomodate the nearest odd integer number of grid points. This value should NOT be the same as that given in the .conf file.
	CHECK(args.courant == 0.3);
	CHECK(args.smooth_n == 1);
	CHECK(args.smooth_rad == 0.25);
	CHECK(strcmp(args.monitor_locs, "(1.0,1.0,1.0)") == 0);
	CHECK(args.post_source_t == 1.0);

	//deallocate memory
	for (_uint i = 0; i < n_args; ++i) free(sim_argv_c[i]);
	free(sim_argv_c);
    }

    free(name_dup);
}

TEST_CASE("Test geometry file reading") {
    parse_ercode er;
    //load settings from the configuration file
    Settings args;
    std::string name = "tests/test.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
    free(name_dup);
    CHECK(ret == 0);

    //make things a little faster because we don't care
    args.pml_thickness = 1.0;
    args.len = 2.0;

    //check the susceptibilities
    bound_geom geometry(args, &er);

    SUBCASE("Test reading of susceptibilities") {
	CHECK(er == E_SUCCESS);
	CompositeObject* root = geometry.problem.get_roots()[0];
	char* dat = strdup(root->fetch_metadata("susceptibilities").c_str());
	std::vector<drude_suscept> sups = geometry.parse_susceptibilities(dat, (int*)(&er));
	free(dat);
	CHECK(sups.size() == 2);
	CHECK(sups[0].omega_0 == 1.0);
	CHECK(sups[0].gamma == 0.48);
	CHECK(sups[0].sigma == 68.5971845);
	CHECK(sups[0].use_denom == false);
	CHECK(sups[1].omega_0 == 8.0);
	CHECK(sups[1].gamma == 0.816);
	CHECK(sups[1].sigma == 452848600800781300);
	CHECK(sups[1].use_denom == false);

	cleanup_settings(&args);
    }

    SUBCASE("Test reading of field sources") {
#ifdef DEBUG_INFO
	CHECK(geometry.sources.size() == 2);
	source_info inf = geometry.sources[0];
	CHECK(inf.type == SRC_GAUSSIAN);
	CHECK(inf.component == meep::Ey);
	CHECK(inf.freq == 1.33);
	CHECK(inf.width == 3.0);
	CHECK(inf.start_time == 0.2);
	CHECK(inf.cutoff == 5.0);
	CHECK(inf.amplitude == 1.0);

	inf = geometry.sources[1];
	CHECK(inf.type == SRC_CONTINUOUS);
	CHECK(inf.component == meep::Hz);
	CHECK(inf.freq == 1.66);
	CHECK(inf.start_time == 0.2);
	CHECK(inf.end_time == 1.2);
	CHECK(inf.width == 0.1);
	CHECK(inf.amplitude == 1.0);
#endif

	cleanup_settings(&args);
    }
}

complex* read_h5_array_raw(const H5::H5File& file, const std::string name, size_t* n_entries) {
    char tname[BUF_SIZE];
    //We need to create an hdf5 data type for complex values
    H5::CompType fieldtype(sizeof(complex));
    //for some reason linking insertMember breaks on the cluster, we do it manually
    hid_t float_member_id = H5::PredType::NATIVE_FLOAT.getId();
    snprintf(tname, BUF_SIZE, "Re");
    herr_t ret_val = H5Tinsert(fieldtype.getId(), tname, HOFFSET(complex, re), float_member_id);
    snprintf(tname, BUF_SIZE, "Im");
    ret_val = H5Tinsert(fieldtype.getId(), tname, HOFFSET(complex, im), float_member_id);

    *n_entries = 0;
    try {
	//find the dataspace for real values
	H5::DataSet dataset = file.openDataSet(name);
	H5::DataType datatype = dataset.getDataType();
	H5::DataSpace dataspace = dataset.getSpace();
	size_t n_pts = dataspace.getSimpleExtentNpoints();
	//allocate memory for storage
	complex* data = new complex[n_pts];
	dataset.read(data, fieldtype);

	*n_entries = n_pts;
	dataspace.close();
	datatype.close();
	dataset.close();
	return data;
    } catch (H5::FileIException error) {
	error.printErrorStack();
	return NULL;
    } catch (H5::GroupIException error) {
	error.printErrorStack();
	return NULL;
    }
}

TEST_CASE("Test running with a very small system") {
    parse_ercode er;

    //load settings from the configuration file
    Settings args;
    std::string name = "tests/run.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
    free(name_dup);
    CHECK(ret == 0);
    
    //make things a little faster because we don't care
    args.pml_thickness = 1.0;
    args.len = 2.0;
    args.resolution = 2.0;

    //try creating the geometry object
    bound_geom geometry(args, &er);
    CHECK(er == E_SUCCESS);

    //make sure that monitor locations were added
    CHECK(geometry.get_monitor_locs().size() > 0);
    geometry.run("/tmp");
    //fetch the field times
    std::vector<data_arr> field_times = geometry.get_field_times();
    CHECK(field_times.size() > 0);

    //check that writing hdf5 files works
    geometry.save_field_times("/tmp");

    //read the h5 file
    //H5::H5std_string fname("/tmp/fields_0.txt");
    H5::H5File file("/tmp/fields_0.h5", H5F_ACC_RDONLY);
    //read the data from the file we opened
    size_t n_f_pts, n_t_pts;
    complex* f_data = read_h5_array_raw(file, "frequency", &n_f_pts);
    complex* t_data = read_h5_array_raw(file, "time", &n_t_pts);
    CHECK(f_data != NULL);
    CHECK(t_data != NULL);
    CHECK(n_f_pts > 0);
    //since the fourier transform should only go to +- the nyquist frequency, it must have fewer elements
    CHECK(n_t_pts >= n_f_pts);
    delete[] f_data;
    delete[] t_data;
    file.close();

    cleanup_settings(&args);
}

int main(int argc, char** argv) {
    meep::initialize mpi(argc, argv);

    doctest::Context context;
    context.applyCommandLine(argc, argv);

    // overrides
    context.setOption("no-breaks", true);             // don't break in the debugger when assertions fail

    int res = context.run(); // run

    if(context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests
    
    int client_stuff_return_code = 0;
    // your program - if the testing framework is integrated in your production code
    
    return res + client_stuff_return_code; // the result from doctest is propagated here as well
}
