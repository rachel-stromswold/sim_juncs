#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>
#include "main_test.hpp"

TEST_CASE("make_dec_str") {
    char buf[BUF_SIZE];
    size_t res = make_dec_str(buf, BUF_SIZE, 0.25, 1, 2, '.');
    CHECK(res == 4);
    CHECK(strcmp(buf, "0.25") == 0);
    res = make_dec_str(buf, BUF_SIZE, 0.25, 2, 2, '_');
    CHECK(res == 5);
    CHECK(strcmp(buf, "00_25") == 0);
    res = make_dec_str(buf, 3, 0.25, 1, 2);
    CHECK(res == -1);
}

TEST_CASE("Fourier transforms") {
    size_t arr_size = 1 << POW_N;

    SUBCASE("a very simple sine wave") {
	data_arr dat;
	make_data_arr_zeros(&dat, 1);
	//we initialize the array to have size 1, so we must initialize the first point
	for (size_t k = 1; k < arr_size; ++k) {
	    complex tmp(sin(4*M_PI*k / arr_size), 0.0);
	    add_point(&dat, tmp);
	}
	CHECK(dat.size == arr_size);

	//fourier transform the data
	data_arr fft_dat = rfft(dat);
	CHECK(fft_dat.size == dat.size/2);
	CHECK(fft_dat[2].re == doctest::Approx(0.0));
	CHECK(fft_dat[2].im == doctest::Approx(-8.0));

	//print out the fourier transform for debugging purposes
	printf("Fourier transform of dat[k]=sin(2*pi*2*k/size): \n");
	for (size_t k = 0; k < fft_dat.size; ++k) printf("%f+%fi ", fft_dat[k].re, fft_dat[k].im);
	printf("\n");

	//take the inverse fourier transform and make sure it's close

	//Only the magnitude is of interest to us for testing
	pw_abs(fft_dat);

	//since we initialize the raw data with dat[k]=sin(2*pi*2*k/size), we expect the second Fourier component to be very large compared to everything else.
	CHECK(abs(fft_dat.buf[2]) > abs(fft_dat.buf[0]));
	CHECK(abs(fft_dat.buf[2]) > abs(fft_dat.buf[1]));
	for (size_t k = 3; k < fft_dat.size; ++k) {
	    CHECK(abs(fft_dat.buf[2]) > abs(fft_dat.buf[k]));
	}

	//perform initial cleanup
	cleanup_data_arr(&fft_dat);
	cleanup_data_arr(&dat);
	CHECK(fft_dat.buf_size == 0);
	CHECK(fft_dat.size == 0);
	CHECK(fft_dat.buf == NULL);
    }

    SUBCASE("the fourier transform of a gaussian is a gaussian") {
	data_arr dat;
	make_data_arr(&dat, arr_size);
	//initialize the gaussian
	const _ftype gauss_sigma_sq = 2;
	const int k_0 = 8;
	for (_uint k = 0; k < arr_size; ++k) {
	    //fill up each index with a pseudo random number between 0 and 1
	    dat.buf[k].re = exp( ((int)k-k_0)*(k_0-(int)k)/(2*gauss_sigma_sq) );
	    dat.buf[k].im = 0.0; 
	}
	dat.size = arr_size;
	//print out the result
	printf("Input Gaussian: \n");
	for (size_t k = 0; k < dat.size; ++k) printf("%f+%fi ", dat[k].re, dat[k].im);
	printf("\n");

	//take the fourier transform and multiply by a normalization
	data_arr rfft_dat = rfft(dat);
	data_arr fft_dat = fft(dat);
	complex scale = {1/(_ftype)sqrt(2*M_PI*gauss_sigma_sq), 0};
	pw_mult_scale(rfft_dat, scale);
	pw_mult_scale(fft_dat, scale);
	//check sizes, the full fourier transform should have twice as many data points
	CHECK(rfft_dat.size == dat.size/2);
	CHECK(fft_dat.size == dat.size);

	//check that multiplying by the conjugate gives the square absolute value
	data_arr fft_dat_conj(fft_dat);
	pw_conj(fft_dat_conj);
	pw_mult(fft_dat_conj, fft_dat);
	data_arr fft_dat_abs(fft_dat);
	pw_abs_sq(fft_dat_abs);
	for (size_t k = 0; k < fft_dat_conj.size; ++k) {
	    CHECK(fft_dat_conj.buf[k].re == doctest::Approx(fft_dat_abs.buf[k].re));
	    CHECK(fft_dat_conj.buf[k].im == doctest::Approx(fft_dat_abs.buf[k].im));
	}
	cleanup_data_arr(&fft_dat_conj);
	cleanup_data_arr(&fft_dat_abs);

	//print out the result
	printf("Fourier transform Gaussian: \n");
	for (size_t k = 0; k < rfft_dat.size; ++k) printf("%f+%fi ", rfft_dat[k].re, rfft_dat[k].im);
	printf("\n");

	//the fourier transform should be a gaussian with width 1/sigma
	_ftype omega_scale = 2*M_PI/dat.size;
	for (size_t k = 0; k < rfft_dat.size; ++k) {
	    _ftype omega = omega_scale*k;
	    //the phase should be -\frac{\sigma^2\omega^2}{2} - i\omega k_0
	    complex expected_phase = {-omega*omega*gauss_sigma_sq/2, -omega*k_0};
	    complex expected = c_exp(expected_phase);
	    CHECK(rfft_dat[k].re == doctest::Approx(expected.re));
	    CHECK(rfft_dat[k].im == doctest::Approx(expected.im));
	    //Since frequency wraps around at zero, the zero frequency part is only stored once, at fft_dat[0]. 
	    if (k > 0) {
		complex tmp = conj(fft_dat[k]);
		CHECK( tmp.re == doctest::Approx(fft_dat[fft_dat.size-k].re) );
		CHECK( tmp.im == doctest::Approx(fft_dat[fft_dat.size-k].im) );
	    }
	}
	cleanup_data_arr(&rfft_dat);
	cleanup_data_arr(&fft_dat);
	cleanup_data_arr(&dat);
    }

    SUBCASE("the inverse fourier transform of a fourier transform is itself") {
	//initialize the lcg
	_uint rand = 314159265;
	complex phi = {0.0, 0.0};
	complex x = {0.0, 0.0};
	//generate a pseudo-random walk
	std::vector<complex> tmp_pts;
	for (size_t k = 0; k < arr_size; ++k) {
	    tmp_pts.push_back(x);
	    //fill up each index with a pseudo random number between 0 and 1
	    //dat.buf[k].re = x.re;
	    //dat.buf[k].im = x.im;
	    rand = lcg(rand);
	    phi.im = 2*M_PI*floatize(rand);
	    x += c_exp(phi);
	}
	data_arr dat(tmp_pts);
	CHECK(dat.size == arr_size);
	_ftype N = (_ftype)dat.size;
	//subtract off constant part and print out the result
	printf("Input walk: \n");
	for (size_t k = 0; k < dat.size; ++k) {
	    printf("%f+%fi ", dat[k].re, dat[k].im);
	}
	printf("\n");

	//take the fourier transform
	data_arr fft_dat = fft(dat);
	//print out the result
	printf("Fourier transform of walk: \n");
	for (size_t k = 0; k < fft_dat.size; ++k) printf("%f+%fi ", fft_dat[k].re, fft_dat[k].im);
	printf("\n");
	//take the inverse fourier transform and scale down by a factor of N-1 (note that ifft(fft(dat)) picks up a factor N-1)
	data_arr new_dat = ifft(fft_dat);
	complex scale = {1/N, 0};
	pw_mult_scale(new_dat, scale);
	//print out the result
	printf("inverse Fourier transform of walk: \n");
	for (size_t k = 0; k < fft_dat.size; ++k) printf("%f+%fi ", new_dat[k].re, new_dat[k].im);
	printf("\n");

	//make sure the sizes are the same and compare elements
	CHECK(dat.size == new_dat.size);
	for (size_t k = 0; k < fft_dat.size; ++k) {
	    CHECK(dat[k].re == doctest::Approx(new_dat[k].re));
	    CHECK(dat[k].im == doctest::Approx(new_dat[k].im));
	}
	cleanup_data_arr(&fft_dat);
	cleanup_data_arr(&new_dat);
	cleanup_data_arr(&dat);
    }
}

TEST_CASE("geometry operations") {
    _uint state = lcg(lcg(TEST_SEED));
    vec3 x_axis(1,0,0);vec3 y_axis(0,1,0);vec3 z_axis(0,0,1);
    //orthogonality and right hand rule
    CHECK(x_axis.dot(y_axis) == 0);
    CHECK(y_axis.dot(z_axis) == 0);
    CHECK(z_axis.dot(x_axis) == 0);
    CHECK(x_axis.cross(y_axis) == z_axis);
    CHECK(y_axis.cross(z_axis) == x_axis);
    CHECK(x_axis.cross(z_axis) == y_axis*-1);
    SUBCASE("cross product identities") {
	double x,y,z;
	for (size_t i = 0; i < SMALL_TEST_N; ++i) {
	    x = floatize(state);state = lcg(state);
	    y = floatize(state);state = lcg(state);
	    z = floatize(state);state = lcg(state);
	    vec3 aa(x, y, z);
	    x = floatize(state);state = lcg(state);
	    y = floatize(state);state = lcg(state);
	    z = floatize(state);state = lcg(state);
	    vec3 bb(x, y, z);
	    x = floatize(state);state = lcg(state);
	    y = floatize(state);state = lcg(state);
	    z = floatize(state);state = lcg(state);
	    vec3 cc(x, y, z);
	    //cross products are orthogonal and have the appropriate signed exchange
	    vec3 bbxcc = bb.cross(cc);
	    CHECK( bbxcc.dot(bb) == doctest::Approx(0) );
	    CHECK( bbxcc.dot(cc) == doctest::Approx(0) );
	    CHECK( bbxcc.dot(cc.cross(bb)) == doctest::Approx(-bbxcc.normsq()) );
	    //dot product commutivity
	    double aa_cc = aa.dot(cc);
	    double aa_bb = aa.dot(bb);
	    CHECK(aa_cc == cc.dot(aa));
	    CHECK(aa_bb == bb.dot(aa));
	    //BAC-CAB
	    vec3 aaxbbxcc = aa.cross(bbxcc);
	    vec3 baccab = bb*aa_cc - cc*aa_bb;
	    CHECK( aaxbbxcc.el[0] == doctest::Approx(baccab.el[0]) );
	    CHECK( aaxbbxcc.el[1] == doctest::Approx(baccab.el[1]) );
	    CHECK( aaxbbxcc.el[2] == doctest::Approx(baccab.el[2]) );
	}
    }
    SUBCASE("rotations") {
	//check a rotation about the z axis
	matrix<3,3> rotpiby6 = make_rotation(M_PI/6, z_axis);
	CHECK(rotpiby6.el[0] == doctest::Approx(sqrt(3)/2));	CHECK(rotpiby6.el[1] == doctest::Approx(-0.5));		CHECK(rotpiby6.el[2] == 0);
	CHECK(rotpiby6.el[3] == doctest::Approx(0.5));		CHECK(rotpiby6.el[4] == doctest::Approx(sqrt(3)/2));	CHECK(rotpiby6.el[5] == 0);
	CHECK(rotpiby6.el[6] == 0);				CHECK(rotpiby6.el[7] == 0);				CHECK(rotpiby6.el[8] == 1);
	vec3 res_mat_pow3 = rotpiby6*rotpiby6*rotpiby6*x_axis;
	CHECK(res_mat_pow3.el[0] == doctest::Approx(0));
	CHECK(res_mat_pow3.el[1] == doctest::Approx(1));
	CHECK(res_mat_pow3.el[2] == doctest::Approx(0));
	//check a rotation about the x axis
	rotpiby6 = make_rotation(M_PI/6, x_axis);
	CHECK(rotpiby6.el[0] == 1);	CHECK(rotpiby6.el[1] == doctest::Approx(0));		CHECK(rotpiby6.el[2] == 0);
	CHECK(rotpiby6.el[3] == 0);	CHECK(rotpiby6.el[4] == doctest::Approx(sqrt(3)/2));	CHECK(rotpiby6.el[5] == doctest::Approx(-0.5));
	CHECK(rotpiby6.el[6] == 0);	CHECK(rotpiby6.el[7] == doctest::Approx(0.5));		CHECK(rotpiby6.el[8] == doctest::Approx(sqrt(3)/2));
	res_mat_pow3 = rotpiby6*rotpiby6*rotpiby6*y_axis;
	CHECK(res_mat_pow3.el[0] == doctest::Approx(0));
	CHECK(res_mat_pow3.el[1] == doctest::Approx(0));
	CHECK(res_mat_pow3.el[2] == doctest::Approx(1));
	//check a rotation about the y axis
	rotpiby6 = make_rotation(M_PI/6, y_axis);
	CHECK(rotpiby6.el[0] == doctest::Approx(sqrt(3)/2));	CHECK(rotpiby6.el[1] == 0);	CHECK(rotpiby6.el[2] == doctest::Approx(0.5));
	CHECK(rotpiby6.el[3] == 0);				CHECK(rotpiby6.el[4] == 1);	CHECK(rotpiby6.el[5] == 0);
	CHECK(rotpiby6.el[6] == doctest::Approx(-0.5));		CHECK(rotpiby6.el[7] == 0);	CHECK(rotpiby6.el[8] == doctest::Approx(sqrt(3)/2));
	res_mat_pow3 = rotpiby6*rotpiby6*rotpiby6*z_axis;
	CHECK(res_mat_pow3.el[0] == doctest::Approx(1));
	CHECK(res_mat_pow3.el[1] == doctest::Approx(0));
	CHECK(res_mat_pow3.el[2] == doctest::Approx(0));
	double x,y,z;
	//make a pi/2 rotation matrix about a random axis
	x = floatize(state);state = lcg(state);
	y = floatize(state);state = lcg(state);
	z = floatize(state);state = lcg(state);
	vec3 rot_axis(x, y, z);
	rot_axis = rot_axis.normalize();
	matrix<3,3> rot = make_rotation(M_PI/2, rot_axis);
	//make random vectors and rotate them
	for (size_t i = 0; i < SMALL_TEST_N; ++i) {
	    x = floatize(state);state = lcg(state);
	    y = floatize(state);state = lcg(state);
	    z = floatize(state);state = lcg(state);
	    vec3 r(x, y, z);
	    vec3 rot_r = rot*r;
	    vec3 rot_cross = rot_r.cross(rot_axis).normalize();
	    double comp_along = r.dot(rot_axis);
	    double comp_ortho = sqrt(r.normsq() - comp_along*comp_along);
	    CHECK(r.dot(rot_cross) == doctest::Approx(comp_ortho));
	    //we rotate by 90 degrees, so only the component along the axis of rotation will be shared
	    CHECK(r.dot(rot_r) == doctest::Approx(comp_along*comp_along));
	    //now rotate again by 90 degrees, so we have the negative of the orthogonal component as well
	    rot_r = rot*rot_r;
	    CHECK(r.dot(rot_r) == doctest::Approx(2*comp_along*comp_along - r.normsq()));
	}
    }
}

TEST_CASE("line reading") {
    FILE* fp = fopen("tests/test_lines.geom", "r");
    char* buf = NULL;
    size_t bsize = 0;
    size_t lineno = 1;
    CHECK(lineno == 1);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, "This should be a line") == 0);
    CHECK(lineno == 2);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, "so should this") == 0);
    CHECK(lineno == 3);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, "    The line should end with curly {"/*}*/) == 0);
    CHECK(lineno == 4);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, "    This line should include curly{"/*}*/) == 0);
    CHECK(lineno == 7);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, /*{*/"    }") == 0);
    CHECK(lineno == 7);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, /*{*/"}") == 0);
    CHECK(lineno == 8);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, "This should be one line") == 0);
    CHECK(lineno == 9);
    CHECK(read_cgs_line(&buf, &bsize, fp, &lineno) > 0);
    CHECK(strcmp(buf, "last_test()") == 0);
    free(buf);
}

TEST_CASE("number writing") {
    char buf[BUF_SIZE];
    //simple numbers
    int res = write_number(buf, BUF_SIZE, 1, 1);
    CHECK(res == 0);
    CHECK(strcmp(buf, "1") == 0);
    res = write_number(buf, BUF_SIZE, 23, 2);
    CHECK(res == 0);
    CHECK(strcmp(buf, "23") == 0);
    res = write_number(buf, BUF_SIZE, 4, 2);
    CHECK(res == 0);
    CHECK(strcmp(buf, "04") == 0);
    res = write_number(buf, BUF_SIZE, 5, 3);
    CHECK(res == 0);
    CHECK(strcmp(buf, "005") == 0);
    //negative numbers
    res = write_number(buf, BUF_SIZE, -6, 2);
    CHECK(res == 0);
    CHECK(strcmp(buf, "-6") == 0);
    res = write_number(buf, BUF_SIZE, -78, 3);
    CHECK(res == 0);
    CHECK(strcmp(buf, "-78") == 0);

    //check invalid conditions
    res = write_number(buf, BUF_SIZE, -1, 1);
    CHECK(res < 0);
    CHECK(strlen(buf) == 1);
    res = write_number(buf, BUF_SIZE, 23, 1);
    CHECK(res < 0);
    CHECK(strlen(buf) == 1);
    res = write_number(buf, BUF_SIZE, 456, 2);
    CHECK(res < 0);
    CHECK(strlen(buf) == 2);
}

TEST_CASE("value makers") {
    //try making a list
    value val_bufs[SMALL_TEST_N];
    std::string tmp_str("foo");
    std::vector<double> tmp_darr(SMALL_TEST_N-VAL_TEST_S);
    for (size_t i = VAL_TEST_S; i < SMALL_TEST_N; ++i) {
	val_bufs[i] = make_val_num((double)i/2);
	tmp_darr[i-VAL_TEST_S] = (double)i/2;
    }
    val_bufs[0] = make_val_str(tmp_str.c_str());
    val_bufs[1] = make_val_std_str(tmp_str);
    val_bufs[2].type = VAL_LIST;
    val_bufs[2].val.l = NULL;
    val_bufs[2].n_els = 0;
    val_bufs[3] = make_val_mat(mat3x3::id());
    val_bufs[4] = make_val_vec3(vec3(0,1,2));
    val_bufs[5] = make_val_array(tmp_darr);
    value lst_val = make_val_list(val_bufs, SMALL_TEST_N);
    //make sure no element is undefined
    value tmp_u = make_val_undef();
    for (size_t i = 0; i < lst_val.n_els; ++i) {
	CHECK(lst_val.val.l[i] != tmp_u);
    }
    //make sure the other values are sane
    CHECK(lst_val.n_els == SMALL_TEST_N);
    CHECK(lst_val.val.l[0] == lst_val.val.l[1]);
    CHECK(lst_val.val.l[2].type == VAL_LIST);
    CHECK(lst_val.val.l[2].val.l == NULL);
    CHECK(lst_val.val.l[2].n_els == 0);
    CHECK(lst_val.val.l[3].type == VAL_MAT);
    CHECK(lst_val.val.l[4].type == VAL_3VEC);
    CHECK(lst_val.val.l[5].type == VAL_ARRAY);
    CHECK(lst_val.val.l[5].n_els == SMALL_TEST_N-VAL_TEST_S);
    for (size_t i = VAL_TEST_S; i < lst_val.n_els; ++i) {
	CHECK(lst_val.val.l[i].type == VAL_NUM);
	CHECK(lst_val.val.l[i].val.x == (double)i/2);
	CHECK(lst_val.val.l[5].val.a[i-VAL_TEST_S] == (double)i/2);
    }
    for (size_t i = 0; i < SMALL_TEST_N; ++i) cleanup_val(val_bufs+i);
    cleanup_val(&lst_val);
}

TEST_CASE("string castings") {
    char buf[BUF_SIZE];
    //test strings
    value v = make_val_std_str("test val");
    int write_n = v.rep_string(buf, BUF_SIZE);
    CHECK((size_t)write_n == strlen("test val"));
    CHECK(strcmp(buf, "test val") == 0);
    cleanup_val(&v);
    //test numbers
    v = make_val_num(1);
    write_n = v.rep_string(buf, BUF_SIZE);
    CHECK(write_n == 1);
    CHECK(strcmp(buf, "1") == 0);
    //test vectors
    v = make_val_vec3(vec3(0.0, 1.0, 20.0));
    write_n = v.rep_string(buf, BUF_SIZE);
    CHECK(write_n == strlen("[0,1,20]"));
    CHECK(strcmp(buf, "[0,1,20]") == 0);
    write_n = v.rep_string(buf, REP_TEST_N);
    CHECK(write_n < REP_TEST_N);
    cleanup_val(&v);
    //test matrices
    mat3x3 tmp_mat;
    v = make_val_mat(tmp_mat);
    write_n = v.rep_string(buf, BUF_SIZE);
    CHECK(write_n == strlen("[[0,0,0],[0,0,0],[0,0,0]]"));
    CHECK(strcmp(buf, "[[0,0,0],[0,0,0],[0,0,0]]") == 0);
    write_n = v.rep_string(buf, REP_TEST_N);
    CHECK(write_n < REP_TEST_N);
    cleanup_val(&v);
    //test arrays
    std::vector<double> sd(2);sd[0] = 0;sd[1] = 1;
    v = make_val_array(sd);
    write_n = v.rep_string(buf, BUF_SIZE);
    CHECK(write_n == strlen("{0,1}"));
    CHECK(strcmp(buf, "{0,1}") == 0);
    write_n = v.rep_string(buf, REP_TEST_N);
    CHECK(write_n < REP_TEST_N);
    cleanup_val(&v);
}

TEST_CASE("function parsing") {
    char buf[BUF_SIZE];

    const char* test_func_1 = "f()";
    const char* test_func_2 = "f(\"a\", \"b\", \"c\", 4)";
    const char* test_func_3 = "foo(vec(1,2,3),\"a\",\"banana\")";
    const char* test_func_4 = "foo(1, \"box(0,1,2,3)\", 4+5)";
    const char* test_func_5 = "foo ( 1 , \"b , c\" )";
    const char* test_func_6 = "f(eps = 3.5)";
    const char* test_func_7 = "f(name = \"bar\")";
    const char* test_func_8 = "f(a, b1, c2)";
    const char* bad_test_func_1 = "foo ( a , b , c";
    const char* bad_test_func_2 = "foo ( \"a\" , \"b\" , \"c\"";

    context sc;
    parse_ercode er;
    //check string 1
    strncpy(buf, test_func_1, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cgs_func cur_func = sc.parse_func(buf, 1, er, NULL);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 0);
    CHECK(strcmp(cur_func.name, "f") == 0);
    cleanup_func(&cur_func);
    //check string 2
    strncpy(buf, test_func_2, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 1, er, NULL);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 4);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "f") == 0);
    INFO("func arg=", cur_func.args[0]);
    CHECK(strcmp(cur_func.args[0].to_c_str(), "a") == 0);
    INFO("func arg=", cur_func.args[1]);
    CHECK(strcmp(cur_func.args[1].to_c_str(), "b") == 0);
    INFO("func arg=", cur_func.args[2]);
    CHECK(strcmp(cur_func.args[2].to_c_str(), "c") == 0);
    CHECK(cur_func.args[3].to_float() == 4);
    cleanup_func(&cur_func);
    //check string 3
    strncpy(buf, test_func_3, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 3, er, NULL);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 3);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "foo") == 0);
    INFO("func arg=", cur_func.args[0]);
    CHECK(cur_func.args[0].get_type() == VAL_3VEC);
    vec3* tmp_vec = cur_func.args[0].get_val().v;
    CHECK(tmp_vec->x() == 1.0);
    CHECK(tmp_vec->y() == 2.0);
    CHECK(tmp_vec->z() == 3.0);
    INFO("func arg=", cur_func.args[1].to_c_str());
    CHECK(strcmp(cur_func.args[1].to_c_str(), "a") == 0);
    INFO("func arg=", cur_func.args[1]);
    CHECK(strcmp(cur_func.args[2].to_c_str(), "banana") == 0);
    cleanup_func(&cur_func);
    //check string 4
    strncpy(buf, test_func_4, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 3, er, NULL);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 3);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "foo") == 0);
    INFO("func arg=", cur_func.args[0].to_float());
    CHECK(cur_func.args[0].to_float() == 1);
    INFO("func arg=", cur_func.args[1].to_c_str());
    CHECK(strcmp(cur_func.args[1].to_c_str(), "box(0,1,2,3)") == 0);
    INFO("func arg=", cur_func.args[2].to_float());
    CHECK(cur_func.args[2].to_float() == 9);
    cleanup_func(&cur_func);
    //check string 5
    strncpy(buf, test_func_5, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 4, er, NULL);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 2);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "foo") == 0);
    INFO("func arg=", cur_func.args[0].to_float());
    CHECK(cur_func.args[0].to_float() == 1);
    INFO("func arg=", cur_func.args[1].to_c_str());
    CHECK(strcmp(cur_func.args[1].to_c_str(), "b , c") == 0);
    cleanup_func(&cur_func);
    //check string 6
    strncpy(buf, test_func_6, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 1, er, NULL);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 1);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "f") == 0);
    INFO("func arg=", cur_func.args[0].to_float());
    CHECK(cur_func.args[0].to_float() == 3.5);
    CHECK(strcmp(cur_func.arg_names[0], "eps") == 0);
    cleanup_func(&cur_func);
    //check string 7
    strncpy(buf, test_func_7, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 1, er, NULL);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 1);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "f") == 0);
    INFO("func arg=", cur_func.args[0].to_c_str());
    CHECK(strcmp(cur_func.args[0].to_c_str(), "bar") == 0);
    CHECK(strcmp(cur_func.arg_names[0], "name") == 0);
    cleanup_func(&cur_func);
    //check string 8 (function declaration parsing)
    strncpy(buf, test_func_8, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 1, er, NULL, true);
    CHECK(er == E_SUCCESS);
    CHECK(cur_func.n_args == 3);
    INFO("func name=", cur_func.name);
    CHECK(strcmp(cur_func.name, "f") == 0);
    CHECK(strcmp(cur_func.arg_names[0], "a") == 0);
    CHECK(strcmp(cur_func.arg_names[1], "b1") == 0);
    CHECK(strcmp(cur_func.arg_names[2], "c2") == 0);
    cleanup_func(&cur_func);
    //check bad string 1
    strncpy(buf, bad_test_func_1, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 4, er, NULL);
    CHECK(er == E_BAD_SYNTAX);
    cleanup_func(&cur_func);
    //check bad string 2
    strncpy(buf, bad_test_func_2, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.parse_func(buf, 4, er, NULL);
    CHECK(er == E_BAD_SYNTAX);
    cleanup_func(&cur_func);
}

TEST_CASE("operations") {
    parse_ercode er;
    context sc;
    char buf[BUF_SIZE];
    SUBCASE("Arithmetic works") {
        //single operations
        strncpy(buf, "1+1.1", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        value tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_NUM);
        CHECK(tmp_val.val.x == 2.1);
        strncpy(buf, "2-1.25", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_NUM);
        CHECK(tmp_val.val.x == 0.75);
        strncpy(buf, "2*1.1", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_NUM);
        CHECK(tmp_val.val.x == 2.2);
        strncpy(buf, "2.2/2", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_NUM);
        CHECK(tmp_val.val.x == 1.1);
        //order of operations
        strncpy(buf, "1+3/2", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_NUM);
        CHECK(tmp_val.val.x == 2.5);
        strncpy(buf, "(1+3)/2", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_NUM);
        CHECK(tmp_val.val.x == 2.0);
        strncpy(buf, "2*9/4*3", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_NUM);
        CHECK(tmp_val.val.x == 1.5);
    }
    SUBCASE("Comparisons work") {
	//create a single true and false, this makes things easier
	value false_res = make_val_num(0);
	value true_res = make_val_num(1);
        strncpy(buf, "2 == 2", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	value tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val == true_res);
        strncpy(buf, "1 == 2", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val == false_res);
        strncpy(buf, "[2, 3] == [2, 3]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val == true_res);
        strncpy(buf, "[2, 3, 4] == [2, 3]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val == false_res);
        strncpy(buf, "[2, 3, 4] == [2, 3, 5]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val == false_res);
        strncpy(buf, "\"apple\" == \"apple\"", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val == true_res);
        strncpy(buf, "\"apple\" == \"banana\"", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val == false_res);
    }
    SUBCASE("String concatenation works") {
        //single operations
        strncpy(buf, "\"foo\"+\"bar\"", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        value tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_STR);
        CHECK(strcmp(tmp_val.val.s, "foobar") == 0);
        cleanup_val(&tmp_val);
    }
    SUBCASE("Ternary operators work") {
	strncpy(buf, "(1 == 2) ? 100 : 200", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	value tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 200);
	strncpy(buf, "(2 == 2) ? 100 : 200", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 100);
	strncpy(buf, "(1 > 2) ? 100 : 200", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 200);
	strncpy(buf, "(1 < 2) ? 100 : 200", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 100);
    }
}

TEST_CASE("builtin functions") {
    char buf[BUF_SIZE];
    value tmp_val;
    context sc;
    parse_ercode er = E_SUCCESS;

    SUBCASE("range()") {
	strncpy(buf, "range(4)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.n_els == 4);
	for (size_t i = 0; i < 4; ++i) {
	    CHECK(tmp_val.val.l[i].type == VAL_NUM);
	    CHECK(tmp_val.val.l[i].val.x == i);
	}
	cleanup_val(&tmp_val);
	strncpy(buf, "range(1,4)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.n_els == 3);
	for (size_t i = 0; i < 3; ++i) {
	    CHECK(tmp_val.val.l[i].type == VAL_NUM);
	    CHECK(tmp_val.val.l[i].val.x == i+1);
	}
	cleanup_val(&tmp_val);
	strncpy(buf, "range(1,4,0.5)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.n_els == 6);
	for (size_t i = 0; i < 6; ++i) {
	    CHECK(tmp_val.val.l[i].type == VAL_NUM);
	    CHECK(tmp_val.val.l[i].val.x == 0.5*i+1);
	}
	cleanup_val(&tmp_val);
        //graceful failure cases
        strncpy(buf, "range()", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_LACK_TOKENS);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
        strncpy(buf, "range(\"1\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_BAD_TYPE);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
        strncpy(buf, "range(0.5,\"1\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_BAD_TYPE);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
        strncpy(buf, "range(0.5,1,\"2\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_BAD_TYPE);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
    }
    SUBCASE("linspace()") {
        strncpy(buf, "linspace(1,2,5)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_ARRAY);
        CHECK(tmp_val.n_els == 5);
        for (size_t i = 0; i < 4; ++i) {
            CHECK(tmp_val.val.a[i] == 1.0+0.25*i);
        }
        cleanup_val(&tmp_val);
        strncpy(buf, "linspace(2,1,5)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_SUCCESS);
        CHECK(tmp_val.type == VAL_ARRAY);
        CHECK(tmp_val.n_els == 5);
        for (size_t i = 0; i < 4; ++i) {
            CHECK(tmp_val.val.a[i] == 2.0-0.25*i);
        }
        cleanup_val(&tmp_val);
        //graceful failure cases
        strncpy(buf, "linspace(2,1)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_LACK_TOKENS);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
        strncpy(buf, "linspace(2,1,1)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_BAD_VALUE);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
        strncpy(buf, "linspace(\"2\",1,1)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_BAD_TYPE);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
        strncpy(buf, "linspace(2,\"1\",1)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_BAD_TYPE);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
        strncpy(buf, "linspace(2,1,\"1\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
        tmp_val = sc.parse_value(buf, er);
        CHECK(er == E_BAD_TYPE);
        CHECK(tmp_val.type == VAL_UNDEF);
        CHECK(tmp_val.n_els == 0);
    }
    SUBCASE("flatten()") {
	strncpy(buf, "flatten([])", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.n_els == 0);
	cleanup_val(&tmp_val);
	strncpy(buf, "flatten([1,2,3])", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.n_els == 3);
	for (size_t i = 0; i < 3; ++i) {
	    CHECK(tmp_val.val.l[i].type == VAL_NUM);
	    CHECK(tmp_val.val.l[i].val.x == i+1);
	}
	cleanup_val(&tmp_val);
	strncpy(buf, "flatten([[1,2,3],[4,5],6])", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.n_els == 6);
	for (size_t i = 0; i < 6; ++i) {
	    CHECK(tmp_val.val.l[i].type == VAL_NUM);
	    CHECK(tmp_val.val.l[i].val.x == i+1);
	}
	cleanup_val(&tmp_val);
    }
    SUBCASE("math functions") {
	strncpy(buf, "sin(3.1415926535/2)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(1.0));
	strncpy(buf, "sin(3.1415926535/6)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(0.5));
	strncpy(buf, "cos(3.1415926535/2)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(0.0));
	strncpy(buf, "cos(3.1415926535/6)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(er == E_SUCCESS);
	strncpy(buf, "sqrt(3)/2", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	value sqrt_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(sqrt_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(sqrt_val.val.x));
	strncpy(buf, "tan(3.141592653589793/4)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(1.0));
	strncpy(buf, "tan(0)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(0.0));
	strncpy(buf, "exp(0)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(1.0));
	strncpy(buf, "exp(1)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == doctest::Approx(2.718281828));
	//failure conditions
	strncpy(buf, "sin()", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_LACK_TOKENS);
	strncpy(buf, "cos()", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_LACK_TOKENS);
	strncpy(buf, "tan()", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_LACK_TOKENS);
	strncpy(buf, "exp()", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_LACK_TOKENS);
	strncpy(buf, "sqrt()", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_LACK_TOKENS);
	strncpy(buf, "sin(\"a\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_BAD_TYPE);
	strncpy(buf, "cos(\"a\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_BAD_TYPE);
	strncpy(buf, "tan(\"a\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_BAD_TYPE);
	strncpy(buf, "exp(\"a\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_BAD_TYPE);
	strncpy(buf, "sqrt(\"a\")", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_BAD_TYPE);
    }
}

TEST_CASE("value parsing") {
    char buf[BUF_SIZE];
    context sc;
    parse_ercode er = E_SUCCESS;
    value tmp_val;

    SUBCASE("Reading numbers to values works") {
	//test integers
	strncpy(buf, "1", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 1);
	cleanup_val(&tmp_val);
	strncpy(buf, "12", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 12);
	cleanup_val(&tmp_val);
	//test floats
	strncpy(buf, ".25", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 0.25);
	cleanup_val(&tmp_val);
	strncpy(buf, "1.25", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 1.25);
	cleanup_val(&tmp_val);
	//test scientific notation
	strncpy(buf, ".25e10", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 0.25e10);
	cleanup_val(&tmp_val);
	strncpy(buf, "1.25e10", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 1.25e10);
	cleanup_val(&tmp_val);
	strncpy(buf, "1.25e-10", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_NUM);
	CHECK(tmp_val.val.x == 1.25e-10);
	cleanup_val(&tmp_val);
    }
    SUBCASE("Reading strings to values works") {
	//test a simple string
	strncpy(buf, "\"foo\"", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_STR);
	CHECK(strcmp(tmp_val.val.s, "foo") == 0);
	cleanup_val(&tmp_val);
	//test a string with whitespace
	strncpy(buf, "\" foo bar \"", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_STR);
	CHECK(strcmp(tmp_val.val.s, " foo bar ") == 0);
	cleanup_val(&tmp_val);
	//test a string with stuff inside it
	strncpy(buf, "\"foo(bar)\"", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_STR);
	CHECK(strcmp(tmp_val.val.s, "foo(bar)") == 0);
	cleanup_val(&tmp_val);
	//test a string with an escaped string
	strncpy(buf, "\"foo\\\"bar\\\" \"", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_STR);
	CHECK(strcmp(tmp_val.val.s, "foo\\\"bar\\\" ") == 0);
	cleanup_val(&tmp_val);
    }
    SUBCASE("Reading lists to values works") {
	//test one element lists
	strncpy(buf, "[\"foo\"]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.val.l != NULL);
	CHECK(tmp_val.n_els == 1);
	CHECK(tmp_val.val.l[0].type == VAL_STR);
	CHECK(strcmp(tmp_val.val.l[0].val.s, "foo") == 0);
	cleanup_val(&tmp_val);
	//test two element lists
	strncpy(buf, "[\"foo\", 1]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.val.l != NULL);
	CHECK(tmp_val.n_els == 2);
	CHECK(tmp_val.val.l[0].type == VAL_STR);
	CHECK(strcmp(tmp_val.val.l[0].val.s, "foo") == 0);
	CHECK(tmp_val.val.l[1].type == VAL_NUM);
	CHECK(tmp_val.val.l[1].val.x == 1);
	cleanup_val(&tmp_val);
	//test lists of lists
	strncpy(buf, "[[1,2,3], [\"4\", \"5\", \"6\"]]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.val.l != NULL);
	CHECK(tmp_val.n_els == 2);
	{
	    //check the first sublist
	    value element = tmp_val.val.l[0];
	    CHECK(element.type == VAL_LIST);
	    CHECK(element.n_els == 3);
	    CHECK(element.val.l != NULL);
	    CHECK(element.val.l[0].type == VAL_NUM);
	    CHECK(element.val.l[0].val.x == 1);
	    CHECK(element.val.l[1].type == VAL_NUM);
	    CHECK(element.val.l[1].val.x == 2);
	    CHECK(element.val.l[2].type == VAL_NUM);
	    CHECK(element.val.l[2].val.x == 3);
	    //check the second sublist
	    element = tmp_val.val.l[1];
	    CHECK(element.type == VAL_LIST);
	    CHECK(element.n_els == 3);
	    CHECK(element.val.l != NULL);
	    CHECK(element.val.l[0].type == VAL_STR);
	    CHECK(strcmp(element.val.l[0].val.s, "4") == 0);
	    CHECK(element.val.l[1].type == VAL_STR);
	    CHECK(strcmp(element.val.l[1].val.s, "5") == 0);
	    CHECK(element.val.l[2].type == VAL_STR);
	    CHECK(strcmp(element.val.l[2].val.s, "6") == 0);
	}
	cleanup_val(&tmp_val);
	//test lists interpretations
	strncpy(buf, "[[i*2 for i in range(2)], [i*2-1 for i in range(1,3)], [x for x in range(1,3,0.5)]]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.val.l != NULL);
	CHECK(tmp_val.n_els == 3);
	{
	    //check the first sublist
	    value element = tmp_val.val.l[0];
	    CHECK(element.type == VAL_LIST);
	    CHECK(element.n_els == 2);
	    CHECK(element.val.l != NULL);
	    CHECK(element.val.l[0].type == VAL_NUM);
	    CHECK(element.val.l[0].val.x == 0);
	    CHECK(element.val.l[1].type == VAL_NUM);
	    CHECK(element.val.l[1].val.x == 2);
	    //check the second sublist
	    element = tmp_val.val.l[1];
	    CHECK(element.type == VAL_LIST);
	    CHECK(element.n_els == 2);
	    CHECK(element.val.l != NULL);
	    CHECK(element.val.l[0].type == VAL_NUM);
	    CHECK(element.val.l[0].val.x == 1);
	    CHECK(element.val.l[1].type == VAL_NUM);
	    CHECK(element.val.l[1].val.x == 3);
	    //check the third sublist
	    element = tmp_val.val.l[2];
	    CHECK(element.type == VAL_LIST);
	    CHECK(element.n_els == 4);
	    CHECK(element.val.l != NULL);
	    CHECK(element.val.l[0].type == VAL_NUM);
	    CHECK(element.val.l[0].val.x == 1);
	    CHECK(element.val.l[1].type == VAL_NUM);
	    CHECK(element.val.l[1].val.x == 1.5);
	    CHECK(element.val.l[2].type == VAL_NUM);
	    CHECK(element.val.l[2].val.x == 2);
	    CHECK(element.val.l[3].type == VAL_NUM);
	    CHECK(element.val.l[3].val.x == 2.5);
	}
	cleanup_val(&tmp_val);
	//test nested list interpretations
	strncpy(buf, "[[x*y for x in range(1,6)] for y in range(5)]", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_LIST);
	CHECK(tmp_val.val.l != NULL);
	CHECK(tmp_val.n_els == 5);
	for (size_t yy = 0; yy < tmp_val.n_els; ++yy) {
	    CHECK(tmp_val.val.l[yy].type == VAL_LIST);
	    CHECK(tmp_val.val.l[yy].n_els == 5);
	    for (size_t xx = 0; xx < tmp_val.val.l[yy].n_els; ++xx) {
		CHECK(tmp_val.val.l[yy].val.l[xx].type == VAL_NUM);
		CHECK(tmp_val.val.l[yy].val.l[xx].val.x == (xx+1)*yy);
	    }
	}
	cleanup_val(&tmp_val);
    }
    SUBCASE("Reading vectors to values works") {
	//test one element lists
	strncpy(buf, "vec(1.2, 3.4,56.7)", BUF_SIZE);buf[BUF_SIZE-1] = 0;
	tmp_val = sc.parse_value(buf, er);
	CHECK(er == E_SUCCESS);
	CHECK(tmp_val.type == VAL_3VEC);
	CHECK(tmp_val.val.v != NULL);
	CHECK(tmp_val.val.v->x() == doctest::Approx(1.2));
	CHECK(tmp_val.val.v->y() == doctest::Approx(3.4));
	CHECK(tmp_val.val.v->z() == doctest::Approx(56.7));
	cleanup_val(&tmp_val);
    }
}

TEST_CASE("line_buffer splitting") {
    const char* lines[] = { "apple; banana;c", ";orange" };
    size_t n_lines = sizeof(lines)/sizeof(char*);
    line_buffer lb(lines, n_lines);
    CHECK(lb.get_n_lines() == 2);
    char* strval = lb.get_line(0);CHECK(strcmp(lines[0], strval) == 0);free(strval);
    strval = lb.get_line(1);CHECK(strcmp(lines[1], strval) == 0);free(strval);
    lb.split(';');
    CHECK(lb.get_n_lines() == 5);
    strval = lb.get_line(0);CHECK(strcmp("apple", strval) == 0);free(strval);
    strval = lb.get_line(1);CHECK(strcmp(" banana", strval) == 0);free(strval);
    strval = lb.get_line(2);CHECK(strcmp("c", strval) == 0);free(strval);
    strval = lb.get_line(3);CHECK(strcmp("", strval) == 0);free(strval);
    strval = lb.get_line(4);CHECK(strcmp("orange", strval) == 0);free(strval);
}

TEST_CASE("test line_buffer it_single") {
    const char* lines[] = { "def test_fun(a)", "{", "if a > 5 {", "return 1", "}", "return 0", "}" };
    size_t n_lines = sizeof(lines)/sizeof(char*);
    //check the lines (curly brace on new line)
    line_buffer lb(lines, n_lines);
    //test it_single
    int depth = 0;
    char* it_single_str = NULL;
    //including braces
    for (size_t i = 0; i < n_lines; ++i) {
	line_buffer_ind start(i, 0);
	line_buffer_ind end(i, 0);
	int it_start = lb.it_single(&it_single_str, '{', '}', &start, &end, &depth, true, false);
	if (i == 1)
	    CHECK(it_start == 0);
	else
	    CHECK(it_start == -1);
	CHECK(strcmp(it_single_str, lines[i]) == 0);
	CHECK(start.line == i);
	CHECK(start.off == 0);
	if (i < n_lines-1)
	    CHECK(end.line == start.line+1);
	else
	    CHECK(end.line == start.line);
	CHECK(end.off == strlen(lines[i]));
	if (i == 0 or i == n_lines-1)
	    CHECK(depth == 0);
	else
	    CHECK(depth >= 1);
    }
    //excluding braces
    const char* lines_exc[] = { "def test_fun(a)", "", "if a > 5 {", "return 1", "}", "return 0", "}" };
    for (size_t i = 0; i < n_lines; ++i) {
	line_buffer_ind start(i, 0);
	line_buffer_ind end(i, 0);
	int it_start = lb.it_single(&it_single_str, '{', '}', &start, &end, &depth, false, false);
	if (i == 1)
	    CHECK(it_start == 0);
	else
	    CHECK(it_start == -1);
	CHECK(strcmp(it_single_str, lines_exc[i]) == 0);
	//check that lines start where expected
	CHECK(start.line == i);
	if (i == 1)
	    CHECK(start.off == 1);
	else
	    CHECK(start.off == 0);
	if (i < n_lines-1)
	    CHECK(end.line == start.line+1);
	else
	    CHECK(end.line == start.line);
	//check that lines end where expected
	if (i == n_lines-1)
	    CHECK(end.off == 0);
	else
	    CHECK(end.off == strlen(lines[i]));
	//CHECK that the depths are correct
	if (i == 0 or i == n_lines-1)
	    CHECK(depth == 0);
	else
	    CHECK(depth >= 1);
    }
}

TEST_CASE("line_buffer get_enclosed") {
    const char* fun_contents[] = {"", "if a > 5 {", "return 1", "}", "return 0", ""};
    const char* if_contents[] = {"", "return 1", ""};
    size_t fun_n = sizeof(fun_contents)/sizeof(char*);
    size_t if_n = sizeof(if_contents)/sizeof(char*);

    line_buffer_ind end_ind;
    char* strval = NULL;

    SUBCASE("open brace on a different line") {
	const char* lines[] = { "def test_fun(a)", "{", "if a > 5 {", "return 1", "}", "return 0", "}" };
	size_t n_lines = sizeof(lines)/sizeof(char*);
	//check the lines (curly brace on new line)
	line_buffer lb(lines, n_lines);
	for (size_t i = 0; i < n_lines; ++i) {
	    strval = lb.get_line(i);CHECK(strcmp(lines[i], strval) == 0);free(strval);
	}
	//test wrapper functions that use it_single
	line_buffer_ind bstart;
	line_buffer b_fun_con = lb.get_enclosed(bstart, &end_ind, '{', '}');
	CHECK(b_fun_con.get_n_lines() == fun_n);
	CHECK(end_ind.line == 6);
	for (size_t i = 0; i < fun_n; ++i) {
	    strval = b_fun_con.get_line(i);CHECK(strcmp(fun_contents[i], strval) == 0);free(strval);
	}
	line_buffer b_if_con = b_fun_con.get_enclosed(bstart, &end_ind, '{', '}');
	CHECK(b_if_con.get_n_lines() == if_n);
	CHECK(end_ind.line == 3);
	for (size_t i = 0; i < if_n; ++i) {
	    strval = b_if_con.get_line(i);CHECK(strcmp(if_contents[i], strval) == 0);free(strval);
	}
	//check jumping
	line_buffer_ind blk_start(0,0);
	line_buffer_ind blk_end_ind = lb.jmp_enclosed(blk_start, '{', '}');
	CHECK(blk_end_ind.line == 6);
	CHECK(blk_end_ind.off == 0);
	blk_end_ind = lb.jmp_enclosed(blk_start, '{', '}', true);
	CHECK(blk_end_ind.line == 6);
	CHECK(blk_end_ind.off == 1);
	//try flattening
	char* fun_flat = b_fun_con.flatten();
	CHECK(strcmp(fun_flat, "if a > 5 {return 1}return 0") == 0);
	free(fun_flat);
	fun_flat = b_fun_con.flatten('|');
	CHECK(strcmp(fun_flat, "if a > 5 {|return 1|}|return 0||") == 0);
	free(fun_flat);
    }

    SUBCASE("open brace on the same line") {
	const char* lines[] = { "def test_fun(a) {", "if a > 5 {", "return 1", "}", "return 0", "}" };
	size_t n_lines = sizeof(lines)/sizeof(char*);
	//check the lines (curly brace on same line)
	line_buffer lb(lines, n_lines);
	for (size_t i = 0; i < n_lines; ++i) {
	    strval = lb.get_line(i);CHECK(strcmp(lines[i], strval) == 0);free(strval);
	}
	//test wrapper functions that use it_single
	line_buffer_ind bstart;
	line_buffer b_fun_con = lb.get_enclosed(bstart, &end_ind, '{', '}');
	CHECK(b_fun_con.get_n_lines() == fun_n);
	CHECK(end_ind.line == 5);
	for (size_t i = 0; i < fun_n; ++i) {
	    strval = b_fun_con.get_line(i);CHECK(strcmp(fun_contents[i], strval) == 0);free(strval);
	}
	line_buffer b_if_con_2 = b_fun_con.get_enclosed(bstart, &end_ind, '{', '}');
	CHECK(b_if_con_2.get_n_lines() == if_n);
	CHECK(end_ind.line == 3);
	for (size_t i = 0; i < if_n; ++i) {
	    strval = b_if_con_2.get_line(i);CHECK(strcmp(if_contents[i], strval) == 0);free(strval);
	}
	//check jumping
	line_buffer_ind blk_start(0,0);
	line_buffer_ind blk_end_ind = b_fun_con.jmp_enclosed(blk_start, '{', '}');
	CHECK(blk_end_ind.line == 3);
	CHECK(blk_end_ind.off == 0);
	blk_end_ind = b_fun_con.jmp_enclosed(blk_start, '{', '}', true);
	CHECK(blk_end_ind.line == 3);
	CHECK(blk_end_ind.off == 1);
	//try flattening
	char* fun_flat = b_fun_con.flatten();
	CHECK(strcmp(fun_flat, "if a > 5 {return 1}return 0") == 0);
	free(fun_flat);
	fun_flat = b_fun_con.flatten('|');
	CHECK(strcmp(fun_flat, "if a > 5 {|return 1|}|return 0||") == 0);
	free(fun_flat);
    }

    SUBCASE("everything on one line") {
	const char* lines[] = { "def test_fun(a) {if a > 5 {return 1}return 0}" };
	size_t n_lines = sizeof(lines)/sizeof(char*);
	//check the lines (curly brace on new line)
	line_buffer lb(lines, n_lines);
	for (size_t i = 0; i < n_lines; ++i) {
	    strval = lb.get_line(i);CHECK(strcmp(lines[i], strval) == 0);free(strval);
	}
	//test wrapper functions that use it_single
	line_buffer_ind bstart;
	line_buffer b_fun_con = lb.get_enclosed(bstart, &end_ind, '{', '}');
	CHECK(b_fun_con.get_n_lines() == 1);
	CHECK(end_ind.line == 0);
	strval = b_fun_con.get_line(0);CHECK(strcmp("if a > 5 {return 1}return 0", strval) == 0);free(strval);
	line_buffer b_if_con = b_fun_con.get_enclosed(bstart, &end_ind, '{', '}');
	CHECK(b_if_con.get_n_lines() == 1);
	CHECK(end_ind.line == 0);
	strval = b_if_con.get_line(0);CHECK(strcmp("return 1", strval) == 0);free(strval);
	//check jumping
	line_buffer_ind blk_start;
	line_buffer_ind blk_end_ind = b_fun_con.jmp_enclosed(blk_start, '{', '}');
	CHECK(blk_end_ind.line == 0);
	CHECK(blk_end_ind.off == 18);
	blk_end_ind = b_fun_con.jmp_enclosed(blk_start, '{', '}', true);
	CHECK(blk_end_ind.line == 0);
	CHECK(blk_end_ind.off == 19);
	//try flattening
	char* fun_flat = b_fun_con.flatten();
	CHECK(strcmp(fun_flat, "if a > 5 {return 1}return 0") == 0);
	free(fun_flat);
	fun_flat = b_fun_con.flatten('|');
	CHECK(strcmp(fun_flat, "if a > 5 {return 1}return 0|") == 0);
	free(fun_flat);
    }
}

value test_fun_call(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    double a = f.args[0].val.x;
    if (a > 5) {
	ret.type = VAL_INST;
	ret.val.c = new context();
	value name_contents = make_val_str("hi");
	ret.val.c->emplace("name", name_contents);
	cleanup_val(&name_contents);
	return ret;
    }
    return f.args[0];
}

value test_fun_gamma(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    double a = f.args[0].val.x;
    return make_val_num(sqrt(1 - a*a));
}

TEST_CASE("context parsing") {
    SUBCASE ("without nesting") {
	const char* lines[] = { "a = 1", "\"b\"", "c = [\"d\", \"e\"]" }; 
	size_t n_lines = sizeof(lines)/sizeof(char*);
	line_buffer b_1(lines, n_lines);
	context c;
	parse_ercode er = c.read_from_lines(b_1);
	CHECK(er == E_SUCCESS);
	//lookup the named values
	value val_a = c.lookup("a");
	CHECK(val_a.type == VAL_NUM);
	CHECK(val_a.val.x == 1);
	value val_c = c.lookup("c");
	CHECK(val_c.type == VAL_LIST);
	CHECK(val_c.n_els == 2);
	CHECK(val_c.val.l[0].type == VAL_STR);
	CHECK(val_c.val.l[1].type == VAL_STR);
	//lookup the stack value
	name_val_pair pair_b = c.peek(2);
	value val_b = pair_b.get_val();
	CHECK(val_b.type == VAL_STR);
	printf("val_b = %s", val_b.val.s);
	CHECK(strcmp(val_b.val.s, "b") == 0);
    }
    SUBCASE ("with nesting") {
	const char* lines[] = { "a = {name = \"apple\", values = [20, 11]}", "b = a.values[0]", "c = a.values[1] + a.values[0] + 1" }; 
	size_t n_lines = sizeof(lines)/sizeof(char*);
	line_buffer b_1(lines, n_lines);
	context c;
	parse_ercode er = c.read_from_lines(b_1);
	CHECK(er == E_SUCCESS);
	//lookup the named values
	value val_a = c.lookup("a");
	CHECK(val_a.type == VAL_INST);
	value val_a_name = val_a.val.c->lookup("name");
	CHECK(val_a_name.type == VAL_STR);
	CHECK(strcmp(val_a_name.val.s, "apple") == 0);
	value val_a_value = val_a.val.c->lookup("values");
	CHECK(val_a_value.type == VAL_LIST);
	CHECK(val_a_value.n_els == 2);
	CHECK(val_a_value.val.l[0].type == VAL_NUM);
	CHECK(val_a_value.val.l[1].type == VAL_NUM);
	CHECK(val_a_value.val.l[0].val.x == 20);
	CHECK(val_a_value.val.l[1].val.x == 11);
	value val_b = c.lookup("b");
	CHECK(val_b.type == VAL_NUM);
	CHECK(val_b.val.x == 20);
	value val_c = c.lookup("c");
	CHECK(val_c.type == VAL_NUM);
	CHECK(val_c.val.x == 32);
    }
    SUBCASE ("user defined functions") {
	const char* fun_name = "test_fun";
	char* tmp_name = strdup(fun_name);

	const char* lines[] = { "a = test_fun(1)", "b=test_fun(10)" };
	size_t n_lines = sizeof(lines)/sizeof(char*);
	line_buffer b_1(lines, n_lines);
	context c;
	value tmp_f = make_val_func("test_fun", 1, &test_fun_call);
	c.emplace("test_fun", tmp_f);
	cleanup_val(&tmp_f);
	parse_ercode er = c.read_from_lines(b_1);
	CHECK(er == E_SUCCESS);
	if (!er) {
	    //make sure that the function is there
	    value val_fun = c.lookup("test_fun");
	    CHECK(val_fun.type == VAL_FUNC);
	    //make sure that the number value a is there
	    value val_a = c.lookup("a");
	    CHECK(val_a.type == VAL_NUM);
	    CHECK(val_a.val.x == 1);
	    value val_b = c.lookup("b");
	    CHECK(val_b.type == VAL_INST);
	    value val_b_name = val_b.val.c->lookup("name");
	    CHECK(val_b_name.type == VAL_STR);
	    CHECK(strcmp(val_b_name.val.s, "hi") == 0);
	    /*cleanup_val(&val_fun);
	    cleanup_val(&val_a);
	    cleanup_val(&val_b);*/
	}
	free(tmp_name);
    }
    SUBCASE ("stress test") {
	//first we add a bunch of arbitrary variables to make searching harder for the parser
	const char* lines1[] = {
	    "Vodkis=1","Pagne=2","Meadaj=3","whis=4","nac4=5","RaKi=6","gyn=7","cid=8","Daiqui=9","Mooshi=10","Magnac=2","manChe=3","tes=4","Bourbu=5","magna=6","sak=7","Para=8","Keffi=9","Guino=10",
	    "y = 2.0", "xs = linspace(0, y, 10000)", "arr1 = [sin(6*x/y) for x in xs]" };
	size_t n_lines1 = sizeof(lines1)/sizeof(char*);
	line_buffer b_1(lines1, n_lines1);
	const char* lines2[] = { "arr2 = [gam(x/y) for x in xs]" };
	size_t n_lines2 = sizeof(lines2)/sizeof(char*);
	line_buffer b_2(lines2, n_lines2);
	context c;
	parse_ercode er = c.read_from_lines(b_1);
	CHECK(er == E_SUCCESS);
	value tmp_f = make_val_func("gam", 1, &test_fun_gamma);
	c.emplace("gam", tmp_f);
	cleanup_val(&tmp_f);
	er = c.read_from_lines(b_2);
	CHECK(er == E_SUCCESS);
    }
}

TEST_CASE("context file parsing") {
    line_buffer lb("tests/context_test.geom");
    CHECK(lb.get_n_lines() == 7);
    context c;
    setup_geometry_context(c);
    size_t init_size = c.size();
    parse_ercode er = c.read_from_lines(lb);
    CHECK(er == E_SUCCESS);
    CHECK(c.size() == init_size+4);
    //look at the Gaussian
    value inst = c.peek_val(4); {
	CHECK(inst.type == VAL_INST);
	CHECK(inst.val.c->size() == 9);
	name_val_pair strval = inst.val.c->peek(inst.val.c->size());
	CHECK(strval.name_matches("__type__"));
	CHECK(strval.get_val().type == VAL_STR);
	CHECK(strcmp(strval.get_val().val.s, "Gaussian_source") == 0);
	value tmp = inst.val.c->lookup("component");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == (double)C_EX);
	tmp = inst.val.c->lookup("wavelength");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == doctest::Approx(1.33));
	tmp = inst.val.c->lookup("amplitude");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 1.0);
	tmp = inst.val.c->lookup("width");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 2.0);
	tmp = inst.val.c->lookup("phase");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 0.0);
	tmp = inst.val.c->lookup("cutoff");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 0.125);
	tmp = inst.val.c->lookup("start_time");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == -1.25);
	CHECK(is_type(inst.val.c->peek_val(), "Box"));
    }
    inst = c.peek_val(3); {
	CHECK(inst.type == VAL_INST);
	CHECK(inst.val.c->size() == 8);
	name_val_pair strval = inst.val.c->peek(inst.val.c->size());
	CHECK(strval.name_matches("__type__"));
	CHECK(strval.get_val().type == VAL_STR);
	CHECK(strcmp(strval.get_val().val.s, "CW_source") == 0);
	value tmp = inst.val.c->lookup("component");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == (double)C_HZ);
	tmp = inst.val.c->lookup("wavelength");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 0.625);
	tmp = inst.val.c->lookup("amplitude");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 0.25);
	tmp = inst.val.c->lookup("start_time");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 0.75);
	tmp = inst.val.c->lookup("end_time");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 2.25);
	tmp = inst.val.c->lookup("slowness");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 12);
	CHECK(is_type(inst.val.c->peek_val(), "Box"));
    }
    inst = c.peek_val(2); {
	CHECK(inst.type == VAL_INST);
	CHECK(inst.val.c->size() == 5);
	name_val_pair strval = inst.val.c->peek(inst.val.c->size());
	CHECK(strval.name_matches("__type__"));
	CHECK(strval.get_val().type == VAL_STR);
	CHECK(strcmp(strval.get_val().val.s, "Composite") == 0);
	value tmp = inst.val.c->lookup("eps");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 3.5);
	tmp = inst.val.c->lookup("alpha");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 1);
	tmp = inst.val.c->lookup("color");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 10);
	value geom = inst.val.c->peek_val();
	CHECK(geom.type == VAL_LIST);
	CHECK(geom.n_els == 2);
	CHECK(is_type(geom.val.l[0], "Box"));
	CHECK(is_type(geom.val.l[1], "Box"));
    }
    inst = c.peek_val(1); {
	CHECK(inst.type == VAL_INST);
	CHECK(inst.val.c->size() == 9);
	name_val_pair strval = inst.val.c->peek(inst.val.c->size());
	CHECK(strval.name_matches("__type__"));
	CHECK(strval.get_val().type == VAL_STR);
	CHECK(strcmp(strval.get_val().val.s, "snapshot") == 0);
	value tmp = inst.val.c->lookup("fname");
	CHECK(tmp.type == VAL_STR);CHECK(strcmp(tmp.val.s, "/tmp/run_alpha.pgm") == 0);
	tmp = inst.val.c->lookup("cam_v");
	CHECK(tmp.type == VAL_3VEC);
	tmp = inst.val.c->lookup("look_v");
	CHECK(tmp.type == VAL_3VEC);
	tmp = inst.val.c->lookup("up_v");
	CHECK(tmp.type == VAL_3VEC);
	tmp = inst.val.c->lookup("scale");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == 1);
	tmp = inst.val.c->lookup("res");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == DEF_IM_RES);
	tmp = inst.val.c->lookup("n_samples");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == DEF_TEST_N);
	tmp = inst.val.c->lookup("step");
	CHECK(tmp.type == VAL_NUM);CHECK(tmp.val.x == WALK_STEP);
    }
}

TEST_CASE("volumes") {
    //initialize random state
    _uint state = lcg(lcg(TEST_SEED));
    double x, y, z;
    //test volumes by sampling from a cube with side lengths unit 1 and testing the fraction that are included
    vec3 center(0.5, 0.5, 0.5);
    vec3 corner_1(1, 1, 0);
    vec3 corner_2(0, 1, 1);
    vec3 corner_3(1, 0, 1);
    sphere test_sphere(center, 0.5);
    cylinder test_cyl(center, 0.5, 0.5, 0.5);
    plane test_plane(corner_1, corner_2, corner_3);
    double sphere_frac = 0;
    double plane_frac = 0;
    double cyl_frac = 0;
    //sample random points
    for (size_t i = 0; i < TEST_N; ++i) {
	state = lcg(state);
	x = floatize(state);
	state = lcg(state);
	y = floatize(state);
	state = lcg(state);
	z = floatize(state);
	state = lcg(state);
	vec3 r(x, y, z);
	if (test_sphere.in(r))
	    sphere_frac += 1.0/TEST_N;
	if (test_plane.in(r))
	    plane_frac += 1.0/TEST_N;
	if (test_cyl.in(r))
	    cyl_frac += 1.0/TEST_N;
    }
    double v_sphere = 1.33333*M_PI*0.125;
    double v_plane = 5.0/6;
    double v_cyl = M_PI*0.125;
    
    CHECK(abs(v_sphere - sphere_frac)/v_sphere < EPSILON);
    CHECK(abs(v_plane - plane_frac)/v_plane < EPSILON);
    CHECK(abs(v_cyl - cyl_frac)/v_cyl < EPSILON);
    //draw test images
    vec3 cam_pos(CAM_X, CAM_Y, CAM_Z);
}

TEST_CASE("object Trees") {
    //declare variables
    char buf[BUF_SIZE];
    object_stack test_stack;
    object* cur_obj;object_type cur_type;
    parse_ercode er;

    //setup a bunch of strings describing objects
    const char* root_obj_str = "Composite(eps = 3.5)";
    const char* l_str = "union()";
    const char* ll_str = "Box(vec(0,0,0), vec(1,1,1))";
    const char* lr_str = "Sphere(vec(2,0,0), 1)";
    const char* r_str = "intersect()";
    const char* rl_str = "Box(vec(0,0,0), vec(1,1,1))";
    const char* rr_str = "Cylinder(vec(2,0,0), 1, 1)";

    scene sc;
    //Insert root object
    strncpy(buf, root_obj_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cgs_func cur_func = sc.get_context().parse_func(buf, (size_t)(strchr(buf, '(')-buf), er, NULL);
    CHECK(er == E_SUCCESS);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    cleanup_func(&cur_func);
    //Insert object 1
    strncpy(buf, l_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.get_context().parse_func(buf, (size_t)(strchr(buf, '(')-buf), er, NULL);
    CHECK(er == E_SUCCESS);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    cleanup_func(&cur_func);
    //Insert left union object
    strncpy(buf, ll_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.get_context().parse_func(buf, (size_t)(strchr(buf, '(')-buf), er, NULL);
    CHECK(er == E_SUCCESS);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    cleanup_func(&cur_func);
    //Insert right union object
    strncpy(buf, lr_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.get_context().parse_func(buf, (size_t)(strchr(buf, '(')-buf), er, NULL);
    CHECK(er == E_SUCCESS);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    cleanup_func(&cur_func);
    //Insert object 2
    strncpy(buf, r_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.get_context().parse_func(buf, (size_t)(strchr(buf, '(')-buf), er, NULL);
    CHECK(er == E_SUCCESS);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    cleanup_func(&cur_func);
    //Insert left union object
    strncpy(buf, rl_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.get_context().parse_func(buf, (size_t)(strchr(buf, '(')-buf), er, NULL);
    CHECK(er == E_SUCCESS);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    cleanup_func(&cur_func);
    //Insert right union object
    strncpy(buf, rr_str, BUF_SIZE);buf[BUF_SIZE-1] = 0;
    cur_func = sc.get_context().parse_func(buf, (size_t)(strchr(buf, '(')-buf), er, NULL);
    CHECK(er == E_SUCCESS);
    er = sc.make_object(cur_func, &cur_obj, &cur_type, 0);
    CHECK(er == E_SUCCESS);
    test_stack.emplace_obj(cur_obj, cur_type);
    cleanup_func(&cur_func);

    //get all composite objects in the tree
    composite_object* root = test_stack.get_root();
    CHECK(root != NULL);
    composite_object* comp_l = (composite_object*)(root->get_child_l());
    composite_object* comp_r = (composite_object*)(root->get_child_r());
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

TEST_CASE("File Parsing") {
    parse_ercode er;
    scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    //look at context variables
    context c = s.get_context();
    value offset = c.lookup("offset");
    value list = c.lookup("list");
    value sum_list = c.lookup("sum_list");
    value prod_list = c.lookup("prod_list");
    value acid_test = c.lookup("acid_test");
    CHECK(offset.type == VAL_NUM);CHECK(offset.val.x == 0.2);
    CHECK(list.type == VAL_LIST);
    CHECK(sum_list.type == VAL_NUM);CHECK(sum_list.val.x == 10);
    CHECK(prod_list.type == VAL_NUM);CHECK(prod_list.val.x == 24.2);
    CHECK(acid_test.type == VAL_NUM);CHECK(acid_test.val.x == 16);

    //test geometric information
    std::vector<composite_object*> roots_vec = s.get_roots();
    CHECK(roots_vec.size() > 0);
    composite_object* root = roots_vec[0];

    CHECK(root != NULL);
    //check that all are not NULL and that types are correct
    CHECK(root->get_child_type_l() == CGS_COMPOSITE);
    CHECK(root->get_child_type_r() == CGS_COMPOSITE);
    composite_object* comp_l = (composite_object*)(root->get_child_l());
    composite_object* comp_r = (composite_object*)(root->get_child_r());
    CHECK(comp_r->get_child_type_l() == CGS_COMPOSITE);
    composite_object* comp_rl = (composite_object*)(comp_r->get_child_l());
    composite_object* comp_rr = (composite_object*)(comp_r->get_child_r());
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
    CHECK(comp_rl->get_combine_type() == CGS_INTERSECT);
    CHECK(comp_rl->get_child_l() != NULL);
    CHECK(comp_rl->get_child_type_l() == CGS_BOX);
    CHECK(comp_rl->get_child_r() != NULL);
    CHECK(comp_rl->get_child_type_r() == CGS_COMPOSITE);
    CHECK(comp_rl->get_child_r() != NULL);
    composite_object* comp_rlr = (composite_object*)(comp_rl->get_child_r());
    CHECK(comp_rlr->get_child_type_l() == CGS_PLANE);
    CHECK(comp_rlr->get_child_l() != NULL);
    CHECK(comp_rlr->get_child_type_r() == CGS_UNDEF);
    CHECK(comp_rlr->get_child_r() == NULL);
    CHECK(comp_r->get_child_type_r() == CGS_CYLINDER);
    CHECK(comp_r->get_child_r() != NULL);
}

TEST_CASE("geometric inclusion") {
    parse_ercode er;
    scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    composite_object* root = s.get_roots()[0];

    CHECK(root->in(vec3(.45,.45,.45)) == 1);
    CHECK(root->in(vec3(.45,.41,.6)) == 1);
    CHECK(root->in(vec3(.65,.45,.41)) == 1);
    CHECK(root->in(vec3(.71,.51,.51)) == 0);
    CHECK(root->in(vec3(.55,.45,.45)) == 0);
    CHECK(root->in(vec3(.55,.41,.85)) == 0);
    vec3 cam_pos(CAM_X, CAM_Y, CAM_Z);
    s.draw("/tmp/test_composite.pgm", cam_pos);
}

uint32_t set_alpha(uint32_t col, uint32_t a) {
    return (col & 0x00ffffff) | (a << 24);
}
TEST_CASE("image saving") {
    //check that blending works
    uint32_t c1 = make_col(255, 0, 0);
    uint32_t c2 = make_col(0, 255, 127);
    uint32_t blend_res = blend(c1, c2);
    CHECK(blend_res == c1);
    c1 = set_alpha(c1, 128);
    blend_res = blend(c1, c2);
    CHECK(get_a(blend_res) == 255);
    CHECK(get_r(blend_res) == 127);
    CHECK(get_g(blend_res) == 127);
    CHECK(get_b(blend_res) == 63);
    c2 = set_alpha(c2, 192);
    blend_res = blend(c1, c2);
    CHECK(get_a(blend_res) == 224);
    CHECK(get_r(blend_res) == 145);
    CHECK(get_g(blend_res) == 109);
    CHECK(get_b(blend_res) == 54);
    //check that 0 opacity results in only one color
    c1 = set_alpha(c1, 255);
    c2 = set_alpha(c2, 0);
    CHECK(blend(c1, c2) == c1);
    CHECK(blend(c2, c1) == c1);
    //set up an image buffer
    size_t res = 255;
    uint32_t c_buf[res*res];
    for (size_t i = 0; i < res; ++i) {
	//transparent to red to transparent gradient on the y direction
	int a_mask = ((int)i - res/2);
	a_mask = 255 - (a_mask*a_mask)/64;
	c1 = set_alpha(c1, a_mask);
	for (size_t j = 0; j < res; ++j) {
	    //green to blue gradient in the x direction
	    c2 = make_col(0, 255, j);
	    c_buf[i*res + j] = blend(c1, c2);
	}
    }
    save_imbuf("/tmp/tst_img.pgm", c_buf, res, res);
    //now load a file with test information
    parse_ercode er;
    scene s("tests/alpha.geom", &er);
}

TEST_CASE("dispersion material volumentric inclusion") {
    parse_ercode er;
    //load settings from the configuration file
    parse_settings args;
    std::string name = "tests/test.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
    free(name_dup);
    CHECK(ret == 0);

    //create the geometry object
    scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    composite_object* root = s.get_roots()[0];
    cgs_material_function mat_func(root);

    //check locations
    meep::vec test_loc_1(.45,.45,.45);
    meep::vec test_loc_2(.45,.41,.6);
    meep::vec test_loc_3(.65,.45,.45);
    meep::vec test_loc_4(.71,.51,.51);
    meep::vec test_loc_5(.55,.45,.45);
    meep::vec test_loc_6(.55,.41,.85);
    CHECK(mat_func.in_bound(test_loc_1) == 3.5);
    CHECK(mat_func.in_bound(test_loc_2) == 3.5);
    CHECK(mat_func.in_bound(test_loc_3) == 3.5);
    CHECK(mat_func.in_bound(test_loc_4) == 1.0);
    CHECK(mat_func.in_bound(test_loc_5) == 1.0);
    CHECK(mat_func.in_bound(test_loc_6) == 1.0);

    cleanup_settings(&args);
}

TEST_CASE("reading of configuration files") {
    std::string name = "tests/test.conf";
    char* name_dup = strdup(name.c_str());

    SUBCASE("Reading just a config file works") {
	parse_settings args;
	parse_conf_file(&args, name_dup);

	CHECK(args.n_dims == 3);
	CHECK(args.pml_thickness == 2.0);
	CHECK(args.len == 16.0);
	CHECK(args.um_scale == 1.0);
	//NOTE that resolutions are rounded up to accomodate the nearest odd integer number of grid points. This value should NOT be the same as that given in the .conf file.
	CHECK(args.resolution == 1.55);
	CHECK(args.courant == 0.3);
	CHECK(args.smooth_n == 1);
	CHECK(args.smooth_rad == 0.25);
	//CHECK(strcmp(args.monitor_locs, "(1.0,1.0,1.0)") == 0);
	CHECK(args.post_source_t == 1.0);
	CHECK(args.field_dump_span == 171);
	CHECK(args.ambient_eps == 1.0);
	CHECK(strcmp(args.geom_fname, "tests/test.geom") == 0);

	cleanup_settings(&args);
    }

    SUBCASE("Command line arguments override defaults") {
	parse_settings args;
	//parse commandline arguments
	std::string sim_argv_cpp[] = { "./test", "--conf-file", "blah.conf", "--geom-file", "blah.geom", "--out-dir", "/blah", "--grid-res", "3.0", "--length", "9.0", "--eps1", "2.0", "--opts", "a = 0.1; b_option = [1,\"blah\"]" };
	//gcc doesn't like using string literals as c-strings, ugh.
	size_t n_args = sizeof(sim_argv_cpp)/sizeof(std::string);
	//need to create two lists since parse_args modifies the list in place. The second list is only used so that we know which addresses to free.
	char** sim_argv_c = (char**)malloc(sizeof(char*)*n_args);
	char** post_sim_argv_c = (char**)malloc(sizeof(char*)*n_args);
	for (_uint i = 0; i < n_args; ++i) {
	    sim_argv_c[i] = strdup(sim_argv_cpp[i].c_str());
	    post_sim_argv_c[i] = sim_argv_c[i];
	}

	//finally we can parse the command line arguments
	int post_n_args = n_args;
	parse_args(&args, &post_n_args, sim_argv_c);
	CHECK(post_n_args == 1);

	//this is used when calling parse_args, so it should be checked before everything else
	CHECK(args.conf_fname != NULL);
	CHECK(strcmp(args.conf_fname, "blah.conf") == 0);

	//read the config file
	parse_conf_file(&args, name_dup);
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
	//CHECK(strcmp(args.monitor_locs, "(1.0,1.0,1.0)") == 0);
	CHECK(args.post_source_t == 1.0);

	//check that contexts are loaded with the correct options
	context con = context_from_settings(args);
	value tmp = con.lookup("pml_thickness");
	CHECK(tmp.type == VAL_NUM);
	CHECK(tmp.val.x == 2.0);
	tmp = con.lookup("length");
	CHECK(tmp.type == VAL_NUM);
	CHECK(tmp.val.x == 9.0+2*2.0);
	tmp = con.lookup("l_per_um");
	CHECK(tmp.type == VAL_NUM);
	CHECK(tmp.val.x == 1.0);
	tmp = con.lookup("a");
	CHECK(tmp.type == VAL_NUM);
	CHECK(tmp.val.x == 0.1);
	tmp = con.lookup("b_option");
	CHECK(tmp.type == VAL_LIST);
	CHECK(tmp.n_els == 2);
	CHECK(tmp.val.l[0].type == VAL_NUM);
	CHECK(tmp.val.l[0].val.x == 1.0);
	CHECK(tmp.val.l[1].type == VAL_STR);
	CHECK(strcmp(tmp.val.l[1].val.s, "blah") == 0);

	//deallocate memory
	for (_uint i = 0; i < n_args; ++i) free(post_sim_argv_c[i]);
	free(post_sim_argv_c);
	free(sim_argv_c);
	cleanup_settings(&args);
    }

    free(name_dup);
}

TEST_CASE("geometry file reading") {
    parse_ercode er;
    //load settings from the configuration file
    parse_settings args;
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

    SUBCASE("reading of susceptibilities") {
	CHECK(er == E_SUCCESS);
	composite_object* root = geometry.problem.get_roots()[0];
	std::vector<drude_suscept> sups = geometry.parse_susceptibilities(root->fetch_metadata("susceptibilities"), (int*)(&er));
	CHECK(er == E_SUCCESS);
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

    SUBCASE("reading of field sources") {
#ifdef DEBUG_INFO
	std::vector<source_info> sources = geometry.get_sources();
	CHECK(sources.size() == 2);
	source_info inf = sources[0];
	CHECK(inf.type == SRC_GAUSSIAN);
	CHECK(inf.component == meep::Ey);
	CHECK(inf.wavelen == 1.333333);
	CHECK(inf.amplitude == 7.0);
	CHECK(inf.width == doctest::Approx(3.0));
	CHECK(inf.phase == 0.75);
	CHECK(inf.start_time == doctest::Approx(5.2));
	CHECK(inf.end_time == doctest::Approx(41.2));

	inf = sources[1];
	CHECK(inf.type == SRC_CONTINUOUS);
	CHECK(inf.component == meep::Hz);
	CHECK(inf.wavelen == 1.66);
	CHECK(inf.amplitude == 8.0);
	CHECK(inf.phase == 0.0);
	CHECK(inf.start_time == 0.2);
	CHECK(inf.end_time == 1.2);
	CHECK(inf.width == 0.1);
#endif

	cleanup_settings(&args);
    }
}

void* read_h5_array_raw(const H5::Group& grp, const H5::DataType& ctype, size_t el_size, const std::string name, size_t* n_entries) {
    *n_entries = 0;
    try {
	//find the dataspace for real values
	H5::DataSet dataset = grp.openDataSet(name.c_str());
	H5::DataType datatype = dataset.getDataType();
	H5::DataSpace dataspace = dataset.getSpace();
	size_t n_pts = dataspace.getSimpleExtentNpoints();
	//allocate memory for storage
	void* data = malloc(el_size*n_pts);
	dataset.read(data, ctype);

	*n_entries = n_pts;
	dataspace.close();
	datatype.close();
	dataset.close();
	return data;
    } catch (H5::FileIException& error) {
	error.printErrorStack();
	return NULL;
    } catch (H5::GroupIException& error) {
	error.printErrorStack();
	return NULL;
    }
}

TEST_CASE("monitor loc spans") {
    parse_ercode er;

    //load settings from the configuration file
    parse_settings args;
    std::string name = "tests/span.conf";
    char* name_dup = strdup(name.c_str());
    parse_conf_file(&args, name_dup);
    free(name_dup);

    //try creating the geometry object
    bound_geom geometry(args, &er);
    const size_t SPAN = 4;
    CHECK(er == E_SUCCESS);
    //(2/0.1)^2 == 400
    CHECK(geometry.get_n_monitor_clusters() == 2);
    std::vector<meep::vec> mon_locs = geometry.get_monitor_locs();
    CHECK(mon_locs.size() == SPAN + SPAN*SPAN);
    double y = 1;
    for (_uint i = 0; i < SPAN; ++i) {
	CHECK(mon_locs[i].x() == doctest::Approx(y));
	y += 0.5;
    }
    y = 1;
    for (_uint i = 0; i < SPAN; ++i) {
	double x = 1;
	for (_uint j = 0; j < SPAN; ++j) {
	    CHECK(mon_locs[SPAN + SPAN*i + j].x() == doctest::Approx(x));
	    CHECK(mon_locs[SPAN + SPAN*i + j].y() == doctest::Approx(y));
	    CHECK(mon_locs[SPAN + SPAN*i + j].z() == doctest::Approx(1));
	    x += 0.5;
	}
	y += 0.5;
    }

    cleanup_settings(&args);
}

TEST_CASE("sources") {
    gaussian_src_time_phase gp_1(1.0, 2.0, 0.0, 0.0, 2.0);
    gaussian_src_time_phase gp_2(1.0, 2.0, 0.0, 1.0, 2.0);
    CHECK(gp_1.frequency().real() == doctest::Approx(1.0));
    CHECK(gp_1.get_omega() == doctest::Approx(6.283185307));
    CHECK(gp_1.get_width() == doctest::Approx(2.0));
    CHECK(gp_1.get_fwidth() == doctest::Approx(0.5));
}

TEST_CASE("running with a very small system") {
    parse_ercode er;
    char name_buf[BUF_SIZE];

    //load settings from the configuration file
    parse_settings args;
    std::string name = "tests/run.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
    free(name_dup);
    CHECK(ret == 0);
    
    //make things a little faster because we don't care
    args.pml_thickness = 1.0;
    args.len = 2.0;
    args.resolution = 5.0;

    //try creating the geometry object
    bound_geom geometry(args, &er);
    CHECK(er == E_SUCCESS);

    //check that the problem was initialized correctly
    std::vector<meep::vec> mon_locs = geometry.get_monitor_locs();
    CHECK(mon_locs.size() == 2);
    CHECK(mon_locs[0].x() == 1.0);
    CHECK(mon_locs[0].y() == 1.0);
    CHECK(mon_locs[0].z() == 1.0);
    CHECK(mon_locs[1].x() == 2.0);
    CHECK(mon_locs[1].y() == 4.0);
    CHECK(mon_locs[1].z() == 1.1);

    //make sure that monitor locations were added
    geometry.run(args.out_dir);
    //fetch the field times
    std::vector<std::vector<complex>> field_times = geometry.get_field_times();
    CHECK(field_times.size() > 0);

    //check that writing hdf5 files works
    CHECK(field_times.size() == mon_locs.size());
    geometry.save_field_times(args.out_dir);
    
    //We need to create an hdf5 data type for complex values
    H5::CompType fieldtype(sizeof(complex));
    //for some reason linking insertMember breaks on the cluster, we do it manually
    hid_t float_member_id = H5_float_type.getId();
    //snprintf(name_buf, BUF_SIZE, "Re");
    herr_t ret_val = H5Tinsert(fieldtype.getId(), "Re", HOFFSET(complex, re), float_member_id);
    CHECK(ret_val == 0);
    //snprintf(name_buf, BUF_SIZE, "Im");
    ret_val = H5Tinsert(fieldtype.getId(), "Im", HOFFSET(complex, im), float_member_id);
    CHECK(ret_val == 0);
    //do the same for location types
    H5::CompType loctype(sizeof(sto_vec));
    //for some reason linking insertMember breaks on the cluster, we do it manually
    snprintf(name_buf, BUF_SIZE, "x");
    ret_val = H5Tinsert(loctype.getId(), name_buf, HOFFSET(sto_vec, x), float_member_id);
    CHECK(ret_val == 0);
    snprintf(name_buf, BUF_SIZE, "y");
    ret_val = H5Tinsert(loctype.getId(), name_buf, HOFFSET(sto_vec, y), float_member_id);
    CHECK(ret_val == 0);
    snprintf(name_buf, BUF_SIZE, "z");
    ret_val = H5Tinsert(loctype.getId(), name_buf, HOFFSET(sto_vec, z), float_member_id);
    CHECK(ret_val == 0);

    //read the h5 file
    snprintf(name_buf, BUF_SIZE, "%s/field_samples.h5", args.out_dir);
    H5::H5File file(name_buf, H5F_ACC_RDONLY);
    //check that the info group is correct
    H5::Group grp = file.openGroup("info");
    size_t n_c_pts, n_l_pts, n_f_pts, n_t_pts;
    _ftype* t_bounds = (_ftype*)read_h5_array_raw(grp, H5_float_type, sizeof(_ftype), "time_bounds", &n_t_pts);
    CHECK(n_t_pts == 3);
    CHECK(t_bounds[0] == 0.0);
    CHECK(t_bounds[1] > 0.0);
    CHECK(t_bounds[2] > 0.0);
    CHECK(t_bounds[2] < t_bounds[1]);
    free(t_bounds);
    //read the cluster data
    hsize_t* clust_data = (hsize_t*)read_h5_array_raw(grp, H5::PredType::NATIVE_HSIZE, sizeof(hsize_t), "n_clusters", &n_c_pts);
    CHECK(n_c_pts == 1);
    size_t n_clusts = *clust_data;
    free(clust_data);
    CHECK(n_clusts == geometry.get_n_monitor_clusters());

    //iterate through each of the specified clusters
    size_t n_group_digits = (size_t)(n_clusts / log(10)) + 1;
    for (size_t i = 0; i < n_clusts; ++i) {
	strncpy(name_buf, CLUSTER_NAME, BUF_SIZE);
	write_number(name_buf + strlen(CLUSTER_NAME), BUF_SIZE-strlen(CLUSTER_NAME), i, n_group_digits);
	printf("Now reading cluster %lu\n", i);
	grp = file.openGroup(name_buf);
	//check that the location data is correct
	sto_vec* l_data = (sto_vec*)read_h5_array_raw(grp, loctype, sizeof(sto_vec), "locations", &n_l_pts);
	CHECK(n_l_pts == mon_locs.size());
	for (_uint i = 0; i < mon_locs.size(); ++i) {
	    CHECK(l_data != NULL);
	    CHECK(l_data[i].x == doctest::Approx(mon_locs[i].x()));
	    CHECK(l_data[i].y == doctest::Approx(mon_locs[i].y()));
	    CHECK(l_data[i].z == doctest::Approx(mon_locs[i].z()));
	}
	free(l_data);
	//iterate through each point and read time series
	size_t n_pt_digits = (size_t)(log(n_l_pts) / log(10)) + 1;
	strncpy(name_buf, POINT_NAME, BUF_SIZE);
	//check that the time series and fourier transforms are correct
	for (_uint j = 0; j < n_l_pts; ++j) {
	    //open the appropriate group
	    write_number(name_buf + strlen(POINT_NAME), BUF_SIZE-strlen(POINT_NAME), j, n_pt_digits);
	    printf("\tNow reading point %d\n", j);
	    H5::Group point_grp = grp.openGroup(name_buf);
      
	    //read the data from the file we opened
	    complex* f_data = (complex*)read_h5_array_raw(point_grp, fieldtype, sizeof(complex), "frequency", &n_f_pts);
	    complex* t_data = (complex*)read_h5_array_raw(point_grp, fieldtype, sizeof(complex), "time", &n_t_pts);
	    CHECK(f_data != NULL);
	    CHECK(t_data != NULL);
	    CHECK(n_f_pts > 0);
	    //since the fourier transform should only go to +- the nyquist frequency, it must have fewer elements
	    CHECK(n_t_pts >= n_f_pts);

	    //check that the stored times match the data in the geometry object
	    CHECK(n_t_pts == field_times[j].size());
	    for (_uint k = 0; k < n_t_pts; ++k) {
            CHECK(t_data[k].re == field_times[j][k].re);
            CHECK(t_data[k].im == field_times[j][k].im);
	    }
	    free(f_data);
	    free(t_data);
	}
    }
    file.close();

    cleanup_settings(&args);
}

int main(int argc, char** argv) {
    meep::initialize mpi(argc, argv);

    doctest::Context context;
    context.applyCommandLine(argc, argv);

    // overrides
    context.setOption("no-breaks", true);             // don't break in the debugger when assertions fail

    std::cout << "starting tests!" << std::endl;
    int res = context.run(); // run

    if(context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests
    
    int client_stuff_return_code = 0;
    // your program - if the testing framework is integrated in your production code
    
    return res + client_stuff_return_code; // the result from doctest is propagated here as well
}
