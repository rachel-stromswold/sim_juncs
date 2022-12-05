#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>
#include "main_test.hpp"

TEST_CASE("Test Fourier transforms") {
    data_arr dat;
    make_data_arr(&dat, 1 << POW_N);

    SUBCASE("Test a very simple sine wave") {
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
    }

    SUBCASE("the fourier transform of a gaussian is a gaussian") {
	const _ftype gauss_sigma_sq = 2;
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
	data_arr fft_dat = rfft(dat);
	complex scale = {1/(_ftype)sqrt(2*M_PI*gauss_sigma_sq), 0};
	pw_mult_scale(fft_dat, scale);
	//print out the result
	printf("Fourier transform Gaussian: \n");
	for (size_t k = 0; k < fft_dat.size; ++k) printf("%f+%fi ", fft_dat[k].re, fft_dat[k].im);
	printf("\n");

	//check that the inverse fft of the fft is approximately the original
	_ftype omega_scale = 2*M_PI/dat.size;
	_ftype k_0_by_sigma = k_0/sqrt(gauss_sigma_sq);
	for (size_t k = 0; k < fft_dat.size; ++k) {
	    _ftype omega = omega_scale*k;
	    //the phase should be -\frac{\sigma^2\omega^2}{2} - i\omega k_0
	    complex expected_phase = {-omega*omega*gauss_sigma_sq/2, -omega*k_0};
	    complex expected = c_exp(expected_phase);
	    CHECK(fft_dat[k].re == doctest::Approx(expected.re));
	    CHECK(fft_dat[k].im == doctest::Approx(expected.im));
	}
	cleanup_data_arr(&fft_dat);
    }

    SUBCASE("the inverse fourier transform of a fourier transform is itself") {
	//initialize the lcg
	_uint rand = 314159265;
	complex phi = {0.0, 0.0};
	complex x = {0.0, 0.0};
	_ftype N = (_ftype)dat.size;
	//generate a pseudo-random walk
	for (int k = 0; k < dat.size; ++k) {
	    //fill up each index with a pseudo random number between 0 and 1
	    dat.buf[k].re = x.re;
	    dat.buf[k].im = x.im;
	    rand = lcg(rand);
	    phi.im = 2*M_PI*(double)rand / (double)LCG_MOD;
	    x += c_exp(phi);
	}
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
    }
    cleanup_data_arr(&dat);
}

TEST_CASE("Test line reading") {
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

TEST_CASE("Check that numbers are written correctly") {
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

TEST_CASE("Test builtin functions") {
    char buf[BUF_SIZE];
    Value tmp_val;
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
}

TEST_CASE("Test value parsing") {
    char buf[BUF_SIZE];
    context sc;
    parse_ercode er = E_SUCCESS;
    Value tmp_val;

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
	    Value element = tmp_val.val.l[0];
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
	    Value element = tmp_val.val.l[0];
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

TEST_CASE("Test function parsing") {
    char buf[BUF_SIZE];

    const char* test_func_1 = "f()";
    const char* test_func_2 = "f(\"a\", \"b\", \"c\", 4)";
    const char* test_func_3 = "foo(vec(1,2,3),\"a\",\"banana\")";
    const char* test_func_4 = "foo(1, \"Box(0,1,2,3)\", 4+5)";
    const char* test_func_5 = "foo ( 1 , \"b , c\" )";
    const char* test_func_6 = "f(eps = 3.5)";
    const char* test_func_7 = "f(name = \"bar\")";
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
    evec3* tmp_vec = cur_func.args[0].get_val().v;
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
    CHECK(strcmp(cur_func.args[1].to_c_str(), "Box(0,1,2,3)") == 0);
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

TEST_CASE("Test Object Trees") {
    //declare variables
    char buf[BUF_SIZE];
    ObjectStack test_stack;
    Object* cur_obj;object_type cur_type;
    parse_ercode er;

    //setup a bunch of strings describing objects
    const char* root_obj_str = "Composite(eps = 3.5)";
    const char* l_str = "union()";
    const char* ll_str = "Box(vec(0,0,0), vec(1,1,1))";
    const char* lr_str = "Sphere(vec(2,0,0), 1)";
    const char* r_str = "intersect()";
    const char* rl_str = "Box(vec(0,0,0), vec(1,1,1))";
    const char* rr_str = "Cylinder(vec(2,0,0), 1, 1)";

    Scene sc;
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
    CHECK(data_vec[0]->has_metadata("num"));
    Value name_val = data_vec[0]->fetch_metadata("name");
    Value ntry_val = data_vec[0]->fetch_metadata("entry");
    Value num_val = data_vec[0]->fetch_metadata("num");
    CHECK(name_val.type == VAL_STR);
    CHECK(ntry_val.type == VAL_STR);
    CHECK(num_val.type == VAL_NUM);
    CHECK(name_val.val.s != NULL);
    CHECK(ntry_val.val.s != NULL);
    CHECK(strcmp(name_val.to_c_str(), "foo") == 0);
    CHECK(strcmp(ntry_val.to_c_str(), "bar,(arr),[blah]") == 0);
    CHECK(num_val.to_float() == 3);

    //test geometric information
    std::vector<CompositeObject*> roots_vec = s.get_roots();
    CHECK(roots_vec.size() > 0);
    CompositeObject* root = roots_vec[0];

    CHECK(root != NULL);
    //check that all are not NULL and that types are correct
    CHECK(root->get_child_type_l() == CGS_COMPOSITE);
    CHECK(root->get_child_type_r() == CGS_COMPOSITE);
    CompositeObject* comp_l = (CompositeObject*)(root->get_child_l());
    CompositeObject* comp_r = (CompositeObject*)(root->get_child_r());
    CHECK(comp_r->get_child_type_l() == CGS_COMPOSITE);
    CompositeObject* comp_rl = (CompositeObject*)(comp_r->get_child_l());
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
    CHECK(comp_rl->get_child_type_r() == CGS_PLANE);
    CHECK(comp_r->get_child_r() != NULL);
    CHECK(comp_r->get_child_type_r() == CGS_CYLINDER);
}

TEST_CASE("Test Geometric Inclusion") {
    parse_ercode er;
    Scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    CompositeObject* root = s.get_roots()[0];

    CHECK(root->in(Eigen::Vector3d(.45,.45,.45)) == 1);
    CHECK(root->in(Eigen::Vector3d(.45,.41,.6)) == 1);
    CHECK(root->in(Eigen::Vector3d(.65,.45,.41)) == 1);
    CHECK(root->in(Eigen::Vector3d(.71,.51,.51)) == 0);
    CHECK(root->in(Eigen::Vector3d(.55,.45,.45)) == 0);
    CHECK(root->in(Eigen::Vector3d(.55,.41,.85)) == 0);
}

TEST_CASE("Test volumes") {
    //initialize random state
    _uint state = lcg(lcg(TEST_SEED));
    double x, y, z;
    //initialize the camera for viewing the volumes
    evec3 cam_pos(CAM_X,CAM_Y,CAM_Z);
    evec3 y_comp(-CAM_X*CAM_Z, -CAM_Y*CAM_Z, CAM_X*CAM_X+CAM_Y*CAM_Y);
    evec3 x_comp = -cam_pos.cross(y_comp);
    emat3 view_mat;
    view_mat << x_comp/x_comp.norm(), y_comp/y_comp.norm(), -cam_pos/cam_pos.norm();
    //figure out where the origin gets mapped. We will make this the center of the screen
    evec3 offset = -view_mat*cam_pos;
    offset.z() = 0;
    //test volumes by sampling from a cube with side lengths unit 1 and testing the fraction that are included
    evec3 center(0.5, 0.5, 0.5);
    evec3 corner_1(1, 1, 0);
    evec3 corner_2(0, 1, 1);
    evec3 corner_3(1, 0, 1);
    Sphere test_sphere(center, 0.5);
    Cylinder test_cyl(center, 0.5, 0.5, 0.5);
    Plane test_plane(corner_1, corner_2, corner_3);
    double sphere_frac = 0;
    double plane_frac = 0;
    double cyl_frac = 0;
    //now load the scene just so we can draw it
    parse_ercode er;
    Scene s("tests/test.geom", &er);
    CHECK(er == E_SUCCESS);
    CompositeObject* test_root = s.get_roots()[0];
    //these are matrices specifying the intensity of each pixel in grayscale
    _uint8 sphere_arr[IM_RES*IM_RES];
    _uint8 plane_arr[IM_RES*IM_RES];
    _uint8 cyl_arr[IM_RES*IM_RES];
    _uint8 scene_arr[IM_RES*IM_RES];
    for (size_t i = 0; i < IM_RES*IM_RES; ++i) { sphere_arr[i]=0xff;plane_arr[i]=0xff;cyl_arr[i]=0xff;scene_arr[i]=0xff; }
    //sample random points
    for (size_t i = 0; i < TEST_N; ++i) {
	state = lcg(state);
	x = (double)state / LCG_MOD;
	state = lcg(state);
	y = (double)state / LCG_MOD;
	state = lcg(state);
	z = (double)state / LCG_MOD;
	state = lcg(state);
	//initialize a point from the random values and compute it's position projected onto the view matrix. The z-buffer will act as color intensity
	evec3 r(x, y, z);
	evec3 view_r = view_mat*(r - cam_pos);
	size_t ind = floor( (1+view_r.x())*IM_RES/2 ) + IM_RES*floor( (1+view_r.y())*IM_RES/2 );
	if (ind >= IM_RES*IM_RES) ind = 0;
	_uint8 depth = floor(IM_DPT*view_r.z());
	if (test_sphere.in(r)) {
	    sphere_frac += 1.0/TEST_N;
	    if (depth < sphere_arr[ind]) sphere_arr[ind] = depth;
	}
	if (test_plane.in(r)) {
	    plane_frac += 1.0/TEST_N;
	    if (depth < plane_arr[ind]) plane_arr[ind] = depth;
	}
	if (test_cyl.in(r)) {
	    cyl_frac += 1.0/TEST_N;
	    if (depth < cyl_arr[ind]) cyl_arr[ind] = depth;
	}
	if (test_root->in(r) && depth < scene_arr[ind]) {
	    scene_arr[ind] = depth;
	}
    }
    //write pgm files to visualize the geometries being drawn
    FILE* fp_sphere = fopen("/tmp/sphere.pgm", "w");
    FILE* fp_plane = fopen("/tmp/plane.pgm", "w");
    FILE* fp_cyl = fopen("/tmp/cyl.pgm", "w");
    FILE* fp_scene = fopen("/tmp/scene.pgm", "w");
    fprintf(fp_sphere, "P2\n%d %d\n%d\n", IM_RES, IM_RES, IM_DPT);
    fprintf(fp_plane, "P2\n%d %d\n%d\n", IM_RES, IM_RES, IM_DPT);
    fprintf(fp_cyl, "P2\n%d %d\n%d\n", IM_RES, IM_RES, IM_DPT);
    fprintf(fp_scene, "P2\n%d %d\n%d\n", IM_RES, IM_RES, IM_DPT);
    for (size_t yy = 0; yy < IM_RES; ++yy) {
	for (size_t xx = 0; xx < IM_RES; ++xx) {
	    size_t ind = yy*IM_RES+xx;
	    fprintf(fp_sphere, "%d ", sphere_arr[ind]);
	    fprintf(fp_plane, "%d ", plane_arr[ind]);
	    fprintf(fp_cyl, "%d ", cyl_arr[ind]);
	    fprintf(fp_scene, "%d ", scene_arr[ind]);
	}
	fprintf(fp_sphere, "\n");
	fprintf(fp_plane, "\n");
	fprintf(fp_cyl, "\n");
	fprintf(fp_scene, "\n");
    }
    fclose(fp_sphere);
    fclose(fp_plane);
    fclose(fp_cyl);
    fclose(fp_scene);

    double v_sphere = 1.33333*M_PI*0.125;
    double v_plane = 5.0/6;
    double v_cyl = M_PI*0.125;
    
    CHECK(abs(v_sphere - sphere_frac)/v_sphere < EPSILON);
    CHECK(abs(v_plane - plane_frac)/v_plane < EPSILON);
    CHECK(abs(v_cyl - cyl_frac)/v_cyl < EPSILON);
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
	//CHECK(strcmp(args.monitor_locs, "(1.0,1.0,1.0)") == 0);
	CHECK(args.post_source_t == 1.0);
	CHECK(args.field_dump_span == 171);
	CHECK(args.ambient_eps == 1.0);
	CHECK(strcmp(args.geom_fname, "tests/test.geom") == 0);

	cleanup_settings(&args);
    }

    SUBCASE("Command line arguments override defaults") {
	Settings args;
	//parse commandline arguments
	std::string sim_argv_cpp[] = { "./test", "--conf-file", "blah.conf", "--geom-file", "blah.geom", "--out-dir", "/blah", "--grid-res", "3.0", "--length", "9.0", "--eps1", "2.0" };
	//gcc doesn't like using string literals as c-strings, ugh.
	int n_args = sizeof(sim_argv_cpp)/sizeof(std::string);
	//need to create two lists since parse_args modifies the list in place. The second list is only used so that we know which addresses to free.
	char** sim_argv_c = (char**)malloc(sizeof(char*)*n_args);
	char** post_sim_argv_c = (char**)malloc(sizeof(char*)*n_args);
	for (_uint i = 0; i < n_args; ++i) {
	    sim_argv_c[i] = strdup(sim_argv_cpp[i].c_str());
	    post_sim_argv_c[i] = sim_argv_c[i];
	}

	//finally we can parse the command line arguments
	int post_n_args = n_args;
	int ret = parse_args(&args, &post_n_args, sim_argv_c);
	CHECK(post_n_args == 1);

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
	//CHECK(strcmp(args.monitor_locs, "(1.0,1.0,1.0)") == 0);
	CHECK(args.post_source_t == 1.0);

	//deallocate memory
	for (_uint i = 0; i < n_args; ++i) free(post_sim_argv_c[i]);
	free(post_sim_argv_c);
	free(sim_argv_c);
	cleanup_settings(&args);
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
	char* dat = strdup(root->fetch_metadata("susceptibilities").to_c_str());
	Value dat_val = make_val_str(dat);
	std::vector<drude_suscept> sups = geometry.parse_susceptibilities(dat_val, (int*)(&er));
	cleanup_val(&dat_val);
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
    std::vector<source_info> sources = geometry.get_sources();
	CHECK(sources.size() == 2);
	source_info inf = sources[0];
	CHECK(inf.type == SRC_GAUSSIAN);
	CHECK(inf.component == meep::Ey);
	CHECK(inf.wavelen == 1.333333);
	CHECK(inf.width == doctest::Approx(3.0));
	CHECK(inf.phase == 0.2);
	CHECK(inf.start_time == 5.0);
	CHECK(inf.end_time == doctest::Approx(30.2));
	CHECK(inf.amplitude == 7.0);

	inf = sources[1];
	CHECK(inf.type == SRC_CONTINUOUS);
	CHECK(inf.component == meep::Hz);
	CHECK(inf.wavelen == 1.66);
	CHECK(inf.phase == 0.0);
	CHECK(inf.start_time == 0.2);
	CHECK(inf.end_time == 1.2);
	CHECK(inf.width == 0.1);
	CHECK(inf.amplitude == 8.0);
#endif

	cleanup_settings(&args);
    }
}

void* read_h5_array_raw(const H5::Group& grp, const H5::DataType& ctype, size_t el_size, const std::string name, size_t* n_entries) {
    *n_entries = 0;
    try {
	//find the dataspace for real values
	H5::DataSet dataset = grp.openDataSet(name);
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
    } catch (H5::FileIException error) {
	error.printErrorStack();
	return NULL;
    } catch (H5::GroupIException error) {
	error.printErrorStack();
	return NULL;
    }
}

TEST_CASE("Test that spans of monitor locations are read correctly") {
    parse_ercode er;

    //load settings from the configuration file
    Settings args;
    std::string name = "tests/span.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
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

TEST_CASE("Test running with a very small system") {
    parse_ercode er;
    char name_buf[BUF_SIZE];

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
    args.resolution = 5.0;

    //try creating the geometry object
    bound_geom geometry(args, &er);
    CHECK(er == E_SUCCESS);

    //make sure that monitor locations were added
    geometry.run(args.out_dir);
    //fetch the field times
    std::vector<data_arr> field_times = geometry.get_field_times();
    CHECK(field_times.size() > 0);

    //check that writing hdf5 files works
    std::vector<meep::vec> mon_locs = geometry.get_monitor_locs();
    CHECK(mon_locs.size() > 0);
    CHECK(mon_locs.size() == field_times.size());
    geometry.save_field_times(args.out_dir);
    
    //We need to create an hdf5 data type for complex values
    H5::CompType fieldtype(sizeof(complex));
    //for some reason linking insertMember breaks on the cluster, we do it manually
    hid_t float_member_id = H5_float_type.getId();
    snprintf(name_buf, BUF_SIZE, "Re");
    herr_t ret_val = H5Tinsert(fieldtype.getId(), name_buf, HOFFSET(complex, re), float_member_id);
    CHECK(ret_val == 0);
    snprintf(name_buf, BUF_SIZE, "Im");
    ret_val = H5Tinsert(fieldtype.getId(), name_buf, HOFFSET(complex, im), float_member_id);
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
	printf("Now reading cluster %d\n", i);
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
	    CHECK(n_t_pts == field_times[j].size);
	    for (_uint k = 0; k < n_t_pts; ++k) {
		CHECK(t_data[k].re == field_times[j].buf[k].re);
		CHECK(t_data[k].im == field_times[j].buf[k].im);
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

    int res = context.run(); // run

    if(context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests
    
    int client_stuff_return_code = 0;
    // your program - if the testing framework is integrated in your production code
    
    return res + client_stuff_return_code; // the result from doctest is propagated here as well
}
