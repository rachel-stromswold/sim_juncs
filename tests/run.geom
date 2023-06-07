tot_len = sim_length + 2*pml_thickness

Gaussian_source("Ey", 1.33, 1.0, 2.0, 0.2, 5.0, Box([0,0,1], [4,4,1]))

monitors(locations = [[1.0,1.0,1.0],vec(um_to_l(1),tot_len,1.1)])

Composite(eps = 3.5, color=90, {
    Box([0,0,0], [tot_len/2-0.1,tot_len,tot_len]);
    Box([tot_len/2+0.1,0,0], [tot_len,tot_len,tot_len]);
    snapshot("/tmp/run.pgm", [3,3,3], look=[-1,-1,-1]);
})
a = print("left box goes until ", tot_len/2-0.1, " right box goes until", tot_len/2+0.1)
