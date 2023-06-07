Gaussian_source(\"Ey\", 1.33, 2.0, 0.2, 5.0, 1.0, {
    Box([0,0,1], [4,4,1]);
})

Composite(eps = 3.5, alpha=254, color=22, {
    Box([0,0,0], [1/2, 1, 1]);
})
Composite(eps = 3.5, alpha=1, color=10, {
    Box([1/2,0,0], [1, 1, 1]);
})
snapshot("/tmp/run_alpha.pgm", [1.2, 1.2, 1.2], look=[-1,-1,-1], scale=1, up=[0,0,-1]);
