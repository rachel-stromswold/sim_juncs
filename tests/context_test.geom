Gaussian_source("Ex", 1.33, 1.0, 2.0, 0.0, 0.0, cutoff=.125, start_time=-1.25, Box([0,0,1], [4,4,1]))//comments should still work here

CW_source("Hz", 0.625, 0.25, 0.75, end_time=2.25, slowness=12, Box([0,0,0], [1/2, 1, 1]))
Composite(eps = 3.5, alpha=1, color=10, [
    Box([1/2,0,0], [1, 1, 1]),
    Box([1/2,0,2], [2, 2, 3])
])

//make sure that comments work
/*even if they're spread out
on multiple lines*/
snapshot("/tmp/run_alpha.pgm", [1.2, 1.2, 1.2], look=[-1,-1,-1], scale=1, up=[0,0,-1]);
