data(name = foo, entry = "bar,(arr),[blah]")

Composite(eps = 3.5, susceptibilities = "(1.0,0.48,68.5971845,drude),(8.0,0.816,452848600800781300,drude)") {
    union() {
	Box([0,0,0], [1,1,1]);
	Sphere([0,0,2], 1);
    }
    intersect() {
	difference() {
	    Box([2,0,0], [2,2,2]);
	    Box([3,1,1], [1,1,1]);
	}
	Cylinder([3,1,0], 2, 1);
    }
}
