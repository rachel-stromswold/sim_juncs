data(name = foo, entry = "bar,(arr),[blah]")

Composite(eps = 3.5, susceptibilities = "(0.0000000001,0.48,6859718450000000000000,drude),(0.065,0.816,6859718450000000000000,drude)") {
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
