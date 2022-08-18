data(name = foo, entry = "bar,(arr),[blah]")
data(type = field_source, envelope = "Gaussian(Ey, 1.33, 3.0, 0.2, 5.0, 1.0)") {
    Box([0,0,2], [4,4,2]);
}
data(type = field_source, envelope = "continuous(Hz, 1.66, 0.2, 1.2, 0.1, 1.0)") {
    Box([0,0,2], [4,4,2]);
}

Composite(eps = 3.5, susceptibilities = "(1.0,0.48,68.5971845,drude),(8.0,0.816,452848600800781300,drude)") {
    union() {
	Box([4,4,4], [5,5,5]);
	Sphere([4,4,6], 1);
    }
    intersect() {
	difference() {
	    Box([6,4,4], [8,6,6]);
	    Box([7,5,5], [8,6,6]);
	}
	Cylinder([7,5,4], 2, 1);
    }
}
