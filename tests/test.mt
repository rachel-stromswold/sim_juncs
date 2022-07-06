Composite(eps = 3.5) {
    union() {
	Box([0,0,0], [1,1,1]);
	Sphere([2,0,0], 1);
    }
    intersect() {
	union() {
	    Box([2,0,0], [1,1,1]);
	    Box([3,0,0], [1,1,1]);
	}
	Cylinder([2,0,0], 1, 1);
    }
}
