Composite(eps = 3.5) {
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
