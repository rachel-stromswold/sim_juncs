offset = 0.2
list = [1, 2, 3, 4]
sum_list = list[0] + list[1] + list[2] + list[3]
prod_list = list[0]*list[1]*list[2]*list[3] + offset
acid_test = 100*offset - ((list[2] < list[1]) ? list[2] : list[1] + (list[1] > list[0]) ? list[1] : list[2])

data(name = "foo", entry = "bar,(arr),[blah]", num = 3)
data(type = "field_source", envelope = ["Gaussian", "Ey", 1.333333, 3.0, 0.2, 5.0, 4.2, 7.0]) {
    Box([0,0,.2], [.4,.4,.2]);
}
data(type = "field_source", envelope = ["continuous", "Hz", 1.66, 0.2, 1.2, 0.1, 8.0]) {
    Box([0,0,.2], [.4,.4,.2]);
}

data(type = "monitor", locations = [[0.5,0.5,0.5]]) {}

Composite(eps = 3.5, susceptibilities = [[1.0, 0.48, 68.5971845, "drude"],[8.0, 0.816, 452848600800781300, "drude"]]) {
    union() {
	Box([.4,.4,.4], [0.5,.5,.5]);
	Sphere([.4,.4,.6], .1);
    }
    intersect() {
	difference() {
	    Box([.6,.4,.4], [.6+offset, .6, 0.4+offset]);
	    Plane([.6,.4,.6], [8/sum_list,.4,.6], [.6,.6,.4]);
	}
	Cylinder([0.7,.5,.4], .2, .1);
    }
}
snapshot("/tmp/test.pgm", [3,3,3], look=[-1,-1,-1]);
