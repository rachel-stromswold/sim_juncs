data(type = "field_source", envelope = "Gaussian(\"Ey\", 1.333333, 3.0, 0.2, 5.0, 7.0)") {
    Box([0,0,2], [4,4,2]);
}

data(type = "monitor", locations = [vec(x, y, 1) for x in range(1,3,0.5)]) {}
data(type = "monitor", locations = flatten([[vec(x, y, 1) for x in range(1,3,0.5)] for y in range(1,3,0.5)])) {}

Composite(eps = 3.5, susceptibilities = "(1.0,0.48,68.5971845,drude),(8.0,0.816,452848600800781300,drude)") {
    Box([4,4,4], [5,5,5]);
}
