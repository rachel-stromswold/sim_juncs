data(type = "field_source", envelope = "Gaussian(\"Ey\", 1.33, 2.0, 0.2, 5.0, 1.0)") {
    Box([0,0,1], [4,4,1]);
}

data(type = "monitor", locations = [[1.0,1.0,1.0],vec(2.0,1.5,1.1)]) {}

Composite(eps = 3.5) {
    Box([0,0,2], [5,5,5]);
}
