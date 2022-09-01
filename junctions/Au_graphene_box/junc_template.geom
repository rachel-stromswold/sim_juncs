data(type = field_source, envelope = "Gaussian(Ex, 0.75, 1.5, 0.0, 3.0, 1.0)") {
    Box([0,0,1], [18,18,1]);
}

data(type = monitor, locations = "(2.0, 9.0, 2.0),($LEFT-0.1, 9.0, 2.0),($LEFT+0.1, 9.0, 2.0),(9.0, 9.0, 2.0),(2.0, 9.0, $TOP-0.1),($LEFT-0.1, 9.0, $TOP-0.1),($LEFT+0.1, 9.0, $TOP-0.1),(9.0, 9.0, $TOP-0.1),(2.0, 9.0, 9.0),($LEFT-0.1, 9.0, 9.0),($LEFT+0.1, 9.0, 9.0),(9.0, 9.0, 9.0)") {
}

//Au
//elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
//wavelength range: 0.2 - 12.4 μm
Composite(eps = 1.0, susceptibilities = "(1e-10, 0.04274738474121455, 4.0314052191361974e21, drude),(0.3347200880680007 0.19437961740816426 11.362935694585572)") {
    Box([0,     0, $TOP], [18, $LEFT, $BOT]);
    Box([0, $RGHT, $TOP], [18,    18, $BOT]);
}

//Graphene
//ref: See comments by John Weiner https://github.com/NanoComp/meep/issues/1230
//wavelength range: 0.25 - 1.77 μm
Composite(make_2d=1, eps=1.0, susceptibilities="(1.0, 0.1, 2.88575e31, lorentz)") {
    Box([0, $LEFT, $DELTA_TOP], [18, $RGHT, $DELTA_BOT]);
}

//SiO2
//ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
//wavelength range: 0.25 - 1.77 μm
Composite(eps = 1.0, susceptibilities = "(9.67865314895427, 0.08065544290795199, 1.12, lorentz)") {
    Box([0, 0, $BOT], [18, 18, 18]);
}
