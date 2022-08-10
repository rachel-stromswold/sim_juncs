data(type = field_source, envelope = "Gaussian(Ex, 1.33, 3.5, 0.0, 10.0, 1.0)") {
    Box([0,0,1], [18,18,1]);
}

//Au
//elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
//wavelength range: 0.2 - 12.4 μm
Composite(eps = 1.0, susceptibilities = "(1e-10, 0.04274738474121455, 4.0314052191361974e21, drude),(0.3347200880680007 0.19437961740816426 11.362935694585572)") {
    Box([1,     1, 7.9], [17, 8.975, 8.1]);
    Box([1, 9.025, 7.9], [17,    17, 8.1]);
}

//SiO2
//ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
//wavelength range: 0.25 - 1.77 μm
Composite(eps = 1.0, susceptibilities = "(9.67865314895427, 0.08065544290795199, 1.12, lorentz)") {
    Box([1, 8.975, 7.9], [16, 9.025, 8.1]);
}