data(type = field_source, envelope = "Gaussian(Ex, 0.7, 1.0, 0.0, 3.5, 1.0)") {
    Box([0,0,1], [18,18,1]);
}

data(type = monitor, locations = "(0.5, 9.0, $TOP+0.1),(1.0, 9.0, $TOP+0.1),(1.5, 9.0, $TOP+0.1),(2.0, 9.0, $TOP+0.1),(2.5, 9.0, $TOP+0.1),(3.0, 9.0, $TOP+0.1),(3.5, 9.0, $TOP+0.1),(4.0, 9.0, $TOP+0.1),(4.5, 9.0, $TOP+0.1),(5.0, 9.0, $TOP+0.1),(5.5, 9.0, $TOP+0.1),(6.0, 9.0, $TOP+0.1),(7.5, 9.0, $TOP+0.1),(8.0, 9.0, $TOP+0.1),(8.5, 9.0, $TOP+0.1),(9.0, 9.0, $TOP+0.1)") {
}

data(type = monitor, locations = "(0.5, 9.0, 9.0),(1.0, 9.0, 9.0),(1.5, 9.0, 9.0),(2.0, 9.0, 9.0),(2.5, 9.0, 9.0),(3.0, 9.0, 9.0),(3.5, 9.0, 9.0),(4.0, 9.0, 9.0),(4.5, 9.0, 9.0),(5.0, 9.0, 9.0),(5.5, 9.0, 9.0),(6.0, 9.0, 9.0),(7.5, 9.0, 9.0),(8.0, 9.0, 9.0),(8.5, 9.0, 9.0),(9.0, 9.0, 9.0)") {
}

data(type = monitor, locations = "(0.5, 9.0, $BOT-0.1),(1.0, 9.0, $BOT-0.1),(1.5, 9.0, $BOT-0.1),(2.0, 9.0, $BOT-0.1),(2.5, 9.0, $BOT-0.1),(3.0, 9.0, $BOT-0.1),(3.5, 9.0, $BOT-0.1),(4.0, 9.0, $BOT-0.1),(4.5, 9.0, $BOT-0.1),(5.0, 9.0, $BOT-0.1),(5.5, 9.0, $BOT-0.1),(6.0, 9.0, $BOT-0.1),(7.5, 9.0, $BOT-0.1),(8.0, 9.0, $BOT-0.1),(8.5, 9.0, $BOT-0.1),(9.0, 9.0, $BOT-0.1)") {
}

//Au
//elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
//wavelength range: 0.2 - 12.4 μm
Composite(eps = 1.0, susceptibilities = "(1e-10, 0.04274738474121455, 4.0314052191361974e21, drude),(0.3347200880680007 0.19437961740816426 11.3629356945$LEFT-0.172)") {
    Box([$LEFT, 0, $TOP], [$LEFT, 18, $BOT]);
    Box([$RGHT, 0, $TOP], [18,    18, $BOT]);
}

//SiO2
//ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
//wavelength range: 0.25 - 1.77 μm
Composite(eps = 1.0, susceptibilities = "(9.67$LEFT-0.1314$LEFT-0.1427, 0.08065544290795199, 1.12, lorentz)") {
    Box([$LEFT, 0, $TOP], [$RGHT, 18, $BOT]);
    Box([0, 0, $BOT], [18, 18, 18]);
}
