//unit amplitude linearly polarized 45 deg.
data(type = "field_source", envelope = ["Gaussian", "Ex", 0.7, 1.0, 1.570796, 0.0, 3.5, 0.707107]) {
    Box([0,0,1], [18,18,1]);
}

data(type = "monitor", locations = [vec(x, 9.0, $TOP-0.1) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.0, $TOP+1.5) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.0, $TOP+3.0) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.0, $TOP+3.5) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.0, $TOP+4.0) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.0, $TOP+4.5) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.0, $TOP+5.0) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.2, $TOP-0.1) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.2, $TOP+1.5) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.2, $TOP+3.0) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.2, $TOP+3.5) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.2, $TOP+4.0) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.2, $TOP+4.5) for x in range(-2,8,0.2)]) {}
data(type = "monitor", locations = [vec(x, 9.2, $TOP+5.0) for x in range(-2,8,0.2)]) {}

//Au
//elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
//wavelength range: 0.2 - 12.4 μm
Composite(eps = 1.0, susceptibilities = [[1e-10, 0.04274738474121455, 4.0314052191361974e21, "drude"],[0.3347200880680007 0.19437961740816426 11.3629356945, "lorentz"]]) {
    union() {
        intersect() {
            Plane([-$LEFT, 0, $TOP], [$LEFT, 9, $BOT], [$LEFT, 9, $TOP]);
            Plane([-$LEFT, 18, $TOP], [$LEFT, 9, $TOP], [$LEFT, 9, $BOT]);
        }
        intersect() {
            Plane([18+$RGHT, 0, $TOP], [$RGHT, 9, $TOP], [$RGHT, 9, $BOT]);
            Plane([18+$RGHT, 18, $TOP], [$RGHT, 9, $BOT], [$RGHT, 9, $TOP]);
        }
    }
    Box([$RGHT, 0, $TOP], [18,    18, $BOT]);
}

//SiO2
//ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
//wavelength range: 0.25 - 1.77 μm
Composite(eps = 1.0, susceptibilities = 9.67865314895427, 0.08065544290795199, 1.12, "lorentz"]]) {
    intersect() {
        union() {
            Plane([-$LEFT, 0, $TOP], [$LEFT, 9, $TOP], [$LEFT, 9, $BOT]);
            Plane([-$LEFT, 18, $TOP], [$LEFT, 9, $BOT], [$LEFT, 9, $TOP]);
        }
        union() {
            Plane([18+$RGHT, 0, $TOP], [$RGHT, 9, $BOT], [$RGHT, 9, $TOP]);
            Plane([18+$RGHT, 18, $TOP], [$RGHT, 9, $TOP], [$RGHT, 9, $BOT]);
        }
    }
}