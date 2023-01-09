meep_thick = um_to_l(thick)
meep_width = um_to_l(width)
mid = length/2
top = (length - meep_thick)/2
bot = (length + meep_thick)/2
left = (length - meep_width)/2
rght = (length + meep_width)/2
print_res = print("top = ", top, ", bot = ", bot, ", left = ", left, ", right = ", rght)

//unit amplitude linearly polarized 45 deg.
data(type = "field_source", envelope = ["Gaussian", "Ex", 0.7, 1.0, 1.570796, 0.0, 3.5, 0.707107]) {
    Box([0,0,1], [18,18,1]);
}

data(type = "monitor", locations = [vec(x, mid, top-0.1*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid, top+0.1*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid, top+0.3*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid, top+0.5*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid, top+0.7*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid, top+0.9*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid+0.2, top-0.1*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid+0.2, top+0.1*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid+0.2, top+0.3*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid+0.2, top+0.5*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid+0.2, top+0.7*meep_thick) for x in linspace(left-2, mid, 20)]) {}
data(type = "monitor", locations = [vec(x, mid+0.2, top+0.9*meep_thick) for x in linspace(left-2, mid, 20)]) {}

//Au
//elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
//wavelength range: 0.2 - 12.4 μm
Composite(eps = 1.0, color=42, susceptibilities = [[1e-10, 0.04274738474121455, 4.0314052191361974e21, "drude"],[0.3347200880680007 0.19437961740816426 11.3629356945, "lorentz"]]) {
    intersect() {
	union() {
	    intersect() {
		Plane([left-mid, 0, 	top], [left, mid, top], [left, mid, bot]);
		Plane([left-mid, length,top], [left, mid, bot], [left, mid, top]);
	    }
	    intersect() {
		Plane([rght+mid, 0,	top], [rght, mid, bot], [rght, mid, top]);
		Plane([rght+mid, length,top], [rght, mid, top], [rght, mid, bot]);
	    }
	}
	Box([0,0,top], [length,length,bot]);
    }
}

//SiO2
//ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
//wavelength range: 0.25 - 1.77 μm
Composite(eps = 1.0, susceptibilities = 9.67865314895427, 0.08065544290795199, 1.12, "lorentz"]]) {
    intersect() {
	intersect() {
	    union() {
		Plane([left-mid, 0,	top], [left, mid, bot], [left, mid, top]);
		Plane([left-mid, length,top], [left, mid, top], [left, mid, bot]);
	    }
	    union() {
		Plane([rght+mid, 0,	top], [rght, mid, top], [rght, mid, bot]);
		Plane([rght+mid, length,top], [rght, mid, bot], [rght, mid, top]);
	    }
	}
	Box([0,0,top], [length,length,bot]);
    }
    Box([0, 0, bot], [length,length,length]);
}

snap_file = out_dir+"/tri_junc.pgm"
print_res = print("saving snapshot to ", snap_file)
snapshot(snap_file, [2,2,2], look=[18,18,18], resolution=512, scale=12)
