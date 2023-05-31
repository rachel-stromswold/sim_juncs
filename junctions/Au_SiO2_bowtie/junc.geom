meep_thick = um_to_l(thick)
meep_width = um_to_l(width)
mid = length/2
top = (length - meep_thick)/2
bot = (inf_thick)? length : ((length + meep_thick)/2)
left = (length - meep_width)/2
rght = (length + meep_width)/2
print("top = ", top, ", bot = ", bot, ", left = ", left, ", right = ", rght)

pwidth = n_cycles*wavelen/0.29979
Gaussian_source("Ex", wavelen, 1.0, pwidth, 0.0, cutoff=4, Box([0,0,1], [length,length,1]))

mon_start = left - um_to_l(0.08)
dx = 0.025
npt = 40
monitors(locations = [vec(x, mid, top-dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+2*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+3*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+4*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+5*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+6*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+7*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+8*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+9*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+10*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+11*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+12*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+13*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+14*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+15*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+16*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+17*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+18*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+19*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+20*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+21*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+22*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+23*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+24*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+25*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+26*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+27*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+28*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+29*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+30*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+31*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+32*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+33*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+34*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+35*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+36*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+37*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+38*dx*meep_thick) for x in linspace(mon_start, mid, npt)])
monitors(locations = [vec(x, mid, top+39*dx*meep_thick) for x in linspace(mon_start, mid, npt)])

//Au
//elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
//wavelength range: 0.2 - 12.4 μm
Composite(eps = 1.0, color=42, susceptibilities = [[1e-10, 0.04274738474121455, 4.0314052191361974e21, "drude"],[0.3347200880680007 0.19437961740816426 11.3629356945, "lorentz"]], [
    Intersect([
	Union([
	    Intersect([
		Plane([left-mid, 0, 	top], [left, mid, top], [left, mid, bot]),
		Plane([left-mid, length,top], [left, mid, bot], [left, mid, top]),
	    ]),
	    Intersect([
		Plane([rght+mid, 0,	top], [rght, mid, bot], [rght, mid, top]),
		Plane([rght+mid, length,top], [rght, mid, top], [rght, mid, bot]),
	    ])
	]),
	Box([0,0,top], [length,length,bot])
    ])
])

//SiO2
//ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
//wavelength range: 0.25 - 1.77 μm
Composite(eps = 1.0, susceptibilities = [[9.67865314895427, 0.08065544290795199, 1.12, "lorentz"]], [
    Intersect([
	Intersect([
	    Union([
		Plane([left-mid, 0,	top], [left, mid, bot], [left, mid, top]),
		Plane([left-mid, length,top], [left, mid, top], [left, mid, bot]),
	    ]),
	    Union([
		Plane([rght+mid, 0,	top], [rght, mid, top], [rght, mid, bot]),
		Plane([rght+mid, length,top], [rght, mid, bot], [rght, mid, top]),
	    ])
	]),
	Box([0,0,top], [length,length,bot])
    ]),
    Box([0, 0, bot], [length,length,length])
])

snap_file = out_dir+"/junc.pgm"
print_res = print("saving snapshot to ", snap_file)
snapshot(snap_file, [2,2,2], look=[18,18,18], resolution=512, scale=12)
