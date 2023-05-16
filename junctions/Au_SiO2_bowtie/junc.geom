meep_thick = um_to_l(thick)
meep_width = um_to_l(width)
mid = length/2
top = (length - meep_thick)/2
bot = (length + meep_thick)/2
left = (length - meep_width)/2
rght = (length + meep_width)/2
print("top = ", top, ", bot = ", bot, ", left = ", left, ", right = ", rght)

//unit amplitude linearly polarized 45 deg.
Gaussian_source("Ex", wavelen, 1.0, 1.0, 0.0, cutoff=5, Box([0,0,1], [length,length,1]))

mon_start = left - um_l(0.08333)
dx = 0.05
monitors(locations = [vec(x, mid, top-dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, mid, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
ofa = mid + um_to_l(0.008333)
monitors(locations = [vec(x, ofa, top-dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])
monitors(locations = [vec(x, ofa, top+dx*meep_thick) for x in linspace(mon_start, mid, dx)])

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
