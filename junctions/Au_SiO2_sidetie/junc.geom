meep_thick = um_to_l(thick)
meep_width = um_to_l(width)
mid = length/2
sharp = 1.0
top = mid*(1.0 - sharp)
bot = mid*(1.0 + sharp)
left = (length - meep_width)/2
rght = (length + meep_width)/2
print_res = print("top = ", top, ", bot = ", bot, ", left = ", left, ", right = ", rght)

//unit amplitude linearly polarized 45 deg.
Gaussian_source("Ex", wavelen, 1.0, 1.0, 0.0, cutoff=3.5, Box([0,0,1], [length,length,1]))

//after some thinking, these work out to the same heights used in the box
mon_start = left - um_l(0.08333)
monitors(locations = [vec(x, mid, mid-0.45*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.40*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.35*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.30*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.25*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.20*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.15*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.10*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid-0.05*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.05*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.10*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.15*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.20*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.25*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.30*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.35*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.40*meep_thick) for x in linspace(mon_start, mid, 20)])
monitors(locations = [vec(x, mid, mid+0.45*meep_thick) for x in linspace(mon_start, mid, 20)])

//Au
//elemental metals from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83 (1998)
//wavelength range: 0.2 - 12.4 μm
Composite(eps = 1.0, color=42, susceptibilities = [[1e-10, 0.04274738474121455, 4.0314052191361974e21, "drude"],[0.3347200880680007 0.19437961740816426 11.3629356945, "lorentz"]], [
    Intersect([
	Union([
	    Intersect([
		Plane([left-mid,length,top], [left,length,mid],	[left,0,mid]),
		Plane([left-mid,0,bot], [left,0,mid],	[left,length,mid]),
	    ]),
	    Intersect([
		Plane([rght+mid,length,top], [rght,0,mid],	[rght,length,mid]),
		Plane([rght+mid,0,bot], [rght,length,mid],	[rght,0,mid]),
	    ]),
	]),
	Box([0,0,0],[length,length,length])
    ])
])

//SiO2
//ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
//wavelength range: 0.25 - 1.77 μm
Composite(eps = 1.0, susceptibilities = [[9.67865314895427, 0.08065544290795199, 1.12, "lorentz"]], [
    Intersect([
	Intersect([
	    Union([
		Plane([left-mid,length,top], [left,0,mid],	[left,length,mid]),
		Plane([left-mid,0,bot], [left,length,mid],	[left,0,mid]),
	    ]),
	    Union([
		Plane([rght+mid,length,top], [rght,length,mid],	[rght,0,mid]),
		Plane([rght+mid,0,bot], [rght,0,mid],	[rght,length,mid]),
	    ])
	]),
	Box([0,0,0],[length,length,length])
    ])
])

snap_file = out_dir+"/junc.pgm"
print_res = print("saving snapshot to ", snap_file)
snapshot(snap_file, [0,0,0], look=[length/2,length/2,length/2], resolution=512, scale=12)
