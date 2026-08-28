import gmsh
import sys

num_threads = 0 # 0 for all available

# geometry
D1 = 0.005
D2 = 0.01
D3 = 0.3
D4 = 0.6
D5 = 3
L1 = 1.25
L2 = 0.7
L3 = 5
L4 = 20

dx = 10
dy = 10
dz = 30
cx = 5

# mesh sizes
base_size = 0.8
supports_size = 0.12
ends_size = 0.003
wire_size = 0.001
domain_size = 0.8

growth_rate = 0.15

gmsh.initialize()
gmsh.model.add('sonda')

occ = gmsh.model.occ

gmsh.option.set_string('Geometry.OCCTargetUnit', 'M')
shapedimtags = occ.import_shapes('sonda.step')
occ.fragment(shapedimtags, [])

occ.synchronize()

eps = 1e-6
base = [dimtag[1] for dimtag in gmsh.model.get_entities_in_bounding_box(-D5/2-eps, -D5/2-eps, 0-eps, D5/2+eps, D5/2+eps, L4+eps, 3)]
supports = [dimtag[1] for dimtag in gmsh.model.get_entities_in_bounding_box(-D5/2-eps, -D5/2-eps, L4-eps, D5/2+eps, D5/2+eps, L4+L3+eps, 3)]
end1 = [dimtag[1] for dimtag in gmsh.model.get_entities_in_bounding_box(-D5/2, -D5/2, L4+L3-D2/2-eps, D5/2, -L1/2+eps, L4+L3+D2/2+eps, 3)]
end2 = [dimtag[1] for dimtag in gmsh.model.get_entities_in_bounding_box(-D5/2, L1/2-eps, L4+L3-D2/2-eps, D5/2, D5/2, L4+L3+D2/2+eps, 3)]
wire = [dimtag[1] for dimtag in gmsh.model.get_entities_in_bounding_box(-D5/2, -L1/2-eps, L4+L3-D1/2-eps, D5/2, L1/2+eps, L4+L3+D1/2+eps, 3)]

all_vols = [v[1] for v in gmsh.model.get_entities(3)]
sonda_vols = base + supports + end1 + end2 + wire
domain = [v for v in all_vols if v not in sonda_vols]

base_surfs = [s[1] for s in gmsh.model.get_boundary([(3, v) for v in base])]
supports_surfs = [s[1] for s in gmsh.model.get_boundary([(3, v) for v in supports])]
ends_surfs = [s[1] for s in gmsh.model.get_boundary([(3, v) for v in end1 + end2])]
wire_surfs = [s[1] for s in gmsh.model.get_boundary([(3, v) for v in wire])]

inlet_surf = [x[1] for x in gmsh.model.get_entities_in_bounding_box(-cx-eps, -dy/2-eps, -eps, -cx+eps, dy/2+eps, dz+eps, 2)]
outlet_surf = [x[1] for x in gmsh.model.get_entities_in_bounding_box(dx-cx-eps, -dy/2-eps, -eps, dx-cx+eps, dy/2+eps, dz+eps, 2)]
bottom_surf = [x[1] for x in gmsh.model.get_entities_in_bounding_box(-cx-eps, -dy/2-eps, -eps, dx-cx+eps, dy/2+eps, eps, 2)]
top_surf = [x[1] for x in gmsh.model.get_entities_in_bounding_box(-cx-eps, -dy/2-eps, dz-eps, dx-cx+eps, dy/2+eps, dz+eps, 2)]
side_ym_surf = [x[1] for x in gmsh.model.get_entities_in_bounding_box(-cx-eps, -dy/2-eps, -eps, dx-cx+eps, -dy/2+eps, dz+eps, 2)]
side_yp_surf = [x[1] for x in gmsh.model.get_entities_in_bounding_box(-cx-eps, dy/2-eps, -eps, dx-cx+eps, dy/2+eps, dz+eps, 2)]

gmsh.model.add_physical_group(2, base_surfs, 100, name='Base_surf')
gmsh.model.add_physical_group(2, supports_surfs, 101, name='Supports_surf')
gmsh.model.add_physical_group(2, ends_surfs, 102, name='Ends_surf')
gmsh.model.add_physical_group(2, wire_surfs, 103, name='Wire_surf')

gmsh.model.add_physical_group(2, inlet_surf, 200, name='Inlet')
gmsh.model.add_physical_group(2, outlet_surf, 201, name='Outlet')
gmsh.model.add_physical_group(2, bottom_surf, 202, name='Bottom_wall')
gmsh.model.add_physical_group(2, top_surf, 203, name='Top_wall')
gmsh.model.add_physical_group(2, side_ym_surf+side_yp_surf, 204, name='Sides')

gmsh.model.add_physical_group(3, base, 300, name='Base_vol')
gmsh.model.add_physical_group(3, supports, 301, name='Supports_vol')
gmsh.model.add_physical_group(3, end1+end2, 302, name='Ends_vol')
gmsh.model.add_physical_group(3, wire, 303, name='Wire_vol')
gmsh.model.add_physical_group(3, domain, 400, name='Domain_vol')

def add_vol_size_field(vols, size):
    f = gmsh.model.mesh.field.add("Constant")
    gmsh.model.mesh.field.set_numbers(f, "VolumesList", vols)
    gmsh.model.mesh.field.set_number(f, "VIn", size)
    gmsh.model.mesh.field.set_number(f, "VOut", 1) 
    return f

def add_growth_field(surfs, size_min, growth_rate):
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.set_numbers(f_dist, "SurfacesList", surfs)
    
    f_eval = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.set_string(f_eval, "F", f"{size_min} + {growth_rate} * F{f_dist}")
    return f_eval

f_base_eval = add_growth_field(base_surfs, base_size, growth_rate)
f_supp_eval = add_growth_field(supports_surfs, supports_size, growth_rate)
f_ends_eval = add_growth_field(ends_surfs, ends_size, growth_rate)
f_wire_eval = add_growth_field(wire_surfs, wire_size, growth_rate)

f_domain = gmsh.model.mesh.field.add("Constant")
gmsh.model.mesh.field.set_number(f_domain, "VIn", domain_size)
gmsh.model.mesh.field.set_number(f_domain, "VOut", domain_size)

f_base_vol = add_vol_size_field(base, base_size)
f_supp_vol = add_vol_size_field(supports, supports_size)
f_ends_vol = add_vol_size_field(end1 + end2, ends_size)
f_wire_vol = add_vol_size_field(wire, wire_size)

f_min = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.set_numbers(f_min, "FieldsList", [
    f_base_eval, f_supp_eval, f_ends_eval, f_wire_eval, f_domain,
    f_base_vol, f_supp_vol, f_ends_vol, f_wire_vol
])
gmsh.model.mesh.field.set_as_background_mesh(f_min)

gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)

gmsh.option.set_number("General.NumThreads", num_threads)
gmsh.model.mesh.generate(3)

gmsh.write('mesh.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()