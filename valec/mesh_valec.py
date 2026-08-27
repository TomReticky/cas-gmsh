import gmsh
import sys
import math

num_threads = 0 # 0 for all available

# geometry
dx = 1.0
dy = 0.5
dz = 0.5
cx = 0.3
h = 0.15
R = 0.05
r = 0.01

# mesh sizes
size_domain = 0.2
size_body = 0.01
size_bnd = 0.05
size_wake = 0.05
dist_min_body = 0.02
dist_max_body = 0.1
dist_min_bnd = 0.05
dist_max_bnd = 0.1

gmsh.initialize()
gmsh.model.add('cylinder')

gmsh.option.set_string('Geometry.OCCTargetUnit', 'M')

occ = gmsh.model.occ

# builds the cylinder
cy = dy / 2.0
p1 = occ.addPoint(cx, cy, 0)
p2 = occ.addPoint(cx + R, cy, 0)
p3 = occ.addPoint(cx + R, cy, h - r)
p4 = occ.addPoint(cx + R - r, cy, h)
p5 = occ.addPoint(cx, cy, h)
pc = occ.addPoint(cx + R - r, cy, h - r)

l1 = occ.addLine(p1, p2)
l2 = occ.addLine(p2, p3)
c1 = occ.addCircleArc(p3, pc, p4)
l3 = occ.addLine(p4, p5)
l4 = occ.addLine(p5, p1)

cl = occ.addCurveLoop([l1, l2, c1, l3, l4])
face1 = occ.addPlaneSurface([cl])
face2_ent = occ.copy([(2, face1)])

v1 = occ.revolve([(2, face1)], cx, cy, 0, 0, 0, 1, math.pi)
v2 = occ.revolve(face2_ent, cx, cy, 0, 0, 0, 1, -math.pi)

occ.remove([(2, face1)] + face2_ent)

# creates the domain
box = occ.addBox(0, 0, 0, dx, dy, dz)
v1_vols = [e for e in v1 if e[0] == 3]
v2_vols = [e for e in v2 if e[0] == 3]

occ.cut([(3, box)], v1_vols + v2_vols)

occ.synchronize()

xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.get_bounding_box(-1, -1)
eps = 1e-4

inlet = []
outlet = []
bottom = []
top = []
sides = []
body = []

for dim, tag in gmsh.model.get_boundary(gmsh.model.get_entities(3), oriented=False):
    s_xmin, s_ymin, s_zmin, s_xmax, s_ymax, s_zmax = gmsh.model.get_bounding_box(dim, tag)
    
    if abs(s_xmin - xmin) < eps and abs(s_xmax - xmin) < eps:
        inlet.append(tag)
    elif abs(s_xmin - xmax) < eps and abs(s_xmax - xmax) < eps:
        outlet.append(tag)
    elif abs(s_ymin - ymin) < eps and abs(s_ymax - ymin) < eps:
        sides.append(tag)
    elif abs(s_ymin - ymax) < eps and abs(s_ymax - ymax) < eps:
        sides.append(tag)
    elif abs(s_zmin - zmin) < eps and abs(s_zmax - zmin) < eps:
        bottom.append(tag)
    elif abs(s_zmin - zmax) < eps and abs(s_zmax - zmax) < eps:
        top.append(tag)
    else:
        body.append(tag)

gmsh.model.add_physical_group(2, inlet, 100, name="Inlet")
gmsh.model.add_physical_group(2, outlet, 200, name="Outlet")
gmsh.model.add_physical_group(2, bottom, 300, name="Bottom")
gmsh.model.add_physical_group(2, top, 400, name="Top")
gmsh.model.add_physical_group(2, sides, 500, name="Sides")
gmsh.model.add_physical_group(2, body, 600, name="Body")

vol_tags = [tag for dim, tag in gmsh.model.get_entities(3)]
gmsh.model.add_physical_group(3, vol_tags, 700, name="Volume")

def add_threshold_field(surfs, size_min, size_max, dist_min, dist_max):
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.set_numbers(f_dist, "FacesList", surfs)

    f_thresh = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.set_number(f_thresh, "InField", f_dist)
    gmsh.model.mesh.field.set_number(f_thresh, "SizeMin", size_min)
    gmsh.model.mesh.field.set_number(f_thresh, "SizeMax", size_max)
    gmsh.model.mesh.field.set_number(f_thresh, "DistMin", dist_min)
    gmsh.model.mesh.field.set_number(f_thresh, "DistMax", dist_max)

    return f_thresh

f_body = add_threshold_field(body, size_body, size_domain, dist_min_body, dist_max_body)
f_bnd = add_threshold_field(bottom+top+sides, size_bnd, size_domain, dist_min_bnd, dist_max_bnd)

# wake refinement
cyl_x = []
cyl_y = []
cyl_z = []

for tag in body:
    b_xmin, b_ymin, b_zmin, b_xmax, b_ymax, b_zmax = gmsh.model.get_bounding_box(2, tag)
    cyl_x.append(b_xmin)
    cyl_x.append(b_xmax)
    cyl_y.append(b_ymin)
    cyl_y.append(b_ymax)
    cyl_z.append(b_zmin)
    cyl_z.append(b_zmax)

cyl_x_center = (min(cyl_x) + max(cyl_x)) / 2.0
cyl_y_center = (min(cyl_y) + max(cyl_y)) / 2.0
cyl_diameter = max(cyl_y) - min(cyl_y)
wake_width = 2 * cyl_diameter

gmsh.model.mesh.field.add("Box", 5)
gmsh.model.mesh.field.set_number(5, "VIn", size_wake)
gmsh.model.mesh.field.set_number(5, "VOut", size_domain)
gmsh.model.mesh.field.set_number(5, "XMin", cyl_x_center)
gmsh.model.mesh.field.set_number(5, "XMax", xmax)
gmsh.model.mesh.field.set_number(5, "YMin", cyl_y_center - wake_width / 2.0)
gmsh.model.mesh.field.set_number(5, "YMax", cyl_y_center + wake_width / 2.0)
gmsh.model.mesh.field.set_number(5, "ZMin", zmin)
gmsh.model.mesh.field.set_number(5, "ZMax", max(cyl_z) * 1.5)
gmsh.model.mesh.field.set_number(5, "Thickness", cyl_diameter * 0.4)

gmsh.model.mesh.field.add("Min", 6)
gmsh.model.mesh.field.set_numbers(6, "FieldsList", [f_body, f_bnd, 5])

gmsh.model.mesh.field.set_as_background_mesh(6)

gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)

# creates quads from prisms by division
gmsh.option.set_number("Mesh.SubdivisionAlgorithm", 2)

gmsh.option.set_number("General.NumThreads", num_threads)

gmsh.model.mesh.generate(3)

gmsh.write('domain.step')
gmsh.write('mesh.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()