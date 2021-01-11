import numpy as np
import pyvista as pv

def addGeometry_to_view(p, geom):
    """
    Add a geometry to the current view
    """
    spheres=[]
    sphere_color=[]
    cylinders=[]
    for at in geom.atoms:
        mesh_sphere = pv.Sphere(radius=0.3, center=[at['x'], at['y'], at['z']])
        if at['label'] == 'C':
            color=[0.4, 0.4, 0.4]
        elif at['label'] == 'O':
            color=[0.9, 0.0, 0.0]
        elif at['label'] == 'H':
            color=[0.9, 0.9, 0.9]
        else:
            color=[1.0, 0.0, 0.0]
        spheres.append(mesh_sphere)
        sphere_color.append(color)
    for e in geom.getEdges():
        i, j = e
        at1 = geom.getAtom(i)
        at2 = geom.getAtom(j)
        pos1 = np.asarray([at1['x'], at1['y'], at1['z']])
        pos2 = np.asarray([at2['x'], at2['y'], at2['z']])
        vect_bond = pos2 - pos1
        middle_bond = 0.5 * (pos1 + pos2)
    
        mesh_cylinder = pv.Cylinder(center=middle_bond, direction=vect_bond, radius=.05, height=np.linalg.norm(vect_bond))
        cylinders.append(mesh_cylinder)
    
    for sphere,col in zip(spheres,sphere_color):
        p.add_mesh(sphere, color=col, show_edges=False)
    for cyl in cylinders:
        p.add_mesh(cyl, color="tan", show_edges=False)

def addVector_to_view(p, vect, color):
    """
    Add a vector to the view
    """
    v = pv.Arrow(direction=vect)
    p.add_mesh(v, color=color)
