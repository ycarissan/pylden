import numpy as np
import scipy.spatial.transform
import random
import pymatgen
import pymatgen.transformations.standard_transformations

class Geometry:
    """
    A class to handle molecular geometries.
    Geometry(filename) reads a xyz file to fill up the class.
    Geometry(filename, orient=True) reorients the geometry along the z axis
    """
    def __init__(self, filename, orient = False):
        lines = open(filename, "r").readlines()
        self.header=(lines[1])
        self.atoms = []
        self.pseudoatoms = []
        self.spherecenters = []

        for l in lines[2:]:
            a = l.split()
            try:
                lbl = a[0].strip().upper()
                position = [float(a[1]), float(a[2]), float(a[3])]
                if lbl=="BQ" or lbl=="X" or lbl=="XX":
                    print("BQ found")
                    self.spherecenters.append( { 'label': "E", 'x': position[0], 'y': position[1], 'z': position[2] } )
                else:
                    self.atoms.append( { 'label': lbl, 'x': position[0], 'y': position[1], 'z': position[2] } )
            except:
                pass

        if orient:
            filename_atoms_only = self.getgeomfilename_Atomsonly()
            self.geom_original_orientation = pymatgen.Molecule.from_file(filename_atoms_only)
            self.rotvec = get_rotation_vector_to_align_along_z(self.geom_original_orientation)
            
            orientation_rotation = scipy.spatial.transform.Rotation.from_rotvec(self.rotvec)

            for i in range(len(self.spherecenters)):
                el = self.spherecenters[i]
                position = [ el['x'], el['y'], el['z'] ]
                position = orientation_rotation.apply(position)
                self.spherecenters[i]['x'] = position[0]
                self.spherecenters[i]['y'] = position[1]
                self.spherecenters[i]['z'] = position[2]
            for i in range(len(self.atoms)):
                el = self.atoms[i]
                position = [ el['x'], el['y'], el['z'] ]
                position = orientation_rotation.apply(position)
                self.atoms[i]['x'] = position[0]
                self.atoms[i]['y'] = position[1]
                self.atoms[i]['z'] = position[2]

    def print(self):
        print("Number of atoms: {}".format(len(self.atoms)))
        for a in self.atoms:
            print(a["label"])

    def getEdges(self):
        """
        Return a list of bounded atoms
        """
        nat = len(self.atoms)
        bonds=[]
        for i in range(nat):
            for j in range(nat):
                if i!=j:
                    if self.getDist(i,j) < 1.6:
                        bonds.append([i,j])
        return bonds

    def getAtom(self, index):
        """
        Returns the ith atom
        """
        return self.atoms[index]

    def getXYZ(self, index):
        """
        Returns the x y z of the ith atom as a list
        """
        at = self.getAtom(index)
        return [ at['x'], at['y'], at['z'] ]

    def getDist(self, i, j):
        """
        Return the distance between 2 atoms
        """
        return np.linalg.norm(np.array(self.getXYZ(i)) - np.array(self.getXYZ(j)))

    def getcoords(self, atomlist):
        """ Return the position of the atoms which determine a cycle """
        coords = []
        for at in atomlist:
            pos = np.asarray(self.getXYZ(at), dtype=np.float64)
            coords.append(pos)
        return coords

    def getBarycenter(self, atomlist):
        """ Return the mass free barycenter of the atoms in the list """
        coords = self.getcoords(atomlist)
        nat = len(atomlist)
        barycenter = np.asarray([0,0,0])
        for at in coords:
            barycenter = barycenter + at/nat
        return barycenter

    def addPseudoAtom(self, coords):
        """
        Add a pseudo atom with given coordinates
        """
        self.pseudoatoms.append( { 'label': 'E', 'x': coords[0], 'y': coords[1], 'z': coords[2] } )

    def getgeomfilename_Atomsonly(self):
        """
        Returns the filename which contains the geometry with no pseudo atoms
        """
        xyztmp_filename = "tmpfile_{:05d}.xyz".format(int(random.uniform(0, 99999)))
        fio = open(xyztmp_filename, "w+")
        fio.write("{}\n\n".format(len(self.atoms)))
        for atom in self.atoms:
            fio.write("{} {} {} {}\n".format(atom['label'], atom['x'], atom['y'], atom['z']))
        fio.close()
        return xyztmp_filename

def get_angle_and_axis(op):
    """Return angle and rotation axis from an symmetry operation"""
    matQ = op.rotation_matrix
    Qxx = matQ[0, 0]
    Qyy = matQ[1, 1]
    Qzz = matQ[2, 2]
    Qzy = matQ[2, 1]
    Qyz = matQ[1, 2]
    Qxz = matQ[0, 2]
    Qzx = matQ[2, 0]
    Qyx = matQ[1, 0]
    Qxy = matQ[0, 1]
    x = Qzy-Qyz
    y = Qxz-Qzx
    z = Qyx-Qxy
    r = np.hypot(x,np.hypot(y,z))
    t = Qxx+Qyy+Qzz
    theta = np.arctan2(r,t-1)
    return theta, np.asarray([x/r, y/r, z/r])


def get_principal_axis(pga):
    theta_min = 2 * np.pi
    axis_min = np.asarray([0, 0, 1])
    for op in pga.get_symmetry_operations():
        theta, axis = get_angle_and_axis(op)
        if theta > np.pi/100 and theta < theta_min:
            theta_min = theta
            axis_min = axis
    return theta_min, axis_min

def get_rotation_vector_to_align_along_z(geom_sym):
    pga = pymatgen.symmetry.analyzer.PointGroupAnalyzer(geom_sym)
    theta, axis = get_principal_axis(pga)
#    print("Principal axis found {0[0]} {0[1]} {0[2]} angle: {1}".format(axis, theta))
    rotation_vector = np.cross([0, 0, 1], axis)
    rotation_vector = rotation_vector / np.linalg.norm(rotation_vector)
    rotation_angle = np.arcsin(np.linalg.norm(rotation_vector)/np.linalg.norm(axis))
    rotation_vector = rotation_vector * rotation_angle
#    rot = pymatgen.transformations.standard_transformations.RotationTransformation(rotation_axis, rotation_angle, angle_in_radians=True)
    return rotation_vector

def main():
    """
    Cannot work as newgeom is not defined
    """

    pga = pymatgen.symmetry.analyzer.PointGroupAnalyzer(newgeom)
    theta, axis = get_principal_axis(pga)
    print("Principal axis found {0[0]} {0[1]} {0[2]} angle: {1}".format(axis, theta))

if __name__=="__main__":
    main()
