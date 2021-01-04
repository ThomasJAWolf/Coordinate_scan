import numpy as np
import math

############################################################################################################
## Classes and functions ###################################################################################
############################################################################################################

class Coordinate_scan(object):
    """
    Class for performing coordinate scans along internal coordinates of a molecule, like bond distances,
    angles, or dihedral angles.
    Arguments:
    startgeo      : Molecular geometry object with start coordinates of the interpolation
    endgeo        : Molecular geometry object with end coordinates of the interpolation
    startnumbering: Numpy array of numbers from zero to the number of atoms in the molecule minus 1 to 
    rearrange the geometry in case the atomic numbering in start and end geometry are not the same.
    endnumbering  : Same as start numbering, but for the end geometry
    connectivity  : list of lists defining the connectivity in the molecule in the style of a z-Matrix.
    The length of the list is equal to the number of atoms in the molecule.
    The first element of the list is an empty list. The second element is a list with one number 
    (number of atom in geometry from 0 to number of atoms in the molecule minus one) defining connectivity 
    (bond) to a next neighbor. The third element in the list is a list containing two numbers defining 
    connectivity to a next neighbor (first number, bond) and an angle between the bond and another atom
    (second number, bond angle). The fourth and all other elements in the list are lists of three numbers
    describing connectivity in terms of bond, bond angle, and dihedral angle to three more atoms.
    Created by Thomas Wolf, 01/03/2021
    """
    def __init__(self, startgeo, endgeo, startnumbering, endnumbering, connectivity):
        self.startnumbering = startnumbering
        self.endnumbering = endnumbering
        self.connectivity = connectivity
        
        # Reorder geometries:
        self.startcoord, self.startelements = self.reordergeo(startgeo.coordinates,startgeo.elements,self.startnumbering)
        self.endcoord, self.endelements = self.reordergeo(endgeo.coordinates,endgeo.elements,self.endnumbering)
        if not np.all(self.startelements==self.endelements):
            print('Warning: Numbering seems to be wrong!')
        
        # Create z-matrices from start and end geometries
        self.startDist,self.startAng,self.startDih = self.createzmat(self.startcoord,self.startelements)
        self.endDist,self.endAng,self.endDih = self.createzmat(self.endcoord,self.endelements)
    
    def reordergeo(self, geom, elements, Numbering):
        """
        Function to reorder a geometry of cartesian coordinates.
        Arguments:
        geom     : Numpy array of cartesian coordinates
        elements : List of elemental names as strings
        Numbering: Numbering for the reordered geometry
        Returns:
        geomreord: Numpy array of cartesian coordinates
        elemreord: List of elemental names as strings
        """
        geomreord = np.zeros_like(geom)
        elemreord = np.zeros_like(elements)
        for i in np.arange(len(elements)):
            geomreord[i,:] = geom[Numbering[i]]
            elemreord[i] = elements[Numbering[i]]
        return geomreord,elemreord
    
    def createzmat(self,geomreord,elemreord):
        """
        Function to generate z-matrices from cartesian coordinates.
        Arguments:
        geomreord: Numpy array of cartesian coordinates
        elemreord: List of elemental names as strings
        Returns:
        Distances: Numpy array of distances
        Angles   : Numpy array of angles
        Dihedrals: Numpy array of dihedral angles
        """
        Distances = np.zeros((len(elemreord),))
        Angles    = np.zeros((len(elemreord),))
        Dihedrals = np.zeros((len(elemreord),))


        # Generate z-matrix
        Distances[0] = np.nan
        Angles[0:2] = np.nan
        Dihedrals[0:3] = np.nan
        Distances[1] = np.sqrt(np.sum(np.square(geomreord[1]-geomreord[0])))
        Distances[2] = np.sqrt(np.sum(np.square(geomreord[2]-geomreord[1])))
        Angles[2] = np.degrees(np.arccos(np.sum((geomreord[2]-geomreord[1])*(geomreord[0]-geomreord[1]))/ \
        (np.sqrt(np.sum(np.square(geomreord[2]-geomreord[1])))* \
         np.sqrt(np.sum(np.square(geomreord[0]-geomreord[1]))))))
        for i in np.arange(3,len(elemreord)):
            Vec1 = geomreord[self.connectivity[i][0]]-geomreord[i]
            Vec2 = geomreord[self.connectivity[i][1]]-geomreord[self.connectivity[i][0]]
            Vec3 = geomreord[self.connectivity[i][2]]-geomreord[self.connectivity[i][1]]
            Norm1 = np.cross(Vec1,Vec2)
            Norm2 = np.cross(Vec2,Vec3)
            Distances[i] = np.sqrt(np.sum(np.square(Vec1)))
            Angles[i] = np.degrees(np.arccos(np.sum(-1*Vec1*Vec2)/ \
                                             (np.sqrt(np.sum(np.square(Vec1)))*np.sqrt(np.sum(np.square(Vec2))))))
            x = np.sum(np.multiply(Norm1,Norm2))
            y = np.sum(np.multiply(np.cross(Norm1,Vec2/np.sqrt(np.sum(np.square(Vec2)))),Norm2))
            Dihedrals[i] = -1*np.degrees(np.arctan2(y,x))
            if Dihedrals[i] <0:
                Dihedrals[i] = 360+Dihedrals[i]
        return Distances,Angles,Dihedrals
    
    def interpolate(self,nsteps):
        """
        Function to interpolate geometries between a start and end geometry in z-matrix format. The resulting
        geometries are returned into a list "self.xyzs" as cartesian coordinates.
        Arguments:
        nsteps: Number of interpolation steps including the start and end geometry
        """
        Distvec = self.endDist - self.startDist
        Anglvec = self.endAng - self.startAng
        Dihvec  = np.zeros_like(Anglvec)
        
        for i in np.arange(len(self.startDih)):
            Dihvec[i] = self.endDih[i] - self.startDih[i]
            if abs(Dihvec[i])>180:
                if Dihvec[i]<0:
                    Dihvec[i]=-1*(abs(Dihvec[i])-360)
                else:
                    Dihvec[i]=(abs(Dihvec[i])-360)
        
        self.xyzs = []
        for i in np.arange(nsteps):
            currdist = self.startDist+Distvec/(nsteps-1)*i
            currangl = self.startAng + Anglvec/(nsteps-1)*i
            currdih = self.startDih +Dihvec/(nsteps-1)*i
            currdih[np.where(currdih>180)[0]] = currdih[np.where(currdih>180)[0]] -360
            self.xyzs.append(self.zmat2xyz(currdist,currangl,currdih))
            
    def zmat2xyz(self,Dist,Ang,Dih):
        """
        Function to generate cartesian coordinates from a z-matrix geometry.
        Arguments:
        Dist: Array of distances
        Ang : Array of angles
        Dih : Array of dihedral angles
        Returns:
        cartc: Array of cartesian coordinates
        """
        cartc = np.zeros((len(Dist),3))
        # Origin of cartesian coordinate system is put into coordinate of first atom. 
        # x axis is the direction of the
        # second atom from the origin, y axis the aberation of the third atom from the x axis, ...
        cartc[1,:] = np.array([Dist[1],0.0,0.0])

        # vec1 is the normalized distance vector between the connectivity atoms, 
        # which define the bond angle with the new atom. 
        vec1 = cartc[self.connectivity[2][1]]-cartc[self.connectivity[2][0]]
        vec1 = vec1/np.sqrt(np.square(vec1).sum())
        # theta1 is the bond angle
        theta1 = Ang[2]/180.0*np.pi
        # vec2 is the unit vector in the direction from the bonding atom to the new atom
        vec2 = np.dot(self.rotation_matrix([0,0,1],theta1),vec1)
        # vec3 is the distance vector for the new bond.
        vec3 = vec2*Dist[2]
        cartc[2,:] = vec3+cartc[self.connectivity[2][0]]

        for i in np.arange(3,len(Dist)):
            # vec1 is the normalized distance vector between the connectivity atoms, 
            # which define the bond angle with the new atom.
            vec1 = cartc[self.connectivity[i][1]]-cartc[self.connectivity[i][0]]
            vec1 = vec1/np.sqrt(np.square(vec1).sum())
            # vec1 is the normalized distance vector between the connectivity atoms, 
            # which define the dihedral angle with the new atom.
            vec2 = cartc[self.connectivity[i][2]]-cartc[self.connectivity[i][1]]
            vec2 = vec2/np.sqrt(np.square(vec2).sum())
            theta1 = Ang[i]/180.0*np.pi
            # rotax 1 is the axis orthogonal to the three connectivity atom vectors.
            rotax1 = np.cross(vec2,-1*vec1)
            # vec3 is the unit vector for the new atom includin already the bond angle but not the dihedral.
            vec3 = np.dot(self.rotation_matrix(rotax1,theta1),vec1)
            # rotax2 is the bond axis between the atoms defining the bond angle
            rotax2 = -1*vec1
            # theta2 is the dihedral angle
            theta2 = Dih[i]/180.0*np.pi
            # vec4 is the unit vector pointing from the bonding atom in the direction of the new atom
            vec4 = np.dot(self.rotation_matrix(rotax2,theta2),vec3)
            vec5 = vec4*Dist[i]
            cartc[i,:] = vec5+cartc[self.connectivity[i][0]]
        return cartc
    
    def rotation_matrix(self,axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
        
    def writexyzs(self,items,filename):
        """
        Function to write interpolated geometries in an xyz file.
        Arguments:
        items: list of cartesian coordinates
        filename: name of the file to write coordinates into.
        """
        filename = filename + '.xyz'
        geomtransstr = []
        
        for item in items:
            xyzarr = self.xyzs[item]
            geomtransstr.append('   ' + str(len(self.startelements)) + '\n')
            geomtransstr.append(' \n')
            for i in np.arange(len(self.startelements)):
                string = self.startelements[i]
                for j in np.arange(3):
                    string += '      '
                    string += str('% 1.10f' % -xyzarr[i,j])
                string += '\n'
                geomtransstr.append(string)
        
        with open(filename,'w') as f:
            for line in geomtransstr:
                f.write(line)
