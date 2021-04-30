#!/usr/bin/env python3
#from __future__ import with_statement
from scipy import io, empty, array, recarray, dtype, linalg
import sys, os#, cPickle
from scipy import *
from math import *
import numpy as np

"""
ATOM      8  C  BASP A   2      76.283  57.361  60.309  0.50 84.80           C  
ATOM  44035  H2  TIP3 9165     -42.898 -14.686  13.233  1.00  0.00      WT8
ATOM   3949  OXT ARG E 118       1.647   7.647  53.839  1.00 68.00           O
AAAAAAIIIII AAAA RRRRCNNNN    XXXXXXXXYYYYYYYYZZZZZZZZOOOOOOBBBBBB      SSSSEECC
                    ?     ?                                             ????
12345678901234567890123456789012345678901234567890123456789012345678901234567890
0        1         2         3         4         5         6         7         8

#I don't quite follow this table, see ? above.
#note: This is in base-1 indexing, but my arrays are base-0

COLUMNS      DATA TYPE        FIELD      DEFINITION
------------------------------------------------------
 1 -  6      Record name      "ATOM    "
 7 - 11      Integer          serial     Atom serial number.
13 - 16      Atom             name       Atom name.
17           Character        altLoc     Alternate location indicator.
18 - 20      Residue name     resName    Residue name.
22           Character        chainID    Chain identifier.
23 - 26      Integer          resSeq     Residue sequence number.
27           AChar            iCode      Code for insertion of residues.
31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angs.
39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angs.
47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angs.
55 - 60      Real(6.2)        occupancy  Occupancy.
61 - 66      Real(6.2)        tempFactor Temperature factor.
77 - 78      LString(2)       element    Element symbol, right-justified.
79 - 80      LString(2)       charge     Charge on the atom.

the atom name field (4 chars) is further subdivided: 1st two chars are element
(right justified, so " C  " is carbon, "CA  " is calcium. Next char is 'distance indicator',
named by greek letters (ABGDEZH), so " CA " is an alpha carbon. Last char is branch number,
so eg " HA3" is the third hydrogen attatched to an atom. Hydrogens are named specially:
The 'distance indicator' is labelled by thich carbon the H is attached to, and if there are
many H attached to an atom, they are given names "1H", "2H" etc as the elemtn name, so eg
"1HA3" would be H #1 attached to alpha carbon, 1st branch (If that is even possible).
"""

residueCode = { "GLY": "G", "PRO": "P", "ALA": "A", "VAL": "V", "LEU": "L", 
                "ILE": "I", "MET":  "M", "CYS": "C", "PHE": "F", "TYR": "Y", 
                "TRP": "W", "HIS": "H", "LYS": "K", "ARG": "R", "GLN": "Q", 
                "ASN": "N", "GLU": "E", "ASP": "D", "SER": "S", "THR": "T" } 

pdbtype = [ ('atomId',    int),
                ('atom',  'U4'),
                ('resName',   'U4'),
                ('chain',     'U1'),
                ('resNum',    int),
                ('coords',    float,3),
                ('occupancy', float),
                ('beta',      float),
                ('segment',   'U4'),
                ('element',   'U2'),
               ('charge',    'U2'),
               ('weight',    float) ]

#to speed up loading time, if a directory 'pickedPDBs' exists the data will
#be stored in pickled form in that directory, so it can be loaded faster the next time.
def loadPDB(filename):	
	#in principle more careful error checking is needed but good enough for now
	if(os.path.exists('./pickledPDBs')):
		pickledName = os.path.join('./pickledPDBs', filename + '.pkl')
		
		if( os.path.exists(pickledName) ):
			with open(pickledName, 'rb') as f:
				pdb = cPickle.load(f)
			return pdb
		else:
			print("creating pickle file")
			pdb = loadTextPDB(filename)
			with open(pickledName, 'wb') as f:
				cPickle.dump(pdb, f, 2)
			return pdb		
	else:
		return loadTextPDB(filename)
		
	                
#loads a pdb from text. Returns a recarray, with entries as in pdbtype
def loadTextPDB(filename):
    with open(filename) as f:
    	data = f.readlines()
    
    atomLines = [line for line in data if line[0:6] == 'ATOM  ']
    return loadPDBsingle(atomLines)

def loadPDBsingle(atomLines):
    pdb = empty(len(atomLines), dtype(pdbtype))
    for n,line in enumerate(atomLines):
        pdb[n] = ( int(line[6:11]), line[12:16],  line[17:21],line[21],int(line[22:26]), ( float(line[30:38]),float(line[38:46]),float(line[46:54]) ),float(line[54:60]),		           float(line[60:67]),line[72:76].strip(),line[76:78].strip(),line[78:80].strip(),  0)
    pdb = pdb.view(recarray)
    weights = [('C', 12.01), ('H', 1.01), ('O', 16.00), ('P', 30.97), ('N', 14.01), ('Na', 22.99), ('Cl', 35.45), ('S', 28.09)]
#    for at,w in weights:
#        pdb.weight[pdb.element == at] = w
    return pdb

def loadPDBs(filename):
    with open(filename) as f:
        data = f.readlines()
    nl = len(data)
    i = 0
    atomlines = []
    pdbs = []
    for i in range(nl):
        if data[i][0:6] == 'MODEL ':
            atomlines = []
        if data[i][0:6] == 'ATOM  ':
            atomlines.append(data[i])
        if data[i][0:6] == 'ENDMDL':
            pdbs.append(loadPDBsingle(atomlines))
            print("loaded")
    return pdbs


def loadLIG(filename):
    with open(filename) as f:
        data = f.readlines()
	
        atomLines = [line for line in data if line[0:6] == 'HETATM' or line[0:4] == 'ATOM']
        pdb = empty(len(atomLines), dtype(pdbtype))

        for n,line in enumerate(atomLines):
            pdb[n] = ( int(line[6:11]),  line[12:16], line[17:21], line[21], int(line[22:26]), ( float(line[30:38]), float(line[38:46]), float(line[46:54]) ), float(line[54:60]),		           float(line[60:67]), line[72:76].strip(), line[76:78].stri(), line[78:80].strip(), 0)

        pdb = pdb.view(recarray)
        weights = [('C', 12.01), ('H', 1.01), ('O', 16.00), ('P', 30.97), 
	           ('N', 14.01), ('Na', 22.99), ('Cl', 35.45), ('S', 28.09)]
    for at,w in weights: 
        pdb.weight[pdb.element == at] = w                           
        return pdb

def writePDB(atoms, filename, renumber=False):
    f =  open(filename, 'wt')
    for a in atoms:
        f.write("ATOM  %5d %4s %3s%1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % ( a.atomId, a.atom, a.resName, a.chain, a.resNum,  a.coords[0], a.coords[1], a.coords[2], a.occupancy, a.beta, a.element, a.charge))
    f.close()

#calculate the distance between 2 atoms or two groups if provided the atom number in pdb
def distance(*args,**kwds):
    if len(args) ==3:
        pdb = args[0]
        atom1 = args[1]
        atom2 = args[2]
        if type(atom1) is int:
            for i in pdb:
                if atom1 == i[0]:
                    at1 = i[5]
                    break
        elif type(atom1) is list:
            at1 = centerOM(pdb,atom1)
        else:
            sys.exit('please provide atom number or atomnumber list')
        if type(atom2) is int:
            for i in pdb:
                if atom2 == i[0]:
                    at2 = i[5]
                    break
        elif type (atom2) is list:
            at2 = centerOM(pdb,atom2)
        else:
            sys.exit('please provide atom number or atomnumber list')
        return distanceBase(at1,at2)
    if len(args) == 5:
        pdb = args[0]
        res1 = args[1]
        res2 = args[2]
        atype1 = args[3]
        atype2 = args[4]
        if type(res1) is int:
            at1 = pdb[(pdb.resNum == int(res1)) & (pdb.atom == atype1)].coords
        elif type(res1) is list:
            atmp = [i[0] for i in pdb if i[4] in res1 and atype1 == i[1].strip() ]
            at1 = centerOM(pdb,atmp)
        else:
            sys.exit('please provide residue number or residue number list')

        if type(res2) is int:
            at2 = pdb[(pdb.resNum == int(res2)) & (pdb.atom == atype2)].coords
        elif type(res2) is list:
            atmp = [i[0] for i in pdb if i[4] in res2 and atype2 == i[1].strip() ]
            at2 = centerOM(pdb,atmp)
        else:
            sys.exit('please provide residue number or residue number list')

        return distanceBase(at1,at2)
    sys.exit('incorrect inputs')
            
def distanceBase(coord1,coord2):
    return linalg.norm(array(coord1)-array(coord2))

def centerOM(pdb,chain,atgroup):
    atoms = zeros(3)
    w = 0.
    for i in pdb[pdb.chain == chain]:
        if i[0] in atgroup:
            atoms = atoms + array(i[5])*i[11]
            w = w + i[11]
    return atoms/w

#def dihedral(*args,**kwds):
#    if len(args) ==2:
#        pdb = args[0]
#        atomnums = args[1]
#    if len(args) ==9:
#        pdb = args[0]
#        atomnums=[]
#        for i in range(4):
##            print pdb[(pdb.resNum == int(args[i*2+1])) & (pdb.atom == args[i*2+2])]
#            atomnums.append(pdb[(pdb.resNum == int(args[i*2+1])) & (pdb.atom == args[i*2+2])].atomId)
#    v = [array(pdb[pdb.atomId == atomnums[x+0]].coords[0]) -  array(pdb[pdb.atomId == atomnums[x+1]].coords[0]) for x in range(len(atomnums)-1)]
#    crossv = [cross(v[x], v[x+1]) for x in range(len(v)-1)]
#    crossv = [x/linalg.norm(x) for x in crossv]
#    m1 = cross(crossv[0],v[1]/linalg.norm(v[1]))
#    x = dot(crossv[0],crossv[1])
#    y = dot(m1,crossv[1]) 
#    return atan2(y,x)
#
#def dihedral_base(coor1,coor2,coor3,coor4):
#    v = [ array(coor1)-array(coor2),array(coor2)-array(coor3),array(coor3)-array(coor4)]
#    crossv = [cross(v[x], v[x+1]) for x in range(len(v)-1)]
#    crossv = [x/linalg.norm(x) for x in crossv]
#    m1 = cross(crossv[0],v[1]/linalg.norm(v[1]))
#    x = dot(crossv[0],crossv[1])
#    y = dot(m1,crossv[1])
#    return atan2(y,x)
#
def angle_base(a,b):
    return dot(a,b)/(linalg.norm(a)*linalg.norm(b))
def dihedral_base(p):
    # Calculate vectors between points, b1, b2, and b3 in the diagram
    b = p[:-1] - p[1:]
    # "Flip" the first vector so that eclipsing vectors have dihedral=0
    b[0] *= -1
    # Use dot product to find the components of b1 and b3 that are not
    # perpendicular to b2. Subtract those components. The resulting vectors
    # lie in parallel planes.
    v = array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    #the relationship between cos and dot product to find the desired angle.
    return degrees(arccos(v[0].dot(v[1])/(linalg.norm(v[0]) * linalg.norm(v[1]))))
                                                                                                                                                            
def dihedral_base2(p):
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2( y, x ))

def dihedral_base3(p):
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x  = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def dihedral(pdb,resnum):
    a = []
    a.append(pdb[(pdb.resNum == resnum-1) &(pdb.atom == ' C  ')].coords[0])
    a.append(pdb[(pdb.resNum == resnum) &(pdb.atom == ' N  ')].coords[0])
    a.append(pdb[(pdb.resNum == resnum) &(pdb.atom == ' CA ')].coords[0])
    a.append(pdb[(pdb.resNum == resnum) &(pdb.atom == ' C  ')].coords[0])
    a.append(pdb[(pdb.resNum == resnum+1) &(pdb.atom == ' N  ')].coords[0])
    a= array(a)
    return dihedral_base2(a[:4]),dihedral_base3(a[1:])


def myPCA(data): # the data should be
    data = array(data)
    m,n = data.shape
    mean = data.mean(axis=0)
    data = data - mean
    R = cov(data,rowvar=False)
    eval,evec = linalg.eigh(R)
    idx = argsort(eval)
    eval = eval[idx]
    evec = evec[:,idx]
    proj = dot(data,evec)
    return eval,evec,proj,mean

if __name__ == '__main__':
	loadPDB(sys.argv[1])



