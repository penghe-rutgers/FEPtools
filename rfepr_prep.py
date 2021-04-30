import argparse
from itertools import *
from sets import Set
import sys
#usage python rfepr_prep.py -p topol.top -c solv.gro -g /path/to/gmx -a atoms 1-5,7
#currently only for amber99sb-ildn in gromacs
#please use this script after pdb2gmx and before solvate and genions
#atoms related to atomtype "N*" is not supported

parser = argparse.ArgumentParser()
parser.add_argument('-p', default='topol.top', type=str, help='original topology file',required=True)
parser.add_argument('-f', type=str, help='force field parameters folder',required=True)  #the directory for FF parameters under the .../share/gromacs/top/
parser.add_argument('-c', type=str, help='1st coordinate file',required=True)
parser.add_argument('-d', type=str, help='2nd coordinate file',required=True)
parser.add_argument('-l', type=str, help='perturbed atom file',required=True)  #include the list of atoms generating the dual topology

class torsion():
    t = {}
    #the atom number in each residue type constitue the dihedral angles from phi psi to chi1-4, this would be used 
    t['ALA'] = [[-2,0,2,8],[0,2,8,10]]
    t['GLY'] = [[-2,0,2,5],[0,2,5,7]]
    t['SER'] = [[-2,0,2,9],[0,2,9,11],[0,2,4,7]]
    t['VAL'] = [[-2,0,2,14],[0,2,14,16],[0,2,4,6]]
    t['CYS'] = [[-2,0,2,9],[0,2,9,11],[0,2,4,7]]
    t['THR'] = [[-2,0,2,12],[0,2,12,14],[0,2,4,10]]
    t['LEU'] = [[-2,0,2,17],[0,2,17,19],[0,2,4,7],[2,4,7,9]]
    t['ILE'] = [[-2,0,2,17],[0,2,17,19],[0,2,4,10],[2,4,10,13]]
    t['ASN'] = [[-2,0,2,12],[0,2,12,14],[0,2,4,7],[2,4,7,9]]
    t['ASP'] = [[-2,0,2,10],[0,2,10,12],[0,2,4,7],[2,4,7,8]]
    t['ASH'] = [[-2,0,2,11],[0,2,11,13],[0,2,4,7],[2,4,7,8]]
    t['HIP'] = [[-2,0,2,16],[0,2,16,18],[0,2,4,7],[2,4,7,8]]
    t['HIS'] = [[-2,0,2,15],[0,2,15,17],[0,2,4,7],[2,4,7,8]]
    t['TRP'] = [[-2,0,2,22],[0,2,22,24],[0,2,4,7],[2,4,7,8]]
    t['PHE'] = [[-2,0,2,18],[0,2,18,20],[0,2,4,7],[2,4,7,8]]
    t['TYR'] = [[-2,0,2,19],[0,2,19,21],[0,2,4,7],[2,4,7,8]]
    t['PRO'] = [[-2,0,10,12],[0,10,12,14],[0,10,7,4],[10,7,4,1]]
    t['GLN'] = [[-2,0,2,15],[0,2,15,17],[0,2,4,7],[2,4,7,10],[4,7,10,11]]
    t['GLU'] = [[-2,0,2,13],[0,2,13,15],[0,2,4,7],[2,4,7,10],[4,7,10,11]]
    t['GLH'] = [[-2,0,2,14],[0,2,14,16],[0,2,4,7],[2,4,7,10],[4,7,10,11]]
    t['MET'] = [[-2,0,2,15],[0,2,15,17],[0,2,4,7],[2,4,7,10],[4,7,10,11]]
    t['LYS'] = [[-2,0,2,20],[0,2,20,22],[0,2,4,7],[2,4,7,10],[4,7,10,13],[7,10,13,16]]
    t['ARG'] = [[-2,0,2,22],[0,2,22,24],[0,2,4,7],[2,4,7,10],[4,7,10,13],[7,10,13,15]]
    n = ['phi','psi','chi1','chi2','chi3','chi4']

def load_top(name): #load the topol.top file to collect atomistic information
    t = {}
    para_seq = []
    with open(name) as f:
        para_name =''
        for i in f:
            if i[0] == '[':
                para_name = i.strip('[').strip('\n').strip(']').strip()  #strip the section title
                if t.get(para_name,0) == 0:
                    para_seq.append(para_name)
                    t[para_name] = []
                continue
            if i[0] == '\n' or i[0] == ';' or i[0] == '#': #skip the commented lines
                continue
            t[para_name].append(i[:-1]) #strip the last line
    return t,para_seq        

def write_top(t,seq):  #write the topol.top file after the modification, t is the topology and seq is the parameter list in sequence
    of = open('newtop.top','w')
    of.write('#include "amber99sb-ildn.ff/forcefield.itp"\n\n')
    for i in seq:
        if i == 'position_restraints':
            of.write('#include "amber99sb-ildn.ff/tip3p.itp"\n\n')
            of.write('#include "amber99sb-ildn.ff/ions.itp"\n\n')
            continue
        of.write('[ {} ]\n'.format(i))
        for j in t[i]:
            of.write('{}\n'.format(j))
    print('topology file wrote to newtop.top')
    of.close()

def load_atoms(fname):    
    with open(fname) as f:
        return [int(x) for x in f]

def modify_atom(t,list_atom,maps):  #add the atom information for the additional atoms
    Ntot = int(t['atoms'][-1].split()[0])
    for j in range(len(t['atoms'])):
        if int(t['atoms'][j].split()[0]) in list_atom:
            tmp = t['atoms'][j].split()[:8]
            at_tmp = 'dd'+tmp[1]
            tmp[1] = 'd' + tmp[1]
            tmp.append(at_tmp)
            tmp.append('0.0')
            mass_tmp = tmp[7]
            tmp.append(mass_tmp)
            t['atoms'][j] = '\t '.join(tmp)
            tmp[8] = 'd2'+at_tmp.lstrip('dd')
            tmp[9] = tmp[6]
            tmp[1] = at_tmp
            tmp[6] = '0.0'
            tmp[0] = str(maps[int(tmp[0])])
            t['atoms'].append('\t '.join(tmp))
    return t

def modify_bond(t,list_atom,maps):#adding the bond information for the additional atoms
    Ntot = int(t['atoms'][-1].split()[0])-len(list_atom)
    for j in range(len(t['bonds'])):
        tmp = [int(x) for x in t['bonds'][j].split()]
        f = 0
        if tmp[0] in list_atom:
            tmp[0] = str(maps[int(tmp[0])])
            f = 1
        if tmp[1] in list_atom:
            tmp[1] = str(maps[int(tmp[1])])
            f = 1
        if f == 1:
            t['bonds'].append('\t'.join([str(x) for x in tmp]))

    for j in range(len(t['pairs'])):
        tmp = [int(x) for x in t['pairs'][j].split()]
        f = 0
        if tmp[0] in list_atom:
            tmp[0]= str(maps[int(tmp[0])])
            f = 1
        if tmp[1] in list_atom:
            tmp[1] = str(maps[int(tmp[1])])
            f = 1
        if f == 1:
            t['pairs'].append('\t'.join([str(x) for x in tmp]))
    return t

def modify_angle(t,list_atom,maps):#adding the angle information for the additional atoms
    Ntot = int(t['atoms'][-1].split()[0])-len(list_atom)
    for j in range(len(t['angles'])):
        tmp = [int(x) for x in t['angles'][j].split()]
        f = 0
        if tmp[0] in list_atom:
            tmp[0] = str(maps[int(tmp[0])])
            f = 1
        if tmp[1] in list_atom:
            tmp[1] = str(maps[int(tmp[1])])
            f = 1
        if tmp[2] in list_atom:
            tmp[2] = str(maps[int(tmp[2])])
            f = 1
        if f == 1:
            t['angles'].append('\t'.join([str(x) for x in tmp]))
    return t                                                       

def modify_atomtypes(t,list_atom,ffnb,seq):   #the none interacting atoms needs a new atomtype here is the definition
    t['atomtypes'] = []
    uniq = []
    at = []
    for j in t['atoms']:
        if int(j.split()[0]) in list_atom and j.split()[1] not in uniq:
            tmp = ffnb[j.split()[1]]
            t['atomtypes'].append('d'+j.split()[1]+'   '+j.split()[1]+'  '+tmp)
            t['atomtypes'].append('d2'+j.split()[1]+'   '+j.split()[1]+'  '+tmp)
            t['atomtypes'].append('dd'+j.split()[1]+'   '+j.split()[1]+'  0  0  A  0  0')
            uniq.append(j.split()[1])
            at.append(j.split()[1])
    t['nonbond_params'] = []
    for i in at:
        t['nonbond_params'].append('d'+i+'   '+'d2'+i+'   1  0  0')        
    for i,j in permutations(at,2):
        t['nonbond_params'].append('d'+i+'   '+'d2'+j+'   1  0  0')
    seq.insert(0,'atomtypes')
    seq.insert(1,'nonbond_params')
    return t,seq

def loadffnb(ffnbf):   #load the nbond information into a dict
    ffnb = {}
    with open(ffnbf) as f:
        for i in f:
            if i[0] != ';' and i[0] != '[':
                ffnb[i.split()[0]] = i[16:-1]
    return ffnb

def loadffb(ffbf):  #load the bond information "dihedral" 
    ffb = []
    multi = {}
    with open(ffbf) as f:
        flag = 0
        for i in f:
            if i[:10] == '[ dihedral':
                flag = 1
            if flag == 1 and i[0] != ';' and i[0] != '#':
                ffb.append(i.split()[:8])
            if i[0] == '#':
                multi[i.split()[1]] = i.split()[2:5]
    return ffb,multi

def modify_dihedrals(t,list_atom,ffb,multi,maps,seq):  #rfepr will turn the dihedral angle off and on for any dihedral related to those atoms
    Ntot = int(t['atoms'][-1].split()[0])-len(list_atom)
    seq.insert(8,'dihedral_restraints')
    t['dihedral_restraints'] = []
    for j in range(len(t['dihedrals'])):
        dat = [x for x in t['dihedrals'][j].split()[:6]]
        for i in range(len(dat)):
            try:
                dat[i] = int(dat[i])
            except:
                pass
        f = 0
        if len(dat) == 6 and dat[5][:6] == 'torsio' and len(Set(dat[:4]) & Set(list_atom)) > 0:
            tmp = dat[5]
            dat[5] = '  '.join(multi[tmp])
            dat.append('{}   0   {} '.format(multi[tmp][0],multi[tmp][2]))
            t['dihedrals'][j] = '  '.join([str(x) for x in dat])
            for i in range(4):
                if dat[i] in list_atom:
                    dat[i]  = maps[dat[i]]
            tmp_a = dat[5]
            dat[5] = dat[6]
            dat[6] = tmp_a
            t['dihedrals'].append('  '.join([str(x) for x in dat]))
        elif len(Set(dat[:4]) & Set(list_atom)) > 0:
            types = []
            for x in dat[:4]:
                for y in t['atoms']:
                    if x == int(y.split()[0]):
                        types.append(y.split()[1].lstrip('d'))
            ff = match(types[:4],ffb,dat[4])
            for i in ff:
                tmp_dat = list(dat)
                tmp_dat.append('  '.join(i))
                tmp_dat.append('{}   0   {} '.format(i[0],i[2]))
                if f == 0:
                    t['dihedrals'][j] = '  '.join([str(x) for x in tmp_dat])
                else:
                    t['dihedrals'].append('  '.join([str(x) for x in tmp_dat]))
                f = 1
                for k in range(4):
                    if tmp_dat[k] in list_atom:
                        tmp_dat[k] = maps[tmp_dat[k]]
                tmp_a = tmp_dat[5]
                tmp_dat[5] = tmp_dat[6]
                tmp_dat[6] = tmp_a
                t['dihedrals'].append('  '.join([str(x) for x in tmp_dat]))
            if dat[4] == 4:
                tmp_dat = list(dat)
                tmp_dat[4] = '  1  '
                tmp_dat.append(' 180  0  0  180  0  1000 ')
                t['dihedral_restraints'].append('  '.join([str(x) for x in tmp_dat]))
                for k in range(4):
                    if tmp_dat[k] in list_atom:
                        tmp_dat[k] = maps[tmp_dat[k]]
                t['dihedral_restraints'].append('  '.join([str(x) for x in tmp_dat]))
    
    Ns = []
    Tor = torsion()
    for j in range(len(t['atoms'])):
        if int(t['atoms'][j].split()[0]) in list_atom and int(t['atoms'][j].split()[0]) not in Ns and t['atoms'][j].split()[4] == 'N':
            Ns.append(int(t['atoms'][j].split()[0]))
            resname = t['atoms'][j].split()[3]
            if resname == 'ASP':
                if t['atoms'][j+12].split()[2] == t['atoms'][j+11].split()[2]:
                    resname = 'ASH'
            if resname == 'HIS':
                if t['atoms'][j+17].split()[2] == t['atoms'][j+16].split()[2]:
                    resname = 'HIP'
            if resname == 'GLU':
                if t['atoms'][j+15].split()[2] == t['atoms'][j+14].split()[2]:
                    resname = 'GLH'

            for i in range(len(Tor.t[resname])):
                tmp_dat = [x + int(t['atoms'][j].split()[0]) for x in Tor.t[resname][i]]
                tmp_dat.append('  1   X  0  0  X  0  1000 ; {} {} {}'.format(resname,t['atoms'][j].split()[2],Tor.n[i]))
                t['dihedral_restraints'].append('  '.join([str(x) for x in tmp_dat]))
                for k in range(4):
                    if tmp_dat[k] in list_atom:
                        tmp_dat[k] = maps[tmp_dat[k]]
                t['dihedral_restraints'].append('  '.join([str(x) for x in tmp_dat]))
    
    return t,seq

def match(l,ffb,funct):  #to find out the dihedral angle information based on different gromacs dihedral function type
    data = []
    for lf in ffb:
        if [ x==y for x,y  in zip(l,lf[:4])].count(True) == 4 or [ x==y for x,y  in zip(l[::-1],lf[:4])].count(True) == 4:
            data.append(lf[5:8])
    if len(data) == 0:
        for lf in ffb:
            if (([ x==y for x,y in zip(l[1:3],lf[1:3])].count(True) == 2 or [ x==y for x,y in zip(l[::-1][1:3],lf[1:3])].count(True) == 2) and lf[0] == lf[3] == 'X') and funct == 9:
                data.append(lf[5:8])
    if len(data) == 0:
        for lf in ffb:
            if (([ x==y for x,y in zip(l[1:4],lf[1:4])].count(True) == 3  and lf[0] == 'X') or ([ x==y for x,y in zip(l[::-1][1:4],lf[1:4])].count(True) == 3 and lf[0]  == 'X') ) and funct == 4:
                data.append(lf[5:8])
    if len(data) == 0:
        for lf in ffb:
            if (([ x==y for x,y in zip(l[2:4],lf[2:4])].count(True) == 2  and lf[0] == lf[1] == 'X') or ([ x==y for x,y in zip(l[::-1][2:4],lf[2:4])].count(True) == 2 and lf[0] == lf[1] == 'X') ) and funct == 4:
                data.append(lf[5:8])
    if len(data)==0:
        print(l)
        sys.exit('cannot find the corresponding dihedral parameters')
    else:
        return data

def modify_coords(c1,c2,list_atom):  #add the atomic coordinates of dual topology atoms from c2 to c1
    Ntot = int(c1[0])
    Nnow = Ntot+1
    maps = {}  ##maps is a dict connecting DT atom information and ST atom information
    if Ntot != int(c2[0]): #check if two gro file has same atomic number
        sys.exit('the gro files should have the same length') 
    for i in sorted(c2[1].keys()):
        if i in list_atom:
            tmp = c2[1][i][:15]+'{:>5}'.format(int(c2[1][i][15:20])+Ntot)+c2[1][i][20:]
            c1[1][Nnow] = tmp
            maps[int(c2[1][i][15:20])] = Nnow
            Nnow += 1
    c1[0] = str(int(c1[0])+len(list_atom))
    return c1,maps

def loadgro(grof):
    c = []  #a list of 3 object, a integer for number of atoms, a dictionary for all atom coordinates with keys of atomic number, a string for box size 
    with open(grof) as f:
        fl = f.readlines()
        c.append(fl[1])
        c.append({})
        for i in range(2,len(fl)-1):
            c[1][int(fl[i][15:20])] = fl[i]
        c.append(fl[-1])
    return c

def write_gro(c):
    of = open('newgro.gro','w')
    of.write('New gro file generated by RFEPR\n')
    of.write(c[0])
    of.write('\n')
    for i in sorted(c[1].keys()):
        of.write(c[1][i])
    of.write(c[2])
    of.close()


args = parser.parse_args()
top,seq = load_top(args.p)
list_atom = []
with open(args.l) as f:
    for i in f:
        list_atom.append(int(i.strip()))
ffb,multi = loadffb(args.f+'/ffbonded.itp')
ffnb = loadffnb(args.f+'/ffnonbonded.itp')
c1 = loadgro(args.c)
c2 = loadgro(args.d)

coords,maps = modify_coords(c1,c2,list_atom)
top,seq = modify_atomtypes(top,list_atom,ffnb,seq)
top = modify_atom(top,list_atom,maps)
top = modify_bond(top,list_atom,maps)
top = modify_angle(top,list_atom,maps)
top,seq = modify_dihedrals(top,list_atom,ffb,multi,maps,seq)
write_gro(coords)
write_top(top,seq)
