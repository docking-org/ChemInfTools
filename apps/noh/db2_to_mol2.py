#!/user/bin/python

#################################################################################################################
##
## This libary for reading in db2 file
##
#################################################################################################################
## Writen by Trent Balius in the Shoichet Lab, UCSF in 2015
#################################################################################################################

# Modified by Matthew Smith on 12/13/19
# I updated this script to run in Python 3.7 using the 2to3 package (welcome
# to the future), commented out all the print statements to save your Terminal
# window, and changed the write_mol2 commands to just return mol2 objects.
# Also, importantly, I copied and renamed the Mol, atom, and bond classes from
# Trent's mol2.py script, which was written in Python 2 and is seriously
# broken after I passed it through 2to3. Maybe this is better, because now we
# don't need to keep track of all the functions and classes in that big script.

# The biggest change was converting the main() function into
# db2_file_like(file_name). Instead of taking a parameter from the commandline
# and writing a .mol2 file, we will pass a file_name as a function argument
# and return a a string buffer (or memory file). This object acts like a
# standard file object but only lives in memory, without creating a new file in
# the current directory. All of the normal methods for file objects will work
# on this string buffer object, so we can call the object "file" like we did
# when we had a real file object. I took this fix from my Make_Rama.py script
# for making Ramachandran plots from .pdb files, which was my iPQB summer
# project. I also took the important parts of the write_mol2 and append_mol2
# functions from mol2.py to get this to work properly.

# We will feed these functions into our method based on the Torsion Library
# after sourcing this script in whatever environment we are working in

# Modified by Shuo Gu on 01/01/20, .db2 and .db2.gz input is ready to use.

import math, sys
import os.path
import cmath
from math import sqrt
import gzip
from io import StringIO #For creating a file-like object as a string buffer

# Don't need mol2
# import mol2
#import mol2_debug as mol2

# Rescued from the now-broken mol2.py file (RIP Python 2):
class mol2_Mol:
    def __init__(self,header,name,atom_list,bond_list,residue_list):
        self.header = str(header)
        self.name = str(name)
        self.atom_list = atom_list
        self.bond_list = bond_list
        self.residue_list = residue_list
class mol2_atom: #I changed the name from just atom
    def __init__(self,X,Y,Z,Q,type,name,num,resnum,resname):
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
        self.Q = float(Q)
        self.heavy_atom = False
        self.type = type
        self.name = name
        self.num  = int(num)
        self.resnum  = int(resnum)
        self.resname = resname
class mol2_bond: #I changed the name from just bond
     def __init__(self,a1_num,a2_num,num,type):
        self.a1_num = int(a1_num)
        self.a2_num = int(a2_num)
        self.num = int(num)
        self.type = type



#################################################################################################################
#################################################################################################################
# data structure to store information about each residue with the docked ligand.
class db2Mol:
    def __init__(self,header,atom_list,bond_list,coord_list,seg_list,conformer_list):
        self.header         = str(header)
        #self.name           = str(name)
        self.atom_list      = atom_list
        self.bond_list      = bond_list
        self.coord_list     = coord_list
        self.seg_list       = seg_list
        self.conformer_list = conformer_list

class db2atom: # A
    def __init__(self,Q,type,name,num):
        #self.X = float(X)
        #self.Y = float(Y)
        #self.Z = float(Z)
        self.Q = float(Q)
        self.heavy_atom = False
        self.type = type
        self.name = name
        self.num  = int(num)
	#self.resnum  = int(resnum)
	#self.resname = resname
class db2bond: # B
     def __init__(self,a1_num,a2_num,num,type):
        self.a1_num = int(a1_num)
        self.a2_num = int(a2_num)
        self.num = int(num)
        self.type = type
class db2coords: # X
    def __init__(self,num,atomnum,segnum,X,Y,Z):
        self.num = int(num)
        self.atomnum = int(atomnum)
        self.segnum = int(segnum)
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
class db2segment: #
    def __init__(self,num,start,stop):
        self.num = int(num)
        self.start = int(start)
        self.stop = int(stop)
class db2conformer: # C
    def __init__(self,num,seglist):
        self.num = int(num)
        self.seglist = seglist



#################################################################################################################
#################################################################################################################
#def read_Mol2_filehandel(filehandel,startline):
#    lines  =  filehandel.readlines()
#def read_Moldb2_lines(lines,startline):
def read_Moldb2_file(file):
    # reads in data from multi-Mol2 file.

#T ## namexxxx (implicitly assumed to be the standard 7)
#M zincname protname #atoms #bonds #xyz #confs #sets #rigid #Mlines #clusters
#M charge polar_solv apolar_solv total_solv surface_area
#M smiles
#M longname
#[M arbitrary information preserved for writing out]
#A stuff about each atom, 1 per line
#B stuff about each bond, 1 per line
#X coordnum atomnum confnum x y z
#R rigidnum color x y z
#C confnum coordstart coordend
#S setnum #lines #confs_total broken hydrogens omega_energy
#S setnum linenum #confs confs [until full column]
#D clusternum setstart setend matchstart matchend #additionalmatching
#D matchnum color x y z
#E


    # reads in data from multi-Mol2 file.

    if input[-4:] == ".db2":
        file1 = open(file, 'r')
    if input[-7:] == ".db2.gz":
        file1 = gzip.open(file, 'r')

    lines = file1.readlines()

    file1.close()

    mol_list = []

    header = ''
    for line in lines:
         if input[-7:] == ".db2.gz":
             line = line.decode('utf-8')

         linesplit = line.split() #split on white space

         if(line[0] == "M"): # ATOM
             header = header + line[1:-1]
             atomlist  = []
             bondlist  = []
             coordlist = []
             seglist   = []
             conflist  = []

         elif(line[0] == "A"): # ATOM
            #print line
            #print line[0]
            atomnum    = linesplit[1]
            atomname   = linesplit[2]
            atomtype   = linesplit[3]
            atomcharge = linesplit[6]
            tempatom = db2atom(atomcharge,atomtype,atomname,atomnum)
            atomlist.append(tempatom)
            #print atomnum, atomname, atomtype, "q = ", atomcharge

         elif(line[0] == "B"): # BOND
            #print line
            #print line[0]
            bondnum  = linesplit[1]
            atom1 = linesplit[2]
            atom2 = linesplit[3]
            bondtype = linesplit[4]
            #print "atom1,atom2,bondnum,bondtype = [" , atom1,atom2,bondnum,bondtype, "]"
            tempbond = db2bond(atom1,atom2,bondnum,bondtype)
            bondlist.append(tempbond)
            #print bondnum, atom1,atom2, bondtype
         elif(line[0] == "X"): # COORDS
            #exit()
            #print line
            #print line[0]
            coordnum = linesplit[1]
            atomnum  = linesplit[2]
            segnum   = linesplit[3]
            X        = linesplit[4]
            Y        = linesplit[5]
            Z        = linesplit[6]
            temp_coord = db2coords(coordnum,atomnum,segnum,X,Y,Z)
            coordlist.append(temp_coord)
            #print coordnum,X,Y,Z
         #elif(line[0] == "R"): # Rigid
         #   print line
         elif(line[0] == "C"): # Segment
            #print line
            #print line[0]
            confnum    = linesplit[1]
            coordstart = linesplit[2]
            coordstop  = linesplit[3]
            #print confnum, coordstart, coordstop
            tempseg = db2segment(confnum, coordstart, coordstop)
            seglist.append(tempseg)
            numold = 1
            fristflag = True
         elif(line[0] == "S"): # set -- Conformer
            #print line
            num = int(linesplit[1])
            num2 = int(linesplit[2])
            #print numold, num
            if (fristflag):
                fristflag = False
                segnum_list = []
            elif (numold != num): # we know when it is a new conf when this number changes.
                #print "new conformation"
                tempconf = db2conformer(num,segnum_list)
                conflist.append(tempconf)
                segnum_list = []
                # This fist line does not contain the segment information
                # The second, and higher lines have more information
            else: # there may be multiple lines for enumerating sagments for one conformer.
                #print "continueing, size of segnum_list = " + str(len(segnum_list))
                numofseg = linesplit[3]
                #print numofseg, len(linesplit)
                for i in range(4,len(linesplit)):
                    segnum_list.append(int(linesplit[i]))
            numold = num
         elif(line[0] == "E"): # ATOM
             #if (len(segnum_list) > 0): # this is to put the last conformation in the the list
             tempconf = db2conformer(num,segnum_list)
             conflist.append(tempconf)

             # print(("atomnum =", len(atomlist)))
             # print(("bondnum =", len(bondlist)))
             # print(("coordnum =", len(coordlist)))
             # print(("segnum =", len(seglist)))
             # print(("confnum =", len(conflist)))
             tempmol = db2Mol(header, atomlist, bondlist, coordlist, seglist, conflist)  # this is an ensomble of conformation
             header = ''
             mol_list.append(tempmol)
             #exit()
         else:
             # print(("Warrning: " + line[0] + " is not found in the if statments. "))
             #exit()
             pass

    return mol_list


#################################################################################################################
#################################################################################################################

def convert_db2_to_mol2(db2mols):

    allmol2s = []
    # loop over each molecule
    for mol in db2mols:
         # loop over each conformer in the molcule
         mol2mols = []
         for conf in mol.conformer_list:
              #print conf.seglist
              # the conformer is defined as a set of segement of the molecule
              mol2atomlist = []
              residue_list = {}
              for segint in conf.seglist:
                  segment =  mol.seg_list[segint-1]
                  # print((segment.num, segment.start, segment.stop))
                  # the segement point to a bunch of coordenates, we know what atom the coordenate coresponds to.
                  #print segment.start,segment.stop, range(segment.start,segment.stop)
                  for coordint in range(segment.start,segment.stop+1):
                      coord = mol.coord_list[coordint-1]
                      # print((coord.num, coord.atomnum, coord.segnum,coord.X,coord.Y,coord.Z))
                      tempatom = mol.atom_list[coord.atomnum-1]
                      #X,Y,Z,Q,type,name,num,resnum,resname):
                      res_num = 1
                      resname = "lig"
                      mol2atom = mol2_atom(coord.X,coord.Y,coord.Z,tempatom.Q,tempatom.type,tempatom.name,tempatom.num,res_num,resname)
                      if res_num in residue_list:
                         residue_list[res_num].append(mol2atom)
                      else:
                         residue_list[res_num] = [mol2atom]
                      #residue_list[res_num] = [tempatom]
                      mol2atomlist.append(mol2atom)
              mol2bondlist = []
              for bond in mol.bond_list:
                  #,a1_num,a2_num,num,type
                  mol2bond = mol2_bond(bond.a1_num,bond.a2_num,bond.num,bond.type)
                  mol2bondlist.append(mol2bond)
              mol2mol = mol2_Mol("","",mol2atomlist,mol2bondlist,residue_list)
              mol2mols.append(mol2mol)
         allmol2s.append(mol2mols)
         #exit()
         #return allmol2s

    return allmol2s
#################################################################################################################
#################################################################################################################

# I changed the main() function into a function that creates a string buffer
# object. We can then feed this object into the Mol2MolSupplier function in
# TL_Functions.py. I needed to use large bits of the write_mol2 and append_mol2
# functions from the mol2.py script to get this to work.
def db2_file_like(file_name):
    db2mols = read_Moldb2_file(file_name)
    allmol2s = convert_db2_to_mol2(db2mols)
    allmol2s_simple = [molecule for mol2mols in allmol2s for molecule in mol2mols]
    # My hacky list comprehension way to unlist this. Maybe just change the
    # convert_db2_to_mol2 function?

    file = StringIO() #Initialize the string buffer
    pose_no = 0
    # Initialize a pose number counter, which will take the place of a molecule
    # name
    for molecule in allmol2s_simple:
        # Below is the crap from the write_mol2 and append_mol2 functions. I
        # needed to change the initializing and appending an outmol2 file
        # object to just appending the "file" string buffer object above. I
        # also needed to fix a lot of indentation errors...

        # define a dictionary for help renumbering
        atom_dic = {}
        resid_dic = {}
        count = 1
        for atom in molecule.atom_list:
            # if not atom_dic.has_key(atom.num):
            # I commented this out, as has_key was removed in Python 3
            if atom.num not in atom_dic:
               atom_dic[atom.num] = count
               #print atom.num, ",", count,",", atom_dic[atom.num]
               count=count+1
        count = 1
        for resnum in molecule.residue_list.keys():
            resid_dic[resnum] = count
            count=count+1

        #print molecule.header
        #file.write("##########                 Name:    " + str(pose_no))
        # This line will allow the db2MolSupplier function in
        # TL_Functions.py to understand that this is the beginning of a
        # molecule
        pose_no += 1 #Increase the count
        file.write(molecule.header)      #dock info after #s
        file.write("@<TRIPOS>MOLECULE\n")      #start the MOLECULE RTI (Record Type Indicator)
        file.write(molecule.name+'\n')         #print MOL2FILE name of the molecule
        file.write("%-5d %-5d %-5d 0     0\n" % (len(molecule.atom_list),
        len(molecule.bond_list), len(molecule.residue_list.keys())))
        # For now, the number of residues is hard-coded to 1. To be fixed.
        file.write("SMALL\n")                  #mol_type
        file.write("USER_CHARGES\n")           #charge_type

        file.write("@<TRIPOS>ATOM\n")      #start the ATOM RTI (Record Type Indicator)
        for j in range(0,len(molecule.atom_list)):
            #print atom_dic[molecule.atom_list[j].num], molecule.atom_list[j].num
            file.write("%-6d %-4s %9.4f %9.4f %9.4f %-5s %4s %6s %9.4f\n" %
            (atom_dic[molecule.atom_list[j].num], molecule.atom_list[j].name, molecule.atom_list[j].X, molecule.atom_list[j].Y,
            molecule.atom_list[j].Z, molecule.atom_list[j].type, resid_dic[molecule.atom_list[j].resnum],
            molecule.atom_list[j].resname, molecule.atom_list[j].Q))

        file.write("@<TRIPOS>BOND\n")
        count = 1
        for m in range(0,len(molecule.bond_list)):
            file.write("%-5d %-5d %-5d %s\n" % (count,
            atom_dic[molecule.bond_list[m].a1_num], atom_dic[molecule.bond_list[m].a2_num], molecule.bond_list[m].type))
            count = count + 1

        file.write("@<TRIPOS>SUBSTRUCTURE\n")
        count = 1
        for resnum in molecule.residue_list.keys():
                #outmol2.write("%-3d %-5s %-5d RESIDUE    1   A     %-5s 1\n" % (resnum,
                file.write("%-3d %-5s %-5d RESIDUE    1   A     %-5s 1\n" % (resid_dic[resnum],
                molecule.residue_list[resnum][0].resname, # residue name
                atom_dic[molecule.residue_list[resnum][0].num], molecule.residue_list[resnum][0].resname[0:3]))   # atom num of first atom in this residue
                count = count + 1
    file.seek(0) #Return the "reading line" to the beginning of the "file"
    return(file)
