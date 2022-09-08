# Matthew Smith, 12/5/19
# Shuo Gu, 02/18/20
# Shoichet Lab: Torsion Strain Rotation Project
# This Python script will contain all of the functions for using the Torsion
# Library method for estimating torsion strain energy. To use these functions,
# source this script in either another Python script or a Python commandline
# session as:
#   exec(open("TL_Functions.py").read())
# when in the same directory.

# We will use TL_2.1_VERSION_6.xml as the Torsion Library.

# Explanation of how this works from TL_Lookup_Test.py:
# Since the Torsion Library contains very specific SMARTS patterns to define
# the torsion patterns (also called torsion rules in the acutal library) and
# there are only 514 torsion patterns in TL_2.1, it is best to just search
# over all of the torsion rules in the library, instead of trying to find
# and classify every substructure in the molecule. This was Ying's suggestion.
# For each torsion rule, we will use the GetSubstructMatches function in
# the package rdkit.Chem to find the indeces of the atoms that match the
# torsion rule. If there are multiple such substructures, then this function
# should return multiple sets of indeces. We will then find the information
# we need (dihedral angle, method, estimate, lower bound, upper bound, and
# flag). Since the same torsion pattern in the molecule can have multiple
# torsion rules in the library applying to it, we need to choose the intended
# torsion rule. The authors of the original Torsion Library ordered the
# torsion rules within each subclass in decreasing specificity, so we will
# keep track of the number of the torsion rule for each matching substructure.
# We will then pick the set of information with the lowest counting index,
# assuming that this is the most specific. Since the first hierarchy class in
# the Torsion Library is the general class, we will need to loop over all of
# the specific classes first. This process is not the most efficient, since
# we are looking up and calculating information that we immediately throw away,
# but we can optimize this procedure later.

# Once we have one set of information for each set of indeces defining the
# bonds, we need to condense these by torsion pattern. We must do this because,
# for a given bond, each atom can be bonded to 2 or 3 other atoms in a torsion
# pattern. This gives between 2*2=4 and 3*3=9 possible lists of 4 atom indeces
# defining the torsion pattern, not including reversing the order of the
# indeces! Therefore we need to find all of the information sets (angle,
# estimate, etc.) with the same atom indeces defining the bond (listed
# forwards and backwards) and decide which one we will keep. Since we are
# using these energy estimates for identifying unstable conformations in
# DOCK-ing campaigns, being conservstive with our hit-picking means being
# liberal in our energy assignments. We thus keep the set of information with
# the highest energy estimate, prioritizing the estimates from specific
# hierarchy classes over ones from the general hierarchy class.

# Shuo Gu, 01/10/20, fixed a bug of input reading
# Shuo Gu, 01/12/20, Introduced interpolation between energy bins
# Shuo Gu, 02/18/20, Changed the initial energy 0 to 100 for the approximate method

import xml.etree.ElementTree as ET #For reading XML files
from rdkit import Chem #The module with the RDKit functions we will need
import os #For the function to read in the molecules from the .mol2 file
import numpy as np #For arrays and NumPy function
from math import sqrt, atan2, pi #For calculating dihedral angles
from math import ceil #Ceiling function


# Import the XML file, using this as a guide:
# https://docs.python.org/3/library/xml.etree.elementtree.html
tree = ET.parse("TL_2.1_VERSION_6.xml")
root = tree.getroot()


# We will turn the procedure for estimating torsion strain energy in
# TL_Lookup_Test.py into a function that we can call on every molecule
# in the list, which will return an object of the following class:
class TP_list(object):
    def __init__(self, indeces, angles, smarts, hc, methods, E,
    CI_l, CI_u, flags):
        self.indeces = indeces #List of lists of indeces
        self.angles = angles #List of angles
        self.smarts = smarts #List of SMARTS strings
        self.hc = hc #List of strings for hierarchy class
        self.methods = methods #List of strings of energy-estimation methods
        self.E = E #List of energy estimates
        self.CI_l = CI_l #List of lower bounds for 95% CI of energy estimates
        self.CI_u = CI_u #List of upper bounds for 95% CI of energy estimates
        self.flags = flags
        # List of Booleans flagging whether an angle whose "method" is
        # "approximate" is not observed
        self.TP_indeces = [j for j in range(len(indeces))]
        # Indeces of the torsion patterns, not the atoms!

    # Getters and setters:
    def get_indeces(self):
        return(self.indeces)
    def set_indeces(self, inds):
        self.indeces = inds
    def get_angles(self):
        return(self.angles)
    def set_angles(self, angs):
        self.angles = angs
    def get_smarts(self):
        return(self.smarts)
    def set_smarts(self, sms):
        self.smarts = sms
    def get_hc(self):
        return(self.hc)
    def set_hc(self, hcs):
        self.hc = hcs
    def get_methods(self):
        return(self.methods)
    def set_methods(self, meths):
        self.methods = meths
    def get_E(self):
        return(self.E)
    def set_E(self, Es):
        self.E = Es
    def get_CI_l(self):
        return(self.CI_l)
    def set_CI_l(self, ls):
        self.CI_l = ls
    def get_CI_u(self):
        return(self.CI_u)
    def set_CI_u(self, us):
        self.CI_u = us
    def get_flags(self):
        return(self.flags)
    def set_flags(self, fs):
        self.flags = fs
    def get_TP_indeces(self): #Don't need a setter for this
        return(self.TP_indeces)

    # A method to return the information for a subset of the torsion patterns
    def get_TPs(self, inds = None):
        # The parameter inds is a list of indeces for the torsion patterns,
        # not the atoms! We default to using all of the torsion patterns
        if inds == None:
            inds = [j for j in range(len(self.indeces))]
        # Create a list of torsion pattern info to be returned:
        tps = [] #Initialize
        for j in inds:
            tps.append([self.TP_indeces[j], self.E[j],
            self.CI_l[j], self.CI_u[j], self.indeces[j], self.angles[j],
            self.smarts[j], self.hc[j], self.methods[j], self.flags[j]])
        return(tps)

    # A method to find the sum of the energy estimates, with the 95% CI
    def sum(self, cutoff = None):
        # If at least one of the torsion patterns is flagged, then we
        # do not calculate the estimate or CI. We return -1 * (number flagged)
        flagged = sum(self.flags) #Treats True = 1 and False = 0
        if flagged > 0:
            return([-1*flagged, 0, 0]) #We do not try to calculate the CI
        else:
            if cutoff == None:
                cutoff == 0 #The default cutoff is 0, or using every angle
            ret = [0] * 3 #Initialize the returned array
            for i in range(len(self.E)):
                if self.E[i] >= cutoff:
                    # If the energy estimate is larger than the cutoff
                    ret[0] += self.E[i]
                    ret[1] += self.CI_l[i]
                    ret[2] += self.CI_u[i]
            return(ret)


# We need a function to read in .mol2 files, which RDKit does not supply to us
# We will use a function from
# https://chem-workflows.com/articles/2019/07/18/building-a-multi-molecule-mol2-reader-for-rdkit/
# Shuo has modified this function significantly to get it to work with the
# .mol2 files that this lab uses, which contain many lines of comments before
# each line "@<TRIPOS>MOLECULE"
def Mol2MolSupplier (file = None):
    names = [] #Make a list to hold the molecule names
    mols = {} #Make a dictionary
    with open(file, 'r') as f:
        fileend = os.fstat(f.fileno()).st_size
        count = 0
        line = f.readline()

        while not f.tell() == fileend:
            if line.startswith("#") or line == '\n':
                line = f.readline()
            if line.startswith("@<TRIPOS>MOLECULE"):
                count += 1
                mol = []
                mol.append(line)
                line = f.readline()
                if line != "\n" and line.split()[0].strip() not in names:
                    name = line.split()[0].strip()
                    print(name)
                else:
                    name = "mol2Number" + str(count)
                    print(name)

                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()

                    if f.tell() == fileend:
                        mol.append(line)
                        break
                block = ",".join(mol).replace(',','')
                m = Chem.rdmolfiles.MolFromMol2Block(block, sanitize=False, removeHs = False)
                names.append(name)
                mols[name] = m
    return(names, mols)

# Here is an updated version to use with the output "file" string buffer object
# created by the db2_file_like function in the db2_to_mol2.py script
def db2MolSupplier(file):
    names = [] #Make a list to hold the molecule names
    mols = {} #Make a dictionary
    with file as f: #file is already opened as a string buffer
        bufferend = len(f.getvalue())
        count = 0
        line = f.readline()
        while not f.tell() == bufferend:
            if line.startswith("#") or line == '\n':
                line = f.readline()
            if line.startswith("@<TRIPOS>MOLECULE"):
                count += 1
                mol = []
                mol.append(line)
                line = f.readline()

                name = "db2Number" + str(count)

                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()

                    if f.tell() == bufferend:
                        mol.append(line)
                        break
                block = ",".join(mol).replace(',','')
                m = Chem.rdmolfiles.MolFromMol2Block(block, sanitize=False, removeHs = False)
                names.append(name)
                mols[name] = m
    return(names, mols)


# The first function we will need to make TP_list objects is the dihedral angle
# calculator. I made this function (and its helper function) for my iPQB
# summer project of writing a script that creates Ramachandran plots
# (Make_Rama.py). I have changed some of the comments explaining the math
# First a helper function for creating unit vectors
def unit(a):
    # The argument should be a NumPy array with 1 axis
    return(a / sqrt(np.dot(a,a))) #Scales a by its norm

def dihedral(a_1, a_2, a_3, a_4):
    # The arguments should all be NumPy arrays with 1 axis and length 3
    # These atoms should be in order, with a_2 and a_3 defining the bond
    # of interest

    # The 3 displacement vectors:
    b_1 = a_2 - a_1
    b_2 = a_3 - a_2
    b_3 = a_4 - a_3

    n_1 = unit(np.cross(b_1, b_2))
    n_2 = unit(np.cross(b_2, b_3))

    # Imagine the first atom (a_1) is above the middle bond (from a_2 to a_3),
    # so that b_1 points downward. Then n_1 points out of the page

    m = unit(np.cross(n_1, b_2))
    # I moved the normalization to be after the cross product. Moving
    # the normalization should not change the end result because the cross
    # product commutes with scalar multiplication and ||n_1|| = 1

    # Looking down b_2, we can consider n_1 to be the x-axis and
    # m to be the y-axis. Then the dihedral angle is the angle that
    # n_2 makes with the x-axis when projected into this plane. Since dihedral
    # angles are measured going clockwise, we need to negate the angle
    # that we get back from atan
    x = np.dot(n_1, n_2) #Project n_2 onto n_1
    y = np.dot(m, n_2) #Project n_2 onto m
    return(-atan2(y,x) * 180 / pi) #Return the angle in degrees


# The next function we will need will calculate angular differences,
# in degrees. This function will calculate theta_1 - theta_2. See my notes
# from 9/19/19
def ang_diff(theta_1, theta_2):
    # (-180,180] -> [0, 360)
    if theta_1 < 0:
        theta_1 += 360
    if theta_2 < 0:
        theta_2 += 360
    del_theta = (theta_1 - theta_2) % 360 #Angular difference
    # [0, 360) -> (-180, 180]
    if del_theta > 180:
        del_theta -= 360
    return(del_theta)
    # Test this works: ang_diff(0, 179), ang_diff(0, -179),
    # and ang_diff(-179, 179)


# This function will allow us to do the matching for each torsion rule
def tp_match(tp, hc, j, mol, pos, bi):
    # tp is a torsion pattern, hc is the type of hierarchyClass ("general" or
    # "specific", and j is the current value for i
    # This function turned the global variables mol, positions, and bond_info
    # in TL_Lookup_Test.py into parameters mol, pos and bi, respectively
    smarts = tp.get("smarts")
    # Create the histograms for energy estimates and bounds of confidence
    # intervals, if available
    hist_E = [] #Initialize for energy estimates
    hist_l = [] #Initialize for lower bounds of CIs
    hist_u = [] #Initialize for upper bounds of CIs
    if tp.get("method") == "exact":
        for bin in tp.find("histogram_converted").findall("bin"):
            hist_E.append(float(bin.get("energy")))
            hist_l.append(float(bin.get("lower")))
            hist_u.append(float(bin.get("upper")))
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    # A list of lists
    for match in matches: #For each match
        # Some of the SMARTS for the torion patterns actually have 5 atoms.
        # We need to ingore these
        if len(match) > 4:
            continue #Go to the next match

        if mol.GetAtomWithIdx(match[0]).GetSymbol()=='H' or mol.GetAtomWithIdx(match[3]).GetSymbol()=='H':
            continue

        # First get the atom locations
        # Based on https://github.com/rdkit/rdkit/issues/1982
        r_1 = np.array(pos[match[0]])
        r_2 = np.array(pos[match[1]])
        r_3 = np.array(pos[match[2]])
        r_4 = np.array(pos[match[3]])

        theta = dihedral(r_1, r_2, r_3, r_4) #Dihedral angle

        # Changed next line from "TP.get" to "tp.get"
        if tp.get("method") == "exact": #If using the exact method
            # First figure out what bin we are in in the histogram.
            # We define the bins by the right endpoints, like we did
            # when we made the plots of the energy profiles
            bin_num = ceil(theta / 10) + 17
            # Starting with the -170 deg bin
            energy = (hist_E[bin_num]-hist_E[(bin_num+35)%36])/10.0*(theta-(bin_num-17)*10)+hist_E[bin_num]
            lower = (hist_l[bin_num]-hist_l[(bin_num+35)%36])/10.0*(theta-(bin_num-17)*10)+hist_l[bin_num]
            upper = (hist_u[bin_num]-hist_u[(bin_num+35)%36])/10.0*(theta-(bin_num-17)*10)+hist_u[bin_num]

            bi.append(
            [
            list(match), #Convert tuple to list
            theta,
            smarts,
            hc, #"general" or "specific"
            "exact",
            energy,
            lower,
            upper,
            False, #This only could apply for the approximate method
            j #We will take this out when we create the final object
            ]
            )

        else: #If using the approximate method
            not_observed = True
            # Initialize the flag for not observing that angle
            energy = 100
            # Initialize the energy, in case we cannot estimate it
            # Loop over all the possible angle peaks
            for angle in tp.find("angleList").findall("angle"):
                theta_0 = float(angle.get("theta_0")) #Peak location
                delta = ang_diff(theta, theta_0) #Angular displacement
                if abs(delta) <= float(angle.get("tolerance2")):
                    # If within the tolerance range for that peak
                    beta_1 = float(angle.get("beta_1"))
                    beta_2 = float(angle.get("beta_2"))
                    # The coefficients for the regression
                    energy = beta_1*(delta ** 2) + beta_2*(delta ** 4)
                    # Using the "not-as-small angle approximation"
                    not_observed = False
                    break
                    # Break the for loop to avoid problems if the
                    # observed angle sits at the border between two
                    # peaks
            bi.append(
            [
            list(match), #Convert tuple to list
            theta,
            smarts,
            hc, #"general" or "specific"
            "approximate",
            energy,
            energy, #Lower bound for CI, which we cannot find for approx. method
            energy, #Upper bound for CI, which we cannot find for approx. method
            not_observed,
            j #We will take this out when we create the final object
            ]
            )


# This most general function will automate what we did in TL_Lookup_Test.py
def TL_lookup(mol): #mol is read in from the .mol2 file
    positions = mol.GetConformer().GetPositions()
    # List of lists of atom coordinates. Luckily RDKit starts indexing at 0
    bond_info = []
    # Initialize an empty list that will hold the information for each bond
    i = 0 #Initialize count of torsion rules

    # Loop over all of the specific hierarchy classes
    for HC in root.findall("hierarchyClass"):
        if HC.get("name") != "GG": #Not the general class
            for TP in HC.iter("torsionRule"): #Loop over each torsion rules
                tp_match(TP, "specific", i, mol, positions, bond_info)
                i += 1 #Increase the count for the torsion rule

    # Now for the general method:
    for TP in root.find("hierarchyClass[@name='GG']").iter("torsionRule"):
        tp_match(TP, "general", i, mol, positions, bond_info)
        i += 1 #Increase the count for the torsion rule

    # Now that we have all of the torsion patterns, we need to be able to find
    # duplicates. The first such way is if the entire pattern is reversed. We
    # can fix this problem by making sure that all of the lists of indeces have
    # the second index (the first in the bond of interest) lower than the third
    # index (the second in the bond of interest)
    for  bond in bond_info: #Loop over every bond
        if bond[0][1] > bond[0][2]:
            bond[0].reverse() #Reverse this list
            bond.append(True) #Mark that we reversed this bond's indeces
        else:
            bond.append(False) #Mark that we did not reverse this bond's
            # indeces. We will remove this marking later

    # Next we condense the bond_info by the lists of 4 atoms defining the bonds.
    # We will pick the entry of bond_info that has the lowest value for i for
    # each match, since the torsion rules in the Torsion Library are arranged
    # (within each hierarchy class or hierarchy subclass) in decreasing
    # specificity, and we loop over all of the specific hierarchy classes
    # before the general one
    bond_info_red = [bond_info[0]] #Initialize a list for the reduced bond info
    # This reduced list needs at least one element for checking subelements
    for j in range(1, len(bond_info)):
        # Skip the first bond, which is already in the reduced list
        atom_0 = bond_info[j][0][0] #First atom index
        atom_1 = bond_info[j][0][1] #Second atom index
        atom_2 = bond_info[j][0][2] #Third bond index
        atom_3 = bond_info[j][0][3] #Fourth bond index

        unmatched = True #Initialize not finding a match
        for k in range(len(bond_info_red)):
            # Check against everything in the growing reduced list
            if bond_info_red[k][0][0] == atom_0 \
            and bond_info_red[k][0][1] == atom_1 \
            and bond_info_red[k][0][2] == atom_2 \
            and bond_info_red[k][0][3] == atom_3:
                # If there is a match in ALL of the atom indeces
                unmatched = False
                if bond_info[j][9] < bond_info_red[k][9]:
                    # The ninth index gives the value for i
                    # If the new bond has a lower value, then we use it
                    # to replace the current one
                    bond_info_red[k] = bond_info[j]
                    break
                    # No need to continue looking for matches, since there
                    # should be no more than 1

        if unmatched: #If no match
            bond_info_red.append(bond_info[j]) #Append the current bond

    # Now we condense the bond_info_red by the 2 atoms actually defining the
    # bond. We will pick the entry of bond_info_red for each match that has
    # the highest energy estimate, prioritizing torsion patterns from
    # "specific" hierarchy classes over ones from the "general" hierarchy class
    b_i_r = [bond_info_red[0]] #Initialize a list for the further reduced bond info
    # I used the name b_i_r over bond_info_red_2
    for j in range(1, len(bond_info_red)):
        # Skip the first bond, which is already in the reduced list
        atom_1 = bond_info_red[j][0][1] #First atom index of the bond
        atom_2 = bond_info_red[j][0][2] #Second atom index of the bond

        unmatched = True #Initialize not finding a match
        for k in range(len(b_i_r)):
            # Check against everything in the growing list
            if b_i_r[k][0][1] == atom_1 and b_i_r[k][0][2] == atom_2:
                # If there is a match in the two atom indeces defining the bond
                unmatched = False
                if bond_info_red[j][3][0] > b_i_r[k][3][0] \
                or (
                bond_info_red[j][5] > b_i_r[k][5] \
                and
                bond_info_red[j][3][0] == b_i_r[k][3][0]
                ):
                    # The third index gives "general" or "specific", and we
                    # access the first char in that string. We use:
                    #   's' > 's' gives False
                    #   'g' > 'g' gives False
                    #   'g' > 's' gives False
                    #   's' > 'g' gives True
                    # to prioritize using the specific classe over the general one.
                    # The fifth index gives the energy estimate
                    b_i_r[k] = bond_info_red[j] #Replace the current bond info
                    break
                    # No need to continue looking for matches, since there
                    # should be no more than 1

        if unmatched: #If no match
            b_i_r.append(bond_info_red[j]) #Append the current bond

    # Now that we have only the bond information that we want, we need to
    # un-do any times we reversed the atom indeces
    for  bond in b_i_r:
        if bond[10]:
            # The tenth index gives whether or not we reversed the indeces
            bond[0].reverse() #Reverse this list

    # Now that we have all the information we need, we can return the
    # object that contains it
    return(
    TP_list([bond[0] for bond in b_i_r], [bond[1] for bond in b_i_r],
    [bond[2] for bond in b_i_r], [bond[3] for bond in b_i_r],
    [bond[4] for bond in b_i_r], [bond[5] for bond in b_i_r],
    [bond[6] for bond in b_i_r], [bond[7] for bond in b_i_r],
    [bond[8] for bond in b_i_r]) #Using Python list comprehension
    )
