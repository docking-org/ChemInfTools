# Matthew Smith, 12/13/19
# Shuo Gu, 02/11/20
# Shoichet Lab: Torsion Strain Rotation Project
# In this Python script, we create a .csv file for the energy estimates and
# confidence intervals for each molecule in a .mol2 or .db2 file over a range
# of cutoff values for inclusion of a torsion pattern energy in the summation
# function. We use 8 values from 0 to 3.5 TEU. We will also include all of
# bond information for all of the torsion patterns in each molecule. The
# structure of this file is based on TL_Test.py (for using a commandline
# argument) and Test_D4_Cutoff.py (for using multiple cutoff values) Python
# scripts. This script should also be used as a model for how to use the
# functions in the implementation of this method in other code, or even in a
# Python commandline environment. The functions for this method are in the
# TL_Functions.py Python script.

# You can use this script like this in the commandline:
#   $ Python Torsion_Strain.py {.mol2 or .db2 file name}

# This script will figure out which file type you fed in :)
# For dealing with .db2 files, we use the functions in the db2_to_mol2.py
# Python script that Trent Belius wrote and I updated. Everything should now
# run in Python 3.7 and above.

# This implementation uses TL_2.1_VERSION_6.xml as the updated Torsion Library.
# Please make sure that file, TL_Functions.py, and db2_to_mol2.py are in the
# same directory as this script

# OUTPUT: This script outputs a .csv file. Each row corresponds to a different
# molecular scructure in the .mol2 or .db2 file input. The first column is the
# name of the structure, either the ZINC ID for each molecule in a .mol2 file
# or the number of the pose for a .db2 file (which should contain multiple
# poses for the same molecule). Following the name, we calculate the energy
# estimate and lower and upper bounds of the 95% confidence interval for a
# range of cutoff values for the sum method for TP_list objects. That is, by
# default, for every i in {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5} TEU, we calculate
# the energy estimate and confidence interval using i as a cutoff for inclusion
# (in the sum) of a torsion pattern's energy estimate (and CI bounds). These
# numbers should take up the next 3 * 8 = 24 columns after the name of the
# structure. For the 0-TEU-cutoff-sum, we recommend using 8 TEU as a cutoff
# for considering the entire structure to be strained. We recommend using
# 2.5 TEU as a cutoff for when using the non-zero-cutoff-sum and then
# considering any structure with a sum greater than 0 to be strained. That
# means we would be considering any structure with at least one torsion pattern
# with an estimated energy over 2.5 TEU to be strained overall. The
# non-zero-cutoff-sum is better when we want a more strict cutoff for
# considering a molecule to be strained, for example when we have fewer
# structures to test. We would use 0-TEU-cutoff-sum when we have more molecules
# to test and can afford more false negatives (true binders considered
# strained).

# OUTPUT CONTINUED: The remaining columns in the outputted .csv file contain
# information for each torsion pattern in that row's structure. For each
# torsion pattern, we present the index (starting from 0) of the torsion
# pattern, the atom indices defining the torsion pattern, the dihedral angle in
# degrees, the SMARTS pattern for the torsion rule in the Torsion Library, the
# type of torsion rule (specific or general), the type of energy calculation
# method (exact or approximate), the energy estimate and 95% confidence
# interval (lower and upper bounds) for just that torsion pattern, and a
# Boolean indicating whether or not the torsion pattern has an approximate
# method and a dihedral angle not within the second tolerance of any peak
# (whether it is "flagged"). These last three numbers for each torsion rule get
# added up for sum, depending on if the estimate is greater than the cutoff
# value for inclusion. Any torsion patterns with the approximate method do not
# have proper confidence intervals, so we report the energy estimate for each
# bound as to not mess up the sums of the bounds (which would happen if we set
# them to 0).

# EDIT: The .db2 implementation is not working yet. I could not yet get the
# db2MolSupplier function in the TL_Functions.py file to work. I tried to
# revise the original Mol2MolSupplier function to work with a string buffer
# object, which acts like a virtual file option. This function currently
# throws errors when trying to use the .fileno() method from the io module
# on a string buffer object. If you know how to fix this, please let me know!

# Shuo Gu, 01/01/20, .db2 and .db2.gz input is ready to use.
# Shuo Gu, 01/20/20, exception handling in calculation.
# Shuo Gu, 02/11/20, sorted every torsion energy from max to min

from sys import argv #For passing commandline arguments
import csv #For writing a .csv file
# We will import the other required modules when we run TL_Functions.py

script, input = argv #Input from commandline argument

exec(open("TL_Functions.py").read()) #Run that script
if input[-5:] == ".mol2": #If fed in a .mol2 file
    names, ms = Mol2MolSupplier(input)
    # A list of names and a dictionary of mol2 objects with the names as keys
    # Make sure every mol2 object in the file has a different name, or not all
    # of them will end up in the final dictionary!
    output_name = input[:len(input)-5] #Everything but the ".mol2"
elif input[-4:] == ".db2": #If fed in a .db2 file
    exec(open("db2_to_mol2.py").read()) #Run that script
    file = db2_file_like(input)
    # Creates a string buffer object called "file" that looks like a .mol2 file
    names, ms = db2MolSupplier(file)
    output_name = input[:len(input)-4] #Everything but the ".db2"
elif input[-7:] == ".db2.gz": #If fed in a .db2.gz file
    exec(open("db2_to_mol2.py").read()) #Run that script
    file = db2_file_like(input)
    # Creates a string buffer object called "file" that looks like a .mol2 file
    names, ms = db2MolSupplier(file)
    output_name = input[:len(input)-7] #Everything but the ".db2.gz"
else:
    print("Error. Please pass a .mol2 or .db2 file.")
    exit()

print(str(len(names)) + " molecules finished reading. Calculating strain energy...")

# These (names and ms) will be our list of names and dictionary of molecules.
# We will loop over every molecule in this dictionary, create a TP_list object
# for each one, and output the information we want as a line in a .csv file.
# We will repeat the sum procedure for every cutoff value that we want
file_name = output_name + "_Torsion_Strain.csv"
with open(file_name, mode = "w") as file_out:
    # Initialize the output file
    file_writer = csv.writer(file_out) #Initialize writing with csv, using
    # https://realpython.com/python-csv/ as a guide
    count = 0
    for name in names: #Loop over every molecule in the list
        mol = ms[name] #Get the molecule with that name
        if mol is not None: #Check to make sure the molecule exists
            try:
                M = TL_lookup(mol) #Create a TP_list function
                mol_est = []
                mol_est += M.sum(0)
                mol_info = M.get_TPs() #The molecule's information
                bond_info = [item for sublist in sorted(mol_info, key = lambda l:l[1], reverse=True) for item in sublist]
                # Unlist the list of lists into just a list of elements
                file_writer.writerow([name] + mol_est + bond_info)
                # Add the list as a row
            except:
                count += 1
                file_writer.writerow([name] + ["NA"])
    file_out.close()

print(str(len(names)-count) + " successful / " + str(count) + " NA. Please check: " + file_name)

# END
