# Writen by Trent Balius in the Shoichet Group
# calculates the tanamoto matrics.
# fingerprint are calculated with chemaxon
# this uses a simular chemaxon command as sea.
# bit comparisons are calculated in python


import os

from rdkit.Chem import AllChem


def str_to_bit(s):
    """this function converts a string of ones and zeros to a bit"""

    # check that the string is all zeros and ones
    if any([not (int(c) == 1 or int(c) == 0) for c in s]):
        print("error: string not zerros and ones")
        exit()

    # note the 2^7 = 10000000
    #          2^6 = 01000000
    #          2^5 = 00100000
    #          ... 
    #      1 = 2^0 = 00000001
    #            0 = 00000000

    i = len(s) - 1  # start at 7
    b_int = int(0)
    for c in s:
        x = int(c)
        if i >= 0: 
            b_int += (2**i) * x
        i -= 1

    b = bin(b_int)

    return b, b_int


def num_bit_ones(bitstring):
    """this function counts the number of ones in the bitstring"""
    return sum([1 for c in bitstring if c == '1'])


def write_smiles_strings_to_smi_file(smiles_strings_list, smi_file_path):
    with open(smi_file_path, 'w') as f:
        for smiles_string in smiles_strings_list:
            f.write(f"{smiles_string}\n")


def tanimoto(fp1, fp2):
    """this function computes the tanimoto between to fingerprints"""

    fp1_bits = fp1.split('|')
    fp2_bits = fp2.split('|')

    if len(fp1_bits) != len(fp2_bits):
        print("ERROR: bits do not agree in lenth")

    or_num_one = 0
    and_num_one = 0
    for i in range(len(fp1_bits)):
        bit1, int1 = str_to_bit(fp1_bits[i])
        bit2, int2 = str_to_bit(fp2_bits[i])
        and_bit = bin(int1 & int2)
        and_num_one = and_num_one + num_bit_ones(and_bit)
        or_bit = bin(int1 | int2)
        or_num_one = or_num_one + num_bit_ones(or_bit)

    tc = float(and_num_one) / float(or_num_one)

    return tc


def tversky_index(fp1, fp2, alpha, beta):
    """(A n B)/ [(A n B) + alpha * (A\B) + beta * (B\A)) where n is the intersection and \ is the inverse complement"""

    fp1_bits = fp1.split('|')  # A
    fp2_bits = fp2.split('|')  # B

    if len(fp1_bits) != len(fp2_bits):
        print("ERROR: bits do not agree in lenth")

    a_bs_b_num_one = 0  # inverse complement
    b_bs_a_num_one = 0  #
    and_num_one = 0  # interaction
    for i in range(len(fp1_bits)):
        bit1, int1 = str_to_bit(fp1_bits[i])
        bit2, int2 = str_to_bit(fp2_bits[i])
        and_bit = bin(int1 & int2)
        and_num_one = and_num_one + num_bit_ones(and_bit)

        b_bs_a = bin(~int1 & int2)
        b_bs_a_num_one = b_bs_a_num_one + num_bit_ones(b_bs_a)

        a_bs_b = bin(int1 & ~int2)
        a_bs_b_num_one = a_bs_b_num_one + num_bit_ones(a_bs_b)

    tv = float(and_num_one) / (float(and_num_one) + alpha * float(a_bs_b_num_one) + beta * float(b_bs_a_num_one))

    return tv


def fingerprint_vec(smiles_string_vec, radius=2, n_bits=1024, use_features=False, use_chirality=False, use_bond_types=True, include_redundant_environments=False):
    """this function call chemaxon and computes a bunch of fingerprints"""

    return [AllChem.GetMorganFingerprintAsBitVect(AllChem.MolFromSmiles(smiles_string), radius, nBits=n_bits, useFeatures=use_features, useChirality=use_chirality, useBondTypes=use_bond_types, includeRedundantEnvironments=include_redundant_environments) for smiles_string in smiles_string_vec]


def molecular_mass(smiles_string, pid):
    """this function call chemaxon and computes the molecular Mass of the molecule"""

    TMPDIR = "scratch"

    # this will get the users name. need to write temp file
    name = os.popen('whoami').readlines()[0].strip()

    generatemd = "/nfs/soft/jchem/current/bin/generatemd"

    # write smiles to file
    temp_smi_file_path = f"/{TMPDIR}/{name}/temp{pid}.smi"
    write_smiles_strings_to_smi_file([smiles_string], temp_smi_file_path)

    command = f"{generatemd} c {temp_smi_file_path} -k Mass"
    print(f"running the command: {command}")
    output = os.popen(command).readlines()
    lastline = output[len(output) - 1]
    mass = lastline.split()[1]

    # remove the temp file. 
    os.system(f"rm -fr {temp_smi_file_path}")

    return mass


def heavy_atoms(smiles_string, pid):
    TMPDIR = "scratch"

    # this will get the users name. need to write temp file
    name = os.popen('whoami').readlines()[0].strip()

    generatemd = "/nfs/soft/jchem/current/bin/generatemd"

    # write smiles to file
    temp_smi_file_path = f"/{TMPDIR}/{name}/temp{pid}.smi"
    with open(temp_smi_file_path, 'w') as f:
        f.write(f"{smiles_string}\n")

    command = f"{generatemd} c {temp_smi_file_path} -k Heavy"
    print(f"running the command: {command}")
    outlines = [line.strip() for line in os.popen(command).readlines()]
    last_line = outlines[-1]
    heavy = last_line.split()[1]
    print(outlines, last_line, heavy)

    # remove the temp file. 
    os.system(f"rm -fr {temp_smi_file_path}")

    return heavy


def get_fp(infile_path, outfile_path):
    """this function reads in smiles and writes out the footprints. it returns a footprint vector."""

    smiles_vec = []
    with open(infile_path, 'r') as f:
        for line in f.readlines():
            split_line = line.split()
            if len(split_line) > 2:
                print("ERROR:len(smiles) > 2")
                exit()
            print(split_line)
            smiles_vec.append(split_line[0])

    with open(outfile_path, 'w') as f:
        fp_vec = fingerprint_vec(smiles_vec)
        for fp in fp_vec:
            f.write(f"fingerprint = {fp.ToBitString()}\n")

    return fp_vec

