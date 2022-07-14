# writen by Trent Balius in 2016

import os

import fire

import tanimoto_tversky_cal_axon_lib_py3 as tc_calc


def main(smiles_file_1, outfile_prefix):
    outfile_f = f"{outfile_prefix}.fp"

    # read in names from smiles file.
    names = []
    with open(smiles_file_1, 'r') as f:
        for line in f.readlines():
            print(line)
            name = line.split()[1]
            names.append(name)

    # if the fp file already exists, just read in the fingerprints.
    # otherwise calculate the fingerprints.
    if os.path.isfile(outfile_f):
        print(f"{outfile_f} exists.")
        exit()
    else:
        fp_vec_1 = tc_calc.get_fp(smiles_file_1, outfile_f)


if __name__ == '__main__':
    fire.Fire(main)


