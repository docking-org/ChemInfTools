
# input:
#   data/substances.smi should have two columns: [smiles, id], no header, space separated
# output:
#   intermediate_data/substances.fp
#   the fp extesion is added by default
#   each line looks like:
#     fingerprint = 00000000|00000000|00000000|00000000|... 
python \
    ../../../teb_chemaxon_cheminf_tools/generate_chemaxon_fingerprints_py3.py \
    data/substances.smi \
    intermediate_data/substances

