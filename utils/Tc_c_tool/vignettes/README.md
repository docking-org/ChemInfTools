# This is a demo of doing best-first-clustering


Run scripts in the src folder from the this direcotry vignettes/best_first_clustering.

The vignette
   1) takes in a list of compounds and their smiles and their substance_id
       => data/substances.smi
   2) uses rdkit to produce ECFP fingerprints
   3) does best-first clustering using tanimoto similarity
   4) produces a list of clusters and cluster heads
       => product/cluster_detailed.list
       => product/cluster_heads.list
