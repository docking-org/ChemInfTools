
# Input:
#    (1) intermediate_data/substances.fp
#    (2) data/substances.smi 
#    (3) 0.5 
#    (4) 2000
#   figureprint file = intermediate_data/substances.fp
#   smiles file = data/substances.smi
#   tc threshold = 0.500000
#   maxium clusters to be generated = 2000
#   number of lines in file: 22
# Output:
#   cluster_details.list
#      columns: [cluster_head_index, substance_index, substance_name, substance_to_cluster_head_tc]
#      all compounds in in each cluster and the tanimoto coefficient ot the cluster head
#   cluster_heads.list
#      columns: [cluster_head_index, substance_index, substance_name, substance_to_cluster_head_tc]
#      just the single head molecule for each cluster
../../best_first_clustering/best_first_clustering \
    intermediate_data/substances.fp \
    data/substances.smi \
    0.5 \
    2000

mv cluster_details.list product/
mv cluster_head.list product/
