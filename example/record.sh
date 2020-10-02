# Step 1.
# For each gene:
# calculate inter-group and intra-group distance
for f in `cat phy_list`; do gene=$(basename $f .phy); python3 /Users/mengguanliang/Documents/GitHub/group_genetic_distance/group_genetic_distance/group_dist.py -msa_phylip $f -group_definition species_level_group_def.txt -delimiter '\=' -b_o $gene.between-group.dist  -i_o $gene.within-group.dist ; done


# Step 2.
# distance-file list
ls | grep phy.csv > pairwise_dist_list

# Step 3.
# within-group_dist_of_multi-gene
python3 /Users/mengguanliang/Documents/GitHub/group_genetic_distance/group_genetic_distance/within-group_dist_of_multi-genes.py -pairwise_dist_list pairwise_dist_list -group_definition species_level_group_def.txt -delimiter '\=' -i_o all-genes.within-group.dist


# Step 4.
# between-group_dist_of_multi-genes
python /Users/mengguanliang/myonedrive/OneDrive/GitHub/group_genetic_distance/group_genetic_distance/between-group_dist_of_multi-genes.py -pairwise_dist_list pairwise_dist_list -group_definition species_level_group_def.txt -delimiter '\=' -i_o all-genes.between-group.dist
