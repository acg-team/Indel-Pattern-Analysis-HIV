from ete3 import PhyloNode
from Bio import AlignIO
import pandas as pd
import glob

# Correct corrupted newick tree from FastML output by comparison with the input tree
def FastML_tree(output_file_tree_FastML, input_file_tree_FastML, output_file_nodes_ancestors, output_file_ancestral_msa, output_file_new_tree):
    guide_tree = PhyloNode(newick = output_file_tree_FastML, format=1)
    original_tree = PhyloNode(newick = input_file_tree_FastML, format=1) #this one is the real original tree taken as input by FastML
    n_rel= pd.read_csv(output_file_nodes_ancestors, sep=r'\s+', skiprows=2,names=['name','parent','child1','child2'],engine='python')
    msa = AlignIO.read(output_file_ancestral_msa,"fasta")
    node_names_msa = [record.id for record in msa._records]
    node_names_tree = [node.name for node in guide_tree.traverse('preorder')]
    node_dictionary = {}

    for node in guide_tree.traverse("preorder"):    
        if not node.name in node_names_msa:

            # get ID of parent node in the tree
            parent = node.up.name

            # get ID of both child of this parent in the relation file
            possible_1 = n_rel[n_rel["parent"] == parent]["name"].iloc[0]
            possible_2 = n_rel[n_rel["parent"] == parent]["name"].iloc[1]

            # check if parent's child ID in node_name_tree
            if node.name == possible_1.replace('-', '_') or (possible_1 not in node_names_tree and possible_2 in node_names_tree):
                node_dictionary[node.name]=possible_1
            elif node.name == possible_2.replace('-', '_') or (possible_1 in node_names_tree and possible_2 not in node_names_tree):
                node_dictionary[node.name]=possible_2


            elif possible_1 not in node_names_tree and possible_2 not in node_names_tree:
                print("NOT AN UNIQUE DECISION NEED THE RIGHT INPUT TREE")

                corrupted_leaves = []
                for node in guide_tree.traverse("preorder"):
                    if node.is_leaf():
                        corrupted_leaves.append(node.name)

                original_leaves = []
                for node in original_tree.traverse("preorder"):
                    if node.is_leaf():
                        original_leaves.append(node.name)

                for index in range(len(corrupted_leaves)):
                    if corrupted_leaves[index] != original_leaves[index]:
                        node_dictionary[corrupted_leaves[index]] = original_leaves[index]
                        print(f"Replaced corrupted node : {corrupted_leaves[index]} by original : {original_leaves[index]}")



    for node in guide_tree.traverse('preorder'):
        if node.name in node_dictionary.keys():
            node.name = node_dictionary[node.name] 

    with open(output_file_new_tree, 'w') as f:
        print(guide_tree.write(format=1).split(';')[0]+"N1;", file=f)

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu'] # 'env','gag','nef','pol','rev','tat','vif','vpr','vpu'
MSA_tools = ['IndelMAP','mafft','prank','ProPIP'] # 'IndelMAP','mafft','prank','ProPIP'

for gene in genes:
    for MSA_tool in MSA_tools:
        
        if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):
            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/subset_*'):
                output_file_tree_FastML = subset + '/tree.newick.txt'
                output_file_nodes_ancestors = subset + '/tree.ancestor.txt'
                output_file_ancestral_msa = subset + '/seq.marginal_IndelAndChars.txt'

                output_file_new_tree = subset + '/tree.rewritten.newick.txt'
                subset_num = subset.split('_')[2]
                input_file_tree_FastML = '/cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'_processed.nwk'

                FastML_tree(output_file_tree_FastML, input_file_tree_FastML, output_file_nodes_ancestors, output_file_ancestral_msa, output_file_new_tree)

        else:
            for split in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/split_*'):
                split_num = split.split('_')[2]

                output_file_tree_FastML = split + '/tree.newick.txt'
                output_file_nodes_ancestors = split + '/tree.ancestor.txt'
                output_file_ancestral_msa = split + '/seq.marginal_IndelAndChars.txt'

                output_file_new_tree = split + '/tree.rewritten.newick.txt'

                input_file_tree_FastML = '/cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment_'+str(split_num)+'_processed.nwk'

                FastML_tree(output_file_tree_FastML, input_file_tree_FastML, output_file_nodes_ancestors, output_file_ancestral_msa, output_file_new_tree)


