from ete3 import Tree
import subprocess as sp
from Bio import SeqIO

# gene = 'env'
# sub_boundaries = ['1:200','201:400','401:600','601:781',] 
# dist_ref = 0.124055

# gene = 'gag'
# sub_boundaries = ['1:216','217:433','434:647'] 
# dist_ref = 0.0456211

# gene = 'nef'
# sub_boundaries = ['1:243','244:485'] 
# dist_ref = 0.0636869

# gene = 'pol'
# sub_boundaries = ['1:248','249:495'] 
# dist_ref = 0.0260966

# gene = 'rev'
# sub_boundaries = ['1:244','245:488'] 
# dist_ref = 0.0773173

# gene = 'tat'
# sub_boundaries = ['1:234','235:468'] 
# dist_ref = 0.0704724

# gene = 'vif'
# sub_boundaries = ['1:217','218:435','436:651'] 
# dist_ref = 0.0240495

# gene = 'vpr'
# sub_boundaries = ['1:233','234:467','468:698'] 
# dist_ref = 0.0247406

# gene = 'vpu'
# sub_boundaries = ['1:272','273:544'] 
# dist_ref = 0.0598503

t = Tree('/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.nwk')

for boundaries in sub_boundaries:
    inf = int(boundaries.split(':')[0])
    sup = int(boundaries.split(':')[1])

    #list leaves names in subset
    subset_leaves = []
    count = 1
    for node in t.traverse("postorder"):
        if node.is_leaf():
            if count >= inf and count <= sup:
                subset_leaves.append(node.name)

            count += 1

    sequences = {}
    # Read the input FASTA file and populate the dictionary
    for record in SeqIO.parse('/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_protein_dated.fasta', 'fasta'):
        sequences[record.id] = record

    # Create and open the output FASTA file for writing
    with open('/cfs/earth/scratch/seppemic/TM/data/' + gene + '_subset/' + gene + '_sub_' + boundaries + '.fasta', 'w') as output_fasta:
        # Iterate through the taxa of interest and write matching sequences
        for leave in subset_leaves:
            if leave in sequences:
                SeqIO.write(sequences[leave], output_fasta, 'fasta')

    # prune tree
    n = Tree('/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.nwk')
    n.prune(subset_leaves, preserve_branch_length=True)

    # add ref sequence to each subtree and fasta subset when necessary
    if not 'B.FR.1983.HXB2-LAI-IIIB-BRU.K03455' in subset_leaves:
        root_node = n.get_tree_root()
        root_node.add_child(name="B.FR.1983.HXB2-LAI-IIIB-BRU.K03455", dist = dist_ref)

        CMD = ['seqkit head -n 1 /cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_protein_dated.fasta >> ' + \
                '/cfs/earth/scratch/seppemic/TM/data/' + gene + '_subset/' + gene + '_sub_' + boundaries + '.fasta']
        sp.run(CMD, shell=True)

    # set reference sequence as outgroup for each subtree
    n.set_outgroup("B.FR.1983.HXB2-LAI-IIIB-BRU.K03455")

    # write new tree into a .nwk file
    n.write(outfile='/cfs/earth/scratch/seppemic/TM/data/' + gene + '_subset/HIV-1_' + gene + '_timetree_subset_' + boundaries + '.nwk')
