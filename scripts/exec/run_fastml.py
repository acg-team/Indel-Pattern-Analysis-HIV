import subprocess as sp
import glob
from ete3 import Tree
import math
from Bio import SeqIO

def subprocess_cmd(command, tool_path):
    '''
    Runs a command through Windows Subsystem for Linux (WSL).
    
    Parameters
    -------
    command :   string
            Unix command line to be run in WSL.
    tool_path : string
            Path to the directory where the command is run.
    
    Returns
    -------
    None

    ''' 
    try:
        sp_output = sp.check_output(command, cwd=tool_path, shell=True, stderr=sp.STDOUT, universal_newlines=True)
        print(sp_output)
    except sp.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

def process_tree(input_tree):
    tree_name = input_tree.split('.')[0]
    t = Tree(input_tree)
    for n in t.traverse():
        if n.dist <= 0.000001:
            n.dist = "%.6f" % 0.000001
        else:
            n.dist = "%.6f" % n.dist
    t.write(format=5, outfile = tree_name + '_processed.nwk')

def append_sequence_to_fasta(sequence_id, source_file, output_file):
    target_sequences = list(SeqIO.parse(output_file, "fasta"))
    source_sequences = list(SeqIO.parse(source_file, "fasta"))

    found = False
    for record in source_sequences:
        if record.id == sequence_id:
            target_sequences.append(record)
            found = True
            break

    if found:
        SeqIO.write(target_sequences, output_file, "fasta")
        print(f"Sequence with ID {sequence_id} appended to {output_file}")
    else:
        print(f"Sequence with ID {sequence_id} not found in {source_file}")

def split_msa(input_tree,input_msa,split_num):

    tree_name = input_tree.split('.')[0]
    msa_name = input_msa.split('.')[0]
    print(tree_name)
    print(msa_name)
    t = Tree(input_tree)

    split_size = math.floor(len(t)/split_num)
    if split_num ==2:
        splits = [(0,split_size),(split_size,len(t))]
    if split_num ==4:
        splits = [(0,split_size),(split_size,2*split_size),(2*split_size,3*split_size),(3*split_size,len(t))]
    elif split_num == 6:
        splits = [(0,split_size),(split_size,2*split_size),(2*split_size,3*split_size),(3*split_size,4*split_size),(4*split_size,5*split_size),(5*split_size,len(t))]

    n_split = 0
    for split in splits:
        n_split += 1
        subset_leaves = []
        count = 1
        for node in t.traverse("postorder"):
            if node.is_leaf():
                if count >= split[0] and count <= split[1]:
                    subset_leaves.append(node.name)
                count += 1

            if node.name == "B.FR.1983.HXB2-LAI-IIIB-BRU.K03455":
                dist_ref = node.dist

        sequences = {}
        # Read the input FASTA file and populate the dictionary
        for record in SeqIO.parse(input_msa, 'fasta'):
            sequences[record.id] = record

        # Create and open the output FASTA file for writing
        with open(msa_name + '.' + str(n_split) + '.fasta', 'w') as output_fasta:
            # Iterate through the taxa of interest and write matching sequences
            for leave in subset_leaves:
                if leave in sequences:
                    SeqIO.write(sequences[leave], output_fasta, 'fasta')

        # prune tree
        n = Tree(input_tree)
        n.prune(subset_leaves, preserve_branch_length=True)

        # add ref sequence to each subtree and fasta subset when necessary
        if not 'B.FR.1983.HXB2-LAI-IIIB-BRU.K03455' in subset_leaves:
            root_node = n.get_tree_root()
            root_node.add_child(name="B.FR.1983.HXB2-LAI-IIIB-BRU.K03455", dist = dist_ref)

            append_sequence_to_fasta("B.FR.1983.HXB2-LAI-IIIB-BRU.K03455", \
                                    input_msa, \
                                    msa_name + '.' + str(n_split) + '.fasta')

        # set reference sequence as outgroup for each subtree
        n.set_outgroup("B.FR.1983.HXB2-LAI-IIIB-BRU.K03455")

        # write new tree into a .nwk file
        n.write(outfile = tree_name + '_' + str(n_split) + '.nwk')

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['IndelMAP', 'mafft', 'prank', 'ProPIP']

for gene in genes:
    for MSA_tool in MSA_tools:

        if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):

            CMD = ['rm -rf /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML']
            subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
            CMD = ['mkdir /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML']
            subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_*.fasta'):
                subset_num = subset.split('_')[4].split('.')[0]

                process_tree('/cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk')

                CMD = ['mkdir /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/subset_'+str(subset_num)]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

                CMD = ['perl /cfs/earth/scratch/seppemic/TM/tools/fastml/FastML.v3.11/www/fastml/FastML_Wrapper.pl ' + \
                        '--MSA_File /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta ' + \
                        '--seqType aa ' + \
                        '--outDir /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/subset_'+str(subset_num)+' ' + \
                        '--Tree /cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'_processed.nwk ' + \
                        '--SubMatrix WAG ' + \
                        '--OptimizeBL no ' + \
                        '--indelReconstruction BOTH ' + \
                        '--UseGamma no'
                        ]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/subset_'+str(subset_num))


        else:

            CMD = ['rm -rf /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML']
            subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
            CMD = ['mkdir /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML']
            subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

            split_num = 4
            split_msa('/cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment.fasta.timetree.nwk', \
                    '/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta', \
                    split_num)

            for split in range(1,split_num + 1):

                CMD = ['mkdir /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/split_'+str(split)]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

                process_tree('/cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment_'+str(split)+'.nwk')

                CMD = ['perl /cfs/earth/scratch/seppemic/TM/tools/fastml/FastML.v3.11/www/fastml/FastML_Wrapper.pl ' + \
                        '--MSA_File /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.'+str(split)+'.fasta ' + \
                        '--seqType aa ' + \
                        '--outDir /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/split_'+str(split)+' ' + \
                        '--Tree /cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment_'+str(split)+'_processed.nwk ' + \
                        '--SubMatrix WAG ' + \
                        '--OptimizeBL no ' + \
                        '--indelReconstruction BOTH ' + \
                        '--UseGamma no'
                        ]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/split_'+str(split))
