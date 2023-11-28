import subprocess as sp
import ete3
from ete3 import Tree, PhyloNode
import math
import glob
import random
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
import os
import pandas as pd
import numpy as np

## read in the ARPIP files and the multiple sequence alignment
def ARPIP_to_indelMaP_format(path_to_ARPIP_ASR, path_to_alignment_file, path_to_ARPIP_nwk, path_to_ARPIP_ml_indelpoints, path_to_ouput_internal_leaves_msa, path_to_output_internal_leaves_events):
    '''
    path_to_ouput_internal_leaves_msa: define path where you want to store the combined msa of leaves and internal nodes in fasta format
    path_to_output_internal_leaves_events: define path where you want to store the combined events of leaves and internal nodes in fasta format
    '''

    ancestral_sequences = os.path.join(path_to_ARPIP_ASR)
    seq_m = os.path.join(path_to_alignment_file)
    ml_indelpoints  = os.path.join(path_to_ARPIP_ml_indelpoints)
    a_seq = AlignIO.read(ancestral_sequences,format='fasta')
    msa = AlignIO.read(seq_m,'fasta')

    n_nodes = []
    for number_node in range(len(a_seq)+len(msa)):
        n_nodes.append('node_'+str(number_node))

    tree = PhyloNode(newick=path_to_ARPIP_nwk, format=1)
    ml_indel = pd.read_csv(ml_indelpoints,sep=';',names=n_nodes, engine='python')      # warning ! engine='python' was added here



    with open(os.path.join(path_to_ouput_internal_leaves_msa),'w') as f:
        for record in a_seq:
            print('>'+record.id, file=f)
            print(record.seq,file=f)
        for record in msa:
            print('>'+record.id, file=f)
            print(record.seq,file=f)

    msa_all = SeqIO.to_dict(SeqIO.parse(os.path.join(path_to_ouput_internal_leaves_msa),'fasta'))

    for site in range(msa.get_alignment_length()):
        res = [ml_indel.iloc[site].values[x] for x in range(len(n_nodes)) if type(ml_indel.iloc[site].values[x]) == str]
        res_dict = {}
        for item in res:
            key = item.split(':')[0]
            value = item.split(':')[1]
            res_dict[key]=value
        inserted = False
        for node in tree.traverse('preorder'):
            if site == 0 and not msa_all[node.name][site]=='*':
                seq = str(msa_all[node.name].seq)
            else:
                seq = msa_all[node.name]
            new_seq = seq
            if (not node.name in res_dict.keys()) and not inserted:
                new_seq = seq[:site]+'*'+ seq[site+1:]
            elif node.name in res_dict.keys() and node.name=='root':
                inserted = True
                msa_all[node.name] = new_seq
            elif node.name in res_dict.keys() and res_dict[node.name] == 'I' and not node.name == 'root':
                new_seq = seq[:site]+seq[site].lower()+ seq[site+1:]
                inserted = True
                tmp_tree = tree.copy()
                detach_node = tmp_tree.search_nodes(name=node.name)[0]
                detach_node.detach()
                msa_all[node.name] = new_seq

                for node in tmp_tree.traverse('preorder'):
                    if site == 0 and not msa_all[node.name][site]=='*':
                        seq = str(msa_all[node.name].seq)
                    else:
                        seq = msa_all[node.name] 
                    new_seq = seq[:site]+'*'+ seq[site+1:]
                    msa_all[node.name] = new_seq
                
            msa_all[node.name] = new_seq

    with open(os.path.join(path_to_output_internal_leaves_events),'w') as f1:
        for item in msa_all.keys():
            print('>'+item,file=f1)
            print(msa_all[item],file=f1)

    print('Files converted and written to :\n')
    print(path_to_output_internal_leaves_events)

def split_msa(input_tree,input_msa,split_num):

    tree_name = input_tree.split('.')[0]
    msa_name = input_msa.split('.')[0]
    t = Tree(input_tree)

    split_size = math.floor(len(t)/split_num)
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
        n.write(outfile = tree_name + '.' + str(n_split) + '.nwk')

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

def LimitBL(path_tree):
    t = Tree(path_tree)

    for n in t.traverse():
        if n.dist <= 0.01:
            n.dist = 0.01

    t.write(outfile= path_tree + '.ModBL')

def replace_ambiguous(seq):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids
    B = 'DN'
    Z = 'EQ'
    J = 'IL'
    replaced_seq = ''
    for char in seq:
        if char == 'X':
            replaced_seq += random.choice(amino_acids)
        elif char == 'B':
            replaced_seq += random.choice(B)
        elif char == 'Z':
            replaced_seq += random.choice(Z)
        elif char == 'J':
            replaced_seq += random.choice(J)
        else:
            replaced_seq += char
    return replaced_seq

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

def process_fasta_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path + '.NoAmb.fasta', 'w') as f:
        for line in lines:
            if line.startswith('>'):
                f.write(line)
            else:
                f.write(replace_ambiguous(line))


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu'] # 'env','gag','nef','pol','rev','tat','vif','vpr','vpu'
MSA_tools = ['IndelMAP', 'mafft', 'prank', 'ProPIP'] # 'IndelMAP', 'mafft', 'prank', 'ProPIP'


for gene in genes:
    for MSA_tool in MSA_tools:

        if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):
            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_*.fasta'):
                subset_num = subset.split('_')[4].split('.')[0]

                process_fasta_file('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta')

                if gene in ['env', 'gag', 'nef']:
                    split_num = 6
                else:
                    split_num = 4

                split_msa('/cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk', \
                        '/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta.NoAmb.fasta', \
                        split_num)

                for split in range(1,split_num + 1):
                    # modify ArPIP parameter file
                    with open('/cfs/earth/scratch/seppemic/TM/scripts/exec/ArPIP_param.txt', 'r') as f:
                        lines = f.readlines()
                    with open('/cfs/earth/scratch/seppemic/TM/scripts/exec/ArPIP_param.txt', 'w') as f:
                        for line in lines:
                            if line.startswith('input.sequence.file='):
                                f.write('input.sequence.file=/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.'+str(split)+'.fasta\n')
                            elif line.startswith('input.tree.file='):
                                f.write('input.tree.file=/cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.'+str(split)+'.nwk\n')

                            elif line.startswith('output.msa.file='):
                                f.write('output.msa.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.fasta\n')
                            elif line.startswith('output.tree.file='):
                                f.write('output.tree.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.nwk\n')
                            elif line.startswith('output.ancestral.file='):
                                f.write('output.ancestral.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.ancestral.fasta\n')
                            elif line.startswith('output.node_rel.file='):
                                f.write('output.node_rel.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.node_rel.txt\n')
                            elif line.startswith('output.mlindelpoints.file='):
                                f.write('output.mlindelpoints.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.mlindelpoints.txt\n')
                            elif line.startswith('output.pipparams.file='):
                                f.write('output.pipparams.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.pipparams.txt\n')
                            else:
                                f.write(line)
                    
                    # execute command
                    CMD = ["/cfs/earth/scratch/seppemic/TM/tools/ArPIP/bpp-arpip/ARPIP params=/cfs/earth/scratch/seppemic/TM/scripts/exec/ArPIP_param.txt"]
                    subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

                    # Convert ArPIP output to a IndelMAP-like output
                    ARPIP_to_indelMaP_format('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.ancestral.fasta', \
                                            '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.fasta', \
                                            '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.nwk', \
                                            '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'.'+str(split)+'.mlindelpoints.txt', \
                                            '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'_internal_ancestral_reconstruction.'+str(split)+'.fas', \
                                            '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_'+subset_num+'_internal_evolutionary_events.'+str(split)+'.fas')

                    CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.'+str(split)+'.fasta']
                    subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
                    CMD = ['rm /cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.'+str(split)+'.nwk']
                    subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

                # remove modified alignments (without amb. char.)
                CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta.NoAmb.fasta']
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

        else:


            process_fasta_file('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta')

            if gene in ['env', 'gag', 'nef']:
                split_num = 6
            else:
                split_num = 4

            split_msa('/cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment.fasta.timetree.nwk', \
                    '/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta.NoAmb.fasta', \
                    split_num)

            for split in range(1,split_num + 1):
                with open('/cfs/earth/scratch/seppemic/TM/scripts/exec/ArPIP_param.txt', 'r') as f:
                    lines = f.readlines()
                with open('/cfs/earth/scratch/seppemic/TM/scripts/exec/ArPIP_param.txt', 'w') as f:
                    for line in lines:
                        if line.startswith('input.sequence.file='):
                            f.write('input.sequence.file=/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.'+str(split)+'.fasta\n')
                        elif line.startswith('input.tree.file='):
                            f.write('input.tree.file=/cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment.'+str(split)+'.nwk\n')

                        elif line.startswith('output.msa.file='):
                            f.write('output.msa.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.fasta\n')
                        elif line.startswith('output.tree.file='):
                            f.write('output.tree.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.nwk\n')
                        elif line.startswith('output.ancestral.file='):
                            f.write('output.ancestral.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.ancestral.fasta\n')
                        elif line.startswith('output.node_rel.file='):
                            f.write('output.node_rel.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.node_rel.txt\n')
                        elif line.startswith('output.mlindelpoints.file='):
                            f.write('output.mlindelpoints.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.mlindelpoints.txt\n')
                        elif line.startswith('output.pipparams.file='):
                            f.write('output.pipparams.file=/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.pipparams.txt\n')
                        else:
                            f.write(line)
                CMD = ["/cfs/earth/scratch/seppemic/TM/tools/ArPIP/bpp-arpip/ARPIP params=/cfs/earth/scratch/seppemic/TM/scripts/exec/ArPIP_param.txt"]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')


                # Convert ArPIP output to a IndelMAP-like output
                ARPIP_to_indelMaP_format('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.ancestral.fasta', \
                                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.fasta', \
                                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.nwk', \
                                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP.'+str(split)+'.mlindelpoints.txt', \
                                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_internal_ancestral_reconstruction.'+str(split)+'.fas', \
                                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_ArPIP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_ArPIP_internal_evolutionary_events.'+str(split)+'.fas')

                # remove split tree and MSA
                CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.'+str(split)+'.fasta']
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
                CMD = ['rm /cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment.'+str(split)+'.nwk']
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

            # remove modified alignments (without amb. char.)
            CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta.NoAmb.fasta']
            subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
