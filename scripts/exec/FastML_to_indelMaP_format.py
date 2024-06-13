import os
from ete3 import PhyloNode,Tree
from Bio import AlignIO,SeqIO
import glob
import subprocess as sp
import pandas as pd

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

def FastML_to_indelmap_format(ASR,tree,output_msa,output_events,gene,tool):

    ancestral_sequences = ASR
    # load tree
    t = PhyloNode(tree, format=1)

    msa = AlignIO.read(ancestral_sequences,'fasta')
    msa_all = SeqIO.to_dict(SeqIO.parse(os.path.join(ancestral_sequences),'fasta'))

    deletion = 0
    insertion = 0
    mutation = 0
    reinsertion = 0 # amino acid that was inserted again after a deletion (same column in MSA)

    for site in range(msa.get_alignment_length()):
        for node in t.traverse('preorder'):

            if site == 0:
                seq = str(msa_all[node.name].seq)
            else:
                seq = msa_all[node.name]

            if node.is_root() and msa_all[node.name][site] == '-':
                new_seq = seq[:site]+'*'+ seq[site+1:]

            elif (not node.is_root()) and (msa_all[node.name][site] != msa_all[node.get_ancestors()[0].name][site]):
                if (msa_all[node.get_ancestors()[0].name][site] == '*') and (msa_all[node.name][site] == '-'):
                    new_seq = seq[:site]+'*'+ seq[site+1:]
                elif msa_all[node.get_ancestors()[0].name][site] == '*' and msa_all[node.name][site] != '-':
                    insertion += 1
                    new_seq = seq[:site] + msa_all[node.name][site].lower() + seq[site+1:]

                elif msa_all[node.name][site] == '-' and msa_all[node.get_ancestors()[0].name][site] != '*' and msa_all[node.get_ancestors()[0].name][site] != '-':
                    deletion += 1
                    new_seq = seq
                elif msa_all[node.get_ancestors()[0].name][site] == '-' and msa_all[node.name][site] != '-' and msa_all[node.name][site] != '*':
                    reinsertion += 1
                    new_seq = seq[:site] + msa_all[node.name][site].lower() + seq[site+1:]
                else:
                    mutation += 1
                    new_seq = seq
            else:
                new_seq = seq

            msa_all[node.name] = new_seq

    print('Conversion done for '+gene+' gene, alignement using '+tool+':')
    print('total of mutations : ' + str(mutation))
    print('total of deletions : ' + str(deletion))
    print('total of insertions : ' + str(insertion))
    print('total of re-insertions : ' + str(reinsertion)+'\n\n\n')


    # print MSA with internal leaves events in new file
    with open(os.path.join(output_events),'w') as f1:
        for item in msa_all.keys():
            print('>'+item,file=f1)
            print(msa_all[item],file=f1)

    return(gene, tool, 'FastML', mutation, deletion, insertion, reinsertion)


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu'] # 'env','gag','nef','pol','rev','tat','vif','vpr','vpu'
MSA_tools = ['IndelMAP', 'mafft', 'prank', 'ProPIP'] # 'IndelMAP', 'mafft', 'prank', 'ProPIP'

         
Indel_count_df = pd.DataFrame(columns=['gene', 'MSA', 'ASR', 'substitution', 'deletion', 'insertion', 'reinsertion'])

for gene in genes:
    for MSA_tool in MSA_tools:

        if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):
            
            indel_list = [gene, MSA_tool, 'FastML', 0,0,0,0]
            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/subset_*'):

                ASR = subset + '/seq.marginal_IndelAndChars.txt'
                tree = subset + '/tree.rewritten.newick.txt'

                output_msa = subset + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_FastML_MSA_all.fasta'
                output_events = subset + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_FastML_all_events.fasta'


                indel_count = FastML_to_indelmap_format(ASR,tree,output_msa,output_events,gene,MSA_tool)
                for i in [3,4,5,6]:
                    indel_list[i] = indel_list[i] + indel_count[i]
            
            Indel_count_df.loc[len(Indel_count_df)] = indel_list
            

        else:
            indel_list = [gene, MSA_tool, 'FastML', 0,0,0,0]
            for split in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_FastML/split_*'):

                ASR = split + '/seq.marginal_IndelAndChars.txt'
                tree = split + '/tree.rewritten.newick.txt'

                output_msa = split + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_FastML_MSA_all.fasta'
                output_events = split + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_FastML_all_events.fasta'


                indel_count = FastML_to_indelmap_format(ASR,tree,output_msa,output_events,gene,MSA_tool)

                for i in [3,4,5,6]:
                    indel_list[i] = indel_list[i] + indel_count[i]

            Indel_count_df.loc[len(Indel_count_df)] = indel_list


output_file = '/cfs/earth/scratch/seppemic/TM/results/others/FastML_indel_count.csv'
with open(output_file, 'w') as f:
    Indel_count_df.to_csv(f, sep='\t')