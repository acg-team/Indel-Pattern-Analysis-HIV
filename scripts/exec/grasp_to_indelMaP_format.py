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

def grasp_to_indelmap_format(ASR,leaves,tree,output_msa,output_events,gene,tool):

    ancestral_sequences = ASR
    leaves_msa = leaves
    path_to_ouput_internal_leaves_msa = output_msa
    path_to_output_internal_leaves_events = output_events
    # load tree
    t = PhyloNode(tree, format=1)



    a_seq = AlignIO.read(ancestral_sequences,format='fasta')
    msa = AlignIO.read(leaves_msa,'fasta')

    # Add initial MSA to the reconstructed ancestral sequences
    with open(os.path.join(path_to_ouput_internal_leaves_msa),'w') as f:
        for record in a_seq:
            print('>'+record.id, file=f)
            print(record.seq,file=f)
        for record in msa:
            print('>'+record.id, file=f)
            print(record.seq,file=f)

    msa_all = SeqIO.to_dict(SeqIO.parse(os.path.join(path_to_ouput_internal_leaves_msa),'fasta'))

    deletion = 0
    insertion = 0
    mutation = 0
    reinsertion = 0

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
    with open(os.path.join(path_to_output_internal_leaves_events),'w') as f1:
        for item in msa_all.keys():
            print('>'+item,file=f1)
            print(msa_all[item],file=f1)

    return(gene, tool, 'grasp', mutation, deletion, insertion, reinsertion)


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['IndelMAP', 'mafft', 'prank', 'ProPIP'] 

Indel_count_df = pd.DataFrame(columns=['gene', 'MSA', 'ASR', 'substitution', 'deletion', 'insertion', 'reinsertion'])

for gene in genes:
    for MSA_tool in MSA_tools:

        if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):
            indel_list = [gene, MSA_tool, 'grasp', 0,0,0,0]
            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_*.fasta'):
                subset_num = subset.split('_')[7].split('.')[0]
                if subset_num not in ['MSA','internal','all']:

                    ASR = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_'+subset_num+'.fasta'
                    leaves = '/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta'
                    tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_'+subset_num+'.nwk'

                    output_msa = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_MSA_all_'+subset_num+'.fasta'
                    output_events = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_all_events_'+subset_num+'.fasta'


                    indel_count = grasp_to_indelmap_format(ASR,leaves,tree,output_msa,output_events,gene,MSA_tool)
                    for i in [3,4,5,6]:
                        indel_list[i] = indel_list[i] + indel_count[i]

            Indel_count_df.loc[len(Indel_count_df)] = indel_list
            
        else:

            ASR = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp.fasta'
            leaves = '/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta'
            tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp.nwk'

            output_msa = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_MSA_all.fas'
            output_events = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_all_events.fas'


            Indel_count_df.loc[len(Indel_count_df)] = grasp_to_indelmap_format(ASR,leaves,tree,output_msa,output_events,gene,MSA_tool)

output_file = '/cfs/earth/scratch/seppemic/TM/results/others/grasp_indel_count.csv'
with open(output_file, 'w') as f:
    Indel_count_df.to_csv(f, sep='\t')

    