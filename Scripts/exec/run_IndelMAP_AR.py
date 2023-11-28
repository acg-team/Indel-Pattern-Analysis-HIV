import subprocess as sp
import glob
from Bio import AlignIO
from ete3 import Tree
import os

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

    t.write(format=5, outfile = tree_name + '_processed.nwk')


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['IndelMAP', 'mafft', 'prank', 'ProPIP'] # 'IndelMAP', 'mafft', 'prank', 'ProPIP'

for gene in genes:
    for MSA_tool in MSA_tools:

        if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):
            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_*.fasta'):
                subset_num = subset.split('_')[4].split('.')[0]

                process_tree('/cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk')

                CMD = ['python /cfs/earth/scratch/seppemic/TM/tools/indelMaP/indelMaP_ASR.py ' + \
                                        '--msa_file /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta ' + \
                                        '--tree_file /cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'_processed.nwk ' + \
                                        '--alphabet Protein ' + \
                                        '--RateMatrix HIVb ' + \
                                        '--output_file /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP_'+subset_num
                                        ]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
 
                # append internal and external nodes to the same MSA
                ancestral_sequences = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP_'+subset_num+'_internal_evolutionary_events.fas'
                leaves_msa = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP_'+subset_num+'_leaves_evolutionary_events.fas'
                ouput_internal_leaves_msa = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP_'+subset_num+'_all_events.fas'
                a_seq = AlignIO.read(ancestral_sequences,format='fasta')
                msa = AlignIO.read(leaves_msa,'fasta')
                with open(os.path.join(ouput_internal_leaves_msa),'w') as f:
                    for record in a_seq:
                        print('>'+record.id, file=f)
                        print(record.seq,file=f)
                    for record in msa:
                        print('>'+record.id, file=f)
                        print(record.seq,file=f)

        else:
            CMD = ['python /cfs/earth/scratch/seppemic/TM/tools/indelMaP/indelMaP_ASR.py ' + \
                                    '--msa_file /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta ' + \
                                    '--tree_file /cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment.fasta.timetree.nwk ' + \
                                    '--alphabet Protein ' + \
                                    '--RateMatrix HIVb ' + \
                                    '--output_file /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP'
                                    ]
            subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

            # append internal and external nodes to the same MSA
            ancestral_sequences = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP_internal_evolutionary_events.fas'
            leaves_msa = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP_leaves_evolutionary_events.fas'
            ouput_internal_leaves_msa = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_IndelMAP/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_IndelMAP_all_events.fas'
            a_seq = AlignIO.read(ancestral_sequences,format='fasta')
            msa = AlignIO.read(leaves_msa,'fasta')
            with open(os.path.join(ouput_internal_leaves_msa),'w') as f:
                for record in a_seq:
                    print('>'+record.id, file=f)
                    print(record.seq,file=f)
                for record in msa:
                    print('>'+record.id, file=f)
                    print(record.seq,file=f)

