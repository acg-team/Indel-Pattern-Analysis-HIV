import subprocess as sp
import glob
import os
import random

# Function to replace ambiguous characters with a random amino acid
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

def process_fasta_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(file_path + '.NoAmb', 'w') as f:
        for line in lines:
            if line.startswith('>'):
                f.write(line)
            else:
                f.write(replace_ambiguous(line))

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

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['IndelMAP','mafft','prank','ProPIP']

for gene in genes:
    for MSA_tool in MSA_tools:
        if MSA_tool == 'ProPIP':
            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_*.fasta'):
                subset_num = subset.split('_')[4].split('.')[0]

                CMD = ['cp /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta']
                sp.run(CMD, shell=True)
                CMD = ['cp /cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk']
                sp.run(CMD, shell=True)

                CMD = ['export PATH=/cfs/earth/scratch/seppemic/TM/tools/GRASP/jdk-17.0.8/bin:$PATH ; ' + \
                        'java -jar /cfs/earth/scratch/seppemic/TM/tools/GRASP/bnkit.jar ' + \
                        '-a HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta ' + \
                        '-n HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk'
                        ]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp')

                # remove MSA and tree from output directory
                CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta']
                sp.run(CMD, shell=True)
                CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk']
                sp.run(CMD, shell=True)

                #rename output files
                CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'_ancestors.fa ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_'+subset_num+'.fasta']
                sp.run(CMD, shell=True)
                CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'_ancestors.nwk ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_'+subset_num+'.nwk']
                sp.run(CMD, shell=True)


        elif MSA_tool == 'prank' and gene == 'env':
            for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_*.fasta'):
                subset_num = subset.split('_')[4].split('.')[0]

                process_fasta_file('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta')

                CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta.NoAmb ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta.NoAmb']
                sp.run(CMD, shell=True)
                CMD = ['cp /cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk']
                sp.run(CMD, shell=True)

                CMD = ['export PATH=/cfs/earth/scratch/seppemic/TM/tools/GRASP/jdk-17.0.8/bin:$PATH ; ' + \
                        'java -jar /cfs/earth/scratch/seppemic/TM/tools/GRASP/bnkit.jar ' + \
                        '-a HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta.NoAmb ' + \
                        '-n HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk'
                        ]
                subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp')

                # remove MSA and tree from output directory
                CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'.fasta.NoAmb']
                sp.run(CMD, shell=True)
                CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_timetree_subset_'+subset_num+'.nwk']
                sp.run(CMD, shell=True)

                #rename output files
                CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'_ancestors.fa ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_'+subset_num+'.fasta']
                sp.run(CMD, shell=True)
                CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_'+subset_num+'_ancestors.nwk ' + \
                        '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp_'+subset_num+'.nwk']
                sp.run(CMD, shell=True)


        else:
            process_fasta_file('/cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta')

            CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/MSA/'+gene+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta.NoAmb ' + \
                    '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta.NoAmb']
            sp.run(CMD, shell=True)

            CMD = ['cp /cfs/earth/scratch/seppemic/TM/data/'+gene+'/HIV-1_'+gene+'_DNA_mafft_alignment.fasta.timetree.nwk ' + \
                    '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_DNA_mafft_alignment.fasta.timetree.nwk']
            sp.run(CMD, shell=True)

            CMD = ['export PATH=/cfs/earth/scratch/seppemic/TM/tools/GRASP/jdk-17.0.8/bin:$PATH ; ' + \
                    'java -jar /cfs/earth/scratch/seppemic/TM/tools/GRASP/bnkit.jar ' + \
                    '-a HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta.NoAmb ' + \
                    '-n HIV-1_'+gene+'_DNA_mafft_alignment.fasta.timetree.nwk'
                    ]
            subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp')

            # remove MSA and tree from output directory
            CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'.fasta.NoAmb']
            sp.run(CMD, shell=True)
            CMD = ['rm /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_DNA_mafft_alignment.fasta.timetree.nwk']
            sp.run(CMD, shell=True)

            #rename output files
            CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_ancestors.fa ' + \
                    '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp.fasta']
            sp.run(CMD, shell=True)
            CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_ancestors.nwk ' + \
                    '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_grasp/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_grasp.nwk']
            sp.run(CMD, shell=True)

