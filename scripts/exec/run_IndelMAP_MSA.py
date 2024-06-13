import subprocess as sp
from ete3 import Tree
import re
from decimal import Decimal

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

def convert_scientific_to_decimal(input_file, output_file):
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        for line in in_file:
            line = re.sub(r"(\d+\.\d+[eE][-+]?\d+(?=[,()])|\d+\.\d+(?=[,()]))", lambda match: format(Decimal(match.group(0)), 'f'), line)
            out_file.write(line)

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']

for gene in genes:

    convert_scientific_to_decimal('/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.nwk','/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment_processed.nwk')

    CMD = ['/cfs/earth/scratch/seppemic/TM/tools/rust-indelMaP/indelMaP/target/release/indelMaP ' + \
                            '--seq-file /cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_protein_dated.fasta ' + \
                            '--tree-file /cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment_processed.nwk ' + \
                            '--model HIVb ' + \
                            '--output-msa-file /cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_IndelMAP.fasta'
                            ]
    subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

    CMD = ['rm /cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment_processed.nwk']
    subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
