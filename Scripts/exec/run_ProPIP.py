import subprocess as sp
import glob

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

genes = ['rev'] # 'env','gag','nef','pol','rev','tat', 'vif','vpr','vpu'

for gene in genes:

    for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/'+gene+'_sub_1:244.fasta'):

        subset_num = subset.split('_')[3].split('.')[0]

        #Modify path to the input and output files in the ProPIP parameter text file
        with open('/cfs/earth/scratch/seppemic/TM/scripts/exec/ProPIP_param.txt', 'r') as f:
            lines = f.readlines()
        with open('/cfs/earth/scratch/seppemic/TM/scripts/exec/ProPIP_param.txt', 'w') as f:
            for line in lines:
                if line.startswith('input.sequence.file='):
                    f.write('input.sequence.file=' + subset + '\n')
                elif line.startswith('input.tree.file='):
                    f.write('input.tree.file=/cfs/earth/scratch/seppemic/TM/data/'+gene+'_subset/HIV-1_' + gene + '_timetree_subset_' + subset_num + '.nwk\n')
                elif line.startswith('output.msa.file='):
                    f.write('output.msa.file=/cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_ProPIP_' + subset_num + '.fasta\n')
                else:
                    f.write(line)

        CMD = ["/cfs/earth/scratch/seppemic/TM/tools/ProPIP/ProPIP params=/cfs/earth/scratch/seppemic/TM/scripts/exec/ProPIP_param.txt"]
        subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

        CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_ProPIP_' + subset_num + '.initial.fasta ' + \
                '/cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_ProPIP_' + subset_num + '.fasta']
        subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')
