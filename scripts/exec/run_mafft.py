import subprocess as sp

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

def rename_seq(input_fasta, input_tree, output_tree):
    '''
    Rewrites the input tree file (fasta) with integers instead of seq ID (order is conserved).
    
    Parameters
    -------
    input_fasta :   string
            Full path to the input file with simulated sequences, corresponding to the tree of interest (fasta).
    input_tree :    string
            Full path to the input tree file (newick).
    output_tree :   string
            Full path to the output tree file (newick).
    
    Returns
    -------
    None
    
    ''' 
    dic = {}
    entry = 1
    # fill dict with integer and corresponding sequence ID
    with open(input_fasta, 'r') as file:
        for line in file:
            line = line.rstrip('\n')
            if line.startswith('>'):
                dic[entry] = line.strip('>')
                entry += 1


    tree_file = open(input_tree)

    # rewrite each sequence ID by its corresponding integer in the input tree newick format
    for tree in tree_file:
    
        f = open(output_tree, "w+")
        if not dic[1].isnumeric() :
            for entry in dic:
                tree = tree.replace(dic[entry], str(entry))
    
        f.write(tree)
        f.close()

    tree_file.close()

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']

for gene in genes:

    # rewrite tree file with integers instead of sequence names
    rename_seq('/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_protein_dated.fasta', \
                '/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.nwk', \
                '/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.num.nwk')

    # convert the newick tree file to a mafft compatible format
    CMD = ['ruby /cfs/earth/scratch/seppemic/TM/scripts/exec/newick2mafft.rb ' + \
            '/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.num.nwk > ' + \
            '/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.mafft'
            ]
    
    subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')


    CMD = ['mafft ' + \
            '--auto ' + \
            '--treein /cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.mafft ' + \
            '--aamatrix /cfs/earth/scratch/seppemic/TM/data/Between10 ' + \
            '/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_protein_dated.fasta > ' + \
            '/cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_mafft.fasta'
            ]

    subprocess_cmd(CMD, '/cfs/earth/scratch/seppemic/TM')

    print(gene + ' done')
