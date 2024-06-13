import subprocess as sp
import glob

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']

for gene in genes:

    if gene == 'env':
        
        for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/data/' + gene + '_subset/' + gene + '_sub_*.fasta'):
            subset_num = subset.split('_')[3].split('.')[0]
            CMD = ['prank -d=/cfs/earth/scratch/seppemic/TM/data/' + gene + '_subset/' + gene + '_sub_' + subset_num + '.fasta ' + \
                    '-t=/cfs/earth/scratch/seppemic/TM/data/' + gene + '_subset/HIV-1_' + gene + '_timetree_subset_' + subset_num + '.nwk ' + \
                    '-o=/cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_prank_' + subset_num + ' ' + \
                    '-protein ' + \
                    '-once ' + \
                    '+F'
                    ]
            sp.run(CMD, shell=True)

            CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_prank_' + subset_num + '.best.fas ' + \
                '/cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_prank_' + subset_num + '.fasta']
            sp.run(CMD, '/cfs/earth/scratch/seppemic/TM')

    else:

        CMD = ['prank -d=/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_protein_dated.fasta ' + \
                '-t=/cfs/earth/scratch/seppemic/TM/data/' + gene + '/HIV-1_' + gene + '_DNA_mafft_alignment.fasta.timetree.nwk ' + \
                '-o=/cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_prank ' + \
                '-protein ' + \
                '-once ' + \
                '+F'
                ]
        sp.run(CMD, shell=True)

        CMD = ['mv /cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_prank.best.fas ' + \
                '/cfs/earth/scratch/seppemic/TM/results/MSA/' + gene + '/HIV-1_' + gene + '_MSA_prank.fasta']
        sp.run(CMD, '/cfs/earth/scratch/seppemic/TM')


