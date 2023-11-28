from helper_functions import *
from write_evolutionary_events import *
from calculate_indel_rates import *
from calculate_indel_lengths import *
import subprocess as sp
import glob

def indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate):

    # compute events dataframe
    indel_df = write_evolutionary_events(path_to_events_all, path_to_nwk, sub_rate, gene)
    output_file = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv'
    with open(output_file, 'a') as f:
        indel_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

    # compute lengths dataframe
    lengths_df = calculate_indel_lengths(indel_df)
    output_length_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv'
    with open(output_length_df, 'a') as f:
        lengths_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

    # compute rates dataframe
    rates_df = calculate_indel_rates(indel_df, path_to_events_all, path_to_nwk, gene)
    output_rates_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_indel_df_rates.csv'
    with open(output_rates_df, 'a') as f:
        rates_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

def root_tree(unrooted_tree, rooted_tree):
    tree = PhyloNode(newick = unrooted_tree, format=1)
    with open(rooted_tree, 'w') as f:
            print(tree.write(format=1).split(';')[0]+"ROOT;", file=f)


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu'] # 'env','gag','nef','pol','rev','tat','vif','vpr','vpu'
MSA_tools = ['IndelMAP','mafft','prank','ProPIP'] # 'IndelMAP','mafft','prank','ProPIP'
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'IndelMAP'] # 'ArPIP', 'FastML', 'grasp', 'IndelMAP'

sub_rates_dic = {'env': 0.0057326, 
                'gag': 0.00226134,
                'nef': 0.00256465,
                'pol': 0.000902068,
                'rev': 0.00546282,
                'tat': 0.00557396,
                'vif': 0.00125458,
                'vpr': 0.00122074,
                'vpu': 0.00208927}

for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        CMD = ['rm -r /cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool]
        sp.run(CMD, shell=True)
        CMD = ['mkdir /cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool]
        sp.run(CMD, shell=True)


# Iterate through all genes, MSA and ASR and append indel informations to output file
for gene in genes:
    sub_rate = sub_rates_dic[gene]
    for MSA_tool in MSA_tools:
        for ASR_tool in ['FastML','grasp']: # 'FastML','grasp'

            if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):

                if ASR_tool == 'FastML':

                    for split in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/subset_*'):

                        path_to_events_all = split + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fasta'
                        path_to_nwk = split + '/tree.rewritten.newick.txt'

                        indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)  

                else :

                    for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events_*.fasta'):
                        subset_num = subset.split('_')[9].split('.')[0]

                        path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events_'+subset_num+'.fasta'
                        path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'.nwk'

                        indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

            elif ASR_tool == 'FastML':

                for split in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/split_*'):

                    path_to_events_all = split + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fasta'
                    path_to_nwk = split + '/tree.rewritten.newick.txt'

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)  
                
            else:

                path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fas'
                path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'.nwk'

                indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

            print( 'done for : ' + gene + ' - ' + MSA_tool + ' - ' + ASR_tool)

## specific for ArPIP
for gene in genes:
    sub_rate = sub_rates_dic[gene]
    for MSA_tool in MSA_tools:
        for ASR_tool in ['ArPIP']:

            if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):

                for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_*_internal_evolutionary_events.*.fas'):
                    subset_num = subset.split('_')[7]
                    split_num = subset.split('.')[1]
                    path_to_events_all = subset
                    path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'.'+split_num+'.nwk'

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)
  
            else:
                for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_internal_evolutionary_events.*.fas'):
                    subset_num = subset.split('.')[1]
                    path_to_events_all = subset
                    path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'.'+subset_num+'.nwk'

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

            print( 'done for : ' + gene + ' - ' + MSA_tool + ' - ' + ASR_tool)


## specific for IndelMAP
for gene in genes:
    sub_rate = sub_rates_dic[gene]
    for MSA_tool in MSA_tools:
        for ASR_tool in ['IndelMAP']:

            if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):

                for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_*_all_events.fas'):
                    subset_num = subset.split('_')[7]

                    unrooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'_tree.nwk'
                    rooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'_tree_rooted.nwk'
                    root_tree(unrooted_tree, rooted_tree)

                    path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'_all_events.fas'
                    path_to_nwk = rooted_tree

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

            else:

                unrooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_tree.nwk'
                rooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_tree_rooted.nwk'
                root_tree(unrooted_tree, rooted_tree)

                path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fas'
                path_to_nwk = rooted_tree

                indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

            print( 'done for : ' + gene + ' - ' + MSA_tool + ' - ' + ASR_tool)

