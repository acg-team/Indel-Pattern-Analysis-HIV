from helper_functions import *
from write_evolutionary_events import *
from calculate_indel_rates import *
from calculate_indel_lengths import *
import subprocess as sp
import glob
from Bio import Phylo, SeqIO
from io import StringIO
import matplotlib.pyplot as plt

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


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['IndelMAP','mafft','prank','ProPIP']
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'IndelMAP']
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

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

            file_paths = []

            if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):

                if ASR_tool == 'FastML':

                    all_sequences = []

                    # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                    reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                    count_df = pd.DataFrame(index=amino_acids, columns=reg)
                    count_df[reg] = 0

                    for split in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/subset_*'):

                        path_to_events_all = split + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fasta'
                        path_to_nwk = split + '/tree.rewritten.newick.txt'

                        indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)  

                        # append subset MSA to general -> For background frequencies only
                        # Open MSA file
                        with open(path_to_events_all, "r") as msa_file:
                            msa_records = list(SeqIO.parse(msa_file, "fasta"))

                        # Open tree file
                        with open(path_to_nwk, "r") as tree_file:
                            tree_content = tree_file.read()
                            tree = Phylo.read(StringIO(tree_content), "newick")

                        # Extract external node names from the tree
                        external_nodes = [terminal.name for terminal in tree.get_terminals()]

                        # Extract sequences corresponding to external nodes from the MSA
                        selected_sequences = [record for record in msa_records if record.id in external_nodes]

                        # Append selected sequences to the list
                        all_sequences.extend(selected_sequences)

                        if gene == 'env':

                            # get specific regions of env and compute background frequencies for those regions
                            start_reg, end_reg = extract_regions(path_to_events_all)

                            for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                                extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                                new_count_df = append_count_aa(extracted_region, count_df, reg[idx])



                    freq_df = compute_background_frequencies(all_sequences)
                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                    if gene == 'env':
                        for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                            total_count = new_count_df[reg[idx]].sum()

                            if total_count !=0 :
                                new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                            else:
                                print('div by 0 !!!')
                                print(new_count_df)

                        output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                        with open(output_freq_df, 'a') as f:
                            new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')


                else :

                    all_sequences = []
                    # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                    reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                    count_df = pd.DataFrame(index=amino_acids, columns=reg)
                    count_df[reg] = 0

                    for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events_*.fasta'):
                        subset_num = subset.split('_')[9].split('.')[0]

                        path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events_'+subset_num+'.fasta'
                        path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'.nwk'

                        indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

                        # append subset MSA to general -> For background frequencies only
                        # Open MSA file
                        with open(path_to_events_all, "r") as msa_file:
                            msa_records = list(SeqIO.parse(msa_file, "fasta"))

                        # Open tree file
                        with open(path_to_nwk, "r") as tree_file:
                            tree_content = tree_file.read()
                            tree = Phylo.read(StringIO(tree_content), "newick")

                        # Extract external node names from the tree
                        external_nodes = [terminal.name for terminal in tree.get_terminals()]

                        # Extract sequences corresponding to external nodes from the MSA
                        selected_sequences = [record for record in msa_records if record.id in external_nodes]

                        # Append selected sequences to the list
                        all_sequences.extend(selected_sequences)

                        if gene == 'env':

                            # get specific regions of env and compute background frequencies for those regions
                            start_reg, end_reg = extract_regions(path_to_events_all)

                            for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                                extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                                new_count_df = append_count_aa(extracted_region, count_df, reg[idx])



                    freq_df = compute_background_frequencies(all_sequences)
                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                    if gene == 'env':
                        for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                            total_count = new_count_df[reg[idx]].sum()

                            if total_count !=0 :
                                new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                            else:
                                print('div by 0 !!!')
                                print(new_count_df)

                        output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                        with open(output_freq_df, 'a') as f:
                            new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')



            elif ASR_tool == 'FastML':

                all_sequences = []
                # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                count_df = pd.DataFrame(index=amino_acids, columns=reg)
                count_df[reg] = 0
                print(str(gene) + ' ' + str(MSA_tool) + ' ' + str(ASR_tool))
                print(count_df)

                for split in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/split_*'):

                    path_to_events_all = split + '/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fasta'
                    path_to_nwk = split + '/tree.rewritten.newick.txt'

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)  

                    # append subset MSA to general -> For background frequencies only
                    # Open MSA file
                    with open(path_to_events_all, "r") as msa_file:
                        msa_records = list(SeqIO.parse(msa_file, "fasta"))

                    # Open tree file
                    with open(path_to_nwk, "r") as tree_file:
                        tree_content = tree_file.read()
                        tree = Phylo.read(StringIO(tree_content), "newick")

                    # Extract external node names from the tree
                    external_nodes = [terminal.name for terminal in tree.get_terminals()]

                    # Extract sequences corresponding to external nodes from the MSA
                    selected_sequences = [record for record in msa_records if record.id in external_nodes]

                    # Append selected sequences to the list
                    all_sequences.extend(selected_sequences)

                    if gene == 'env':

                        # get specific regions of env and compute background frequencies for those regions
                        start_reg, end_reg = extract_regions(path_to_events_all)
                        print(start_reg)
                        print(end_reg)
                        print(external_nodes)

                        for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                            extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                            print(extracted_region)
                            new_count_df = append_count_aa(extracted_region, count_df, reg[idx])


                freq_df = compute_background_frequencies(all_sequences)
                output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                with open(output_freq_df, 'a') as f:
                    freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                if gene == 'env':
                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                        print('pas ok ??')
                        print(reg[idx])
                        print(new_count_df)
                        total_count = new_count_df.loc[:,reg[idx]].sum()
                        print('ok... Total Count : ' + str(total_count))

                        if total_count !=0 :
                            new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                        else:
                            print('div by 0 !!!')
                            print(new_count_df)

                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                
            else:

                # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                count_df = pd.DataFrame(index=amino_acids, columns=reg)
                count_df[reg] = 0

                path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fas'
                path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'.nwk'

                indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

                # append subset MSA to general -> For background frequencies only
                # Open MSA file
                with open(path_to_events_all, "r") as msa_file:
                    msa_records = list(SeqIO.parse(msa_file, "fasta"))

                # Open tree file
                with open(path_to_nwk, "r") as tree_file:
                    tree_content = tree_file.read()
                    tree = Phylo.read(StringIO(tree_content), "newick")

                # Extract external node names from the tree
                external_nodes = [terminal.name for terminal in tree.get_terminals()]

                # Extract sequences corresponding to external nodes from the MSA
                all_sequences = [record for record in msa_records if record.id in external_nodes]

                if gene == 'env':

                    # get specific regions of env and compute background frequencies for those regions
                    start_reg, end_reg = extract_regions(path_to_events_all)

                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                        extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                        new_count_df = append_count_aa(extracted_region, count_df, reg[idx])


                freq_df = compute_background_frequencies(all_sequences)
                output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                with open(output_freq_df, 'a') as f:
                    freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                if gene == 'env':
                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                        total_count = new_count_df[reg[idx]].sum()

                        if total_count !=0 :
                            new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                        else:
                            print('div by 0 !!!')
                            print(new_count_df)

                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')


            print( 'done for : ' + gene + ' - ' + MSA_tool + ' - ' + ASR_tool)

## specific for ArPIP
for gene in genes:
    sub_rate = sub_rates_dic[gene]
    for MSA_tool in MSA_tools:
        for ASR_tool in ['ArPIP']:

            file_paths = []

            if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):

                all_sequences = []
                # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                count_df = pd.DataFrame(index=amino_acids, columns=reg)
                count_df[reg] = 0

                for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_*_internal_evolutionary_events.*.fas'):
                    subset_num = subset.split('_')[7]
                    split_num = subset.split('.')[1]
                    path_to_events_all = subset
                    path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'.'+split_num+'.nwk'

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

                    # Open MSA file
                    with open(path_to_events_all, "r") as msa_file:
                        msa_records = list(SeqIO.parse(msa_file, "fasta"))

                    # Open tree file
                    with open(path_to_nwk, "r") as tree_file:
                        tree_content = tree_file.read()
                        tree = Phylo.read(StringIO(tree_content), "newick")

                    # Extract external node names from the tree
                    external_nodes = [terminal.name for terminal in tree.get_terminals()]

                    # Extract sequences corresponding to external nodes from the MSA
                    selected_sequences = [record for record in msa_records if record.id in external_nodes]

                    # Append selected sequences to the list
                    all_sequences.extend(selected_sequences)

                    if gene == 'env':

                        # get specific regions of env and compute background frequencies for those regions
                        start_reg, end_reg = extract_regions(path_to_events_all)

                        for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                            extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                            new_count_df = append_count_aa(extracted_region, count_df, reg[idx])



                freq_df = compute_background_frequencies(all_sequences)
                output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                with open(output_freq_df, 'a') as f:
                    freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                if gene == 'env':
                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                        total_count = new_count_df[reg[idx]].sum()

                        if total_count !=0 :
                            new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                        else:
                            print('div by 0 !!!')
                            print(new_count_df)

                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

  
            else:

                all_sequences = []
                # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                count_df = pd.DataFrame(index=amino_acids, columns=reg)
                count_df[reg] = 0

                for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_internal_evolutionary_events.*.fas'):
                    subset_num = subset.split('.')[1]
                    path_to_events_all = subset
                    path_to_nwk = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'.'+subset_num+'.nwk'

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

                    # append subset MSA to general -> For background frequencies only
                    # Open MSA file
                    with open(path_to_events_all, "r") as msa_file:
                        msa_records = list(SeqIO.parse(msa_file, "fasta"))

                    # Open tree file
                    with open(path_to_nwk, "r") as tree_file:
                        tree_content = tree_file.read()
                        tree = Phylo.read(StringIO(tree_content), "newick")

                    # Extract external node names from the tree
                    external_nodes = [terminal.name for terminal in tree.get_terminals()]

                    # Extract sequences corresponding to external nodes from the MSA
                    selected_sequences = [record for record in msa_records if record.id in external_nodes]

                    # Append selected sequences to the list
                    all_sequences.extend(selected_sequences)

                    if gene == 'env':

                        # get specific regions of env and compute background frequencies for those regions
                        start_reg, end_reg = extract_regions(path_to_events_all)

                        for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                            extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                            new_count_df = append_count_aa(extracted_region, count_df, reg[idx])



                freq_df = compute_background_frequencies(all_sequences)
                output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                with open(output_freq_df, 'a') as f:
                    freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                if gene == 'env':
                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                        total_count = new_count_df[reg[idx]].sum()

                        if total_count !=0 :
                            new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                        else:
                            print('div by 0 !!!')
                            print(new_count_df)

                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')


            print( 'done for : ' + gene + ' - ' + MSA_tool + ' - ' + ASR_tool)


## specific for IndelMAP
for gene in genes:
    sub_rate = sub_rates_dic[gene]
    for MSA_tool in MSA_tools:
        for ASR_tool in ['IndelMAP']:

            file_paths = []

            if MSA_tool == 'ProPIP' or (MSA_tool == 'prank' and gene == 'env'):

                all_sequences = []
                # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                count_df = pd.DataFrame(index=amino_acids, columns=reg)
                count_df[reg] = 0

                for subset in glob.glob('/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_*_all_events.fas'):
                    subset_num = subset.split('_')[7]

                    unrooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'_tree.nwk'
                    rooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'_tree_rooted.nwk'
                    root_tree(unrooted_tree, rooted_tree)

                    path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_'+subset_num+'_all_events.fas'
                    path_to_nwk = rooted_tree

                    indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

                    # append subset MSA to general -> For background frequencies only
                    # Open MSA file
                    with open(path_to_events_all, "r") as msa_file:
                        msa_records = list(SeqIO.parse(msa_file, "fasta"))

                    # Open tree file
                    with open(path_to_nwk, "r") as tree_file:
                        tree_content = tree_file.read()
                        tree = Phylo.read(StringIO(tree_content), "newick")

                    # Extract external node names from the tree
                    external_nodes = [terminal.name for terminal in tree.get_terminals()]

                    # Extract sequences corresponding to external nodes from the MSA
                    selected_sequences = [record for record in msa_records if record.id in external_nodes]

                    # Append selected sequences to the list
                    all_sequences.extend(selected_sequences)

                    if gene == 'env':

                        # get specific regions of env and compute background frequencies for those regions
                        start_reg, end_reg = extract_regions(path_to_events_all)

                        for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                            extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                            new_count_df = append_count_aa(extracted_region, count_df, reg[idx])



                freq_df = compute_background_frequencies(all_sequences)
                output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                with open(output_freq_df, 'a') as f:
                    freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                if gene == 'env':
                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                        total_count = new_count_df[reg[idx]].sum()

                        if total_count !=0 :
                            new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                        else:
                            print('div by 0 !!!')
                            print(new_count_df)

                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')


            else:

                # reg = ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']
                reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
                count_df = pd.DataFrame(index=amino_acids, columns=reg)
                count_df[reg] = 0

                unrooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_tree.nwk'
                rooted_tree = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_tree_rooted.nwk'
                root_tree(unrooted_tree, rooted_tree)

                path_to_events_all = '/cfs/earth/scratch/seppemic/TM/results/AR/'+gene+'/'+MSA_tool+'_'+ASR_tool+'/HIV-1_'+gene+'_MSA_'+MSA_tool+'_AR_'+ASR_tool+'_all_events.fas'
                path_to_nwk = rooted_tree

                indel_info_extraction(path_to_events_all,path_to_nwk,gene,MSA_tool,ASR_tool,sub_rate)

                # append subset MSA to general -> For background frequencies only
                # Open MSA file
                with open(path_to_events_all, "r") as msa_file:
                    msa_records = list(SeqIO.parse(msa_file, "fasta"))

                # Open tree file
                with open(path_to_nwk, "r") as tree_file:
                    tree_content = tree_file.read()
                    tree = Phylo.read(StringIO(tree_content), "newick")

                # Extract external node names from the tree
                external_nodes = [terminal.name for terminal in tree.get_terminals()]

                # Extract sequences corresponding to external nodes from the MSA
                all_sequences = [record for record in msa_records if record.id in external_nodes]

                if gene == 'env':

                    # get specific regions of env and compute background frequencies for those regions
                    start_reg, end_reg = extract_regions(path_to_events_all)

                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):
                        extracted_region = extract_region_from_msa(path_to_events_all, start_reg[env_reg], end_reg[env_reg], external_nodes)
                        new_count_df = append_count_aa(extracted_region, count_df, reg[idx])


                freq_df = compute_background_frequencies(all_sequences)
                output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv'
                with open(output_freq_df, 'a') as f:
                    freq_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')

                if gene == 'env':
                    for idx, env_reg in enumerate([0,2,4,6,8,10,12,14,16,18,20,22]):

                        total_count = new_count_df[reg[idx]].sum()

                        if total_count !=0 :
                            new_count_df[reg[idx] + ' freq'] = new_count_df[reg[idx]] / total_count
                        else:
                            print('div by 0 !!!')
                            print(new_count_df)

                    output_freq_df = '/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv'
                    with open(output_freq_df, 'a') as f:
                        new_count_df.to_csv(f, mode='a', header=f.tell()==0, sep='\t')


            print( 'done for : ' + gene + ' - ' + MSA_tool + ' - ' + ASR_tool)
