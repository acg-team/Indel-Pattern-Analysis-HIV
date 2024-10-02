import os
import pandas as pd
import numpy as np
from scipy.stats import wasserstein_distance
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import shapiro
from scipy.stats import ttest_ind, mannwhitneyu

def compare_msa_asr_differences(distance_matrix):
    # Extract MSA and ASR tools from the combination labels
    combinations = distance_matrix.index

    # msa_tools = ['Historian','IndelMAP','ProPIP','mafft','prank']
    # asr_tools = ['ArPIP', 'FastML', 'Historian', 'IndelMAP','grasp']

    # Initialize lists to store distances
    msa_diff_distances = []
    asr_diff_distances = []

    # Loop over all pairs of combinations
    for i, comb1 in enumerate(combinations):
        for j, comb2 in enumerate(combinations):

            if i < j:  # Avoid duplicates and diagonal (i==j)
                # Extract MSA and ASR tools for both combinations
                split1 = comb1.split('-')
                split2 = comb2.split('-')

                msa1, asr1 = split1[0], split1[1]
                msa2, asr2 = split2[0], split2[1]
                # Get the Wasserstein distance between the two combinations
                dist = distance_matrix.loc[comb1, comb2]

                # Check if MSA tool is different but ASR is the same
                if msa1 != msa2 and asr1 == asr2:
                    msa_diff_distances.append(dist)
                
                # Check if ASR tool is different but MSA is the same
                if asr1 != asr2 and msa1 == msa2:
                    asr_diff_distances.append(dist)

    # Convert to arrays for statistical testing
    msa_diff_distances = pd.Series(msa_diff_distances)
    asr_diff_distances = pd.Series(asr_diff_distances)

    # Compute mean distances
    mean_msa_diff = msa_diff_distances.mean()
    mean_asr_diff = asr_diff_distances.mean()

    data = msa_diff_distances
    stat, p_value = shapiro(data)
    alpha = 0.05  # significance level
    if p_value > alpha:
        print("The data is normally distributed, Perform a t-test")
        t_stat, p_value = ttest_ind(msa_diff_distances, asr_diff_distances, equal_var=False)

        return {
            'mean_msa_diff': mean_msa_diff,
            'mean_asr_diff': mean_asr_diff,
            't_stat': t_stat,
            'p_value': p_value
        }
    else:
        print("The data is not normally distributed, Perform a Mann-Whitney U test")
        u_stat, p_value = mannwhitneyu(msa_diff_distances, asr_diff_distances)

        return {
            'mean_msa_diff': mean_msa_diff,
            'mean_asr_diff': mean_asr_diff,
            'u_stat': u_stat,
            'p_value': p_value
        }

def plot_distance_matrix(distance_matrix, path_to_output_folder):
    # Set up the matplotlib figure
    plt.figure(figsize=(10, 8))
    
    # Create a heatmap with seaborn
    sns.heatmap(distance_matrix, annot=False, fmt=".3f", cmap="coolwarm", cbar=True, square=True, 
                linewidths=.5, xticklabels=True, yticklabels=True)
    
    # Add title and labels
    plt.title("Wasserstein Distance Matrix", fontsize=16)
    plt.xlabel("Combination", fontsize=12)
    plt.ylabel("Combination", fontsize=12)
    
    # Adjust layout to avoid cutting off labels
    plt.tight_layout()
    
    # Save the plot if save_path is provided
    if path_to_output_folder:
        plt.savefig(path_to_output_folder+'/test_wasserstein.png', format='png', dpi=300)  # Save at 300 dpi for high-quality output
        print(f"Plot saved at {path_to_output_folder}")

def compute_wasserstein_distance_matrix(df):
    # Group by 'combination' and get 'position in ref' for each group
    grouped = df.groupby('combination')['position in ref'].apply(list)
    
    # Get unique combinations
    combinations = grouped.index
    
    # Initialize the distance matrix
    distance_matrix = pd.DataFrame(np.zeros((len(combinations), len(combinations))), 
                                   index=combinations, columns=combinations)
    
    # Compute pairwise Wasserstein distances
    for i, comb1 in enumerate(combinations):
        for j, comb2 in enumerate(combinations):
            if i <= j:  # Only compute for upper triangle and diagonal (symmetric matrix)
                dist = wasserstein_distance(grouped[comb1], grouped[comb2])
                distance_matrix.loc[comb1, comb2] = dist
                distance_matrix.loc[comb2, comb1] = dist  # Symmetric
    
    return distance_matrix

def compute_wassertein(hotspots_df, path_to_output_folder, order, reg, env_specific):
    """Plot indel length for each protein or DNA region specified in order.

    Args:
        path_to_csv_length (string): path to csv file with indel lengths
        alphabet (string): 'protein' or 'DNA'
        order (list of strings): list with the proteins or DNA regions that should be plotted. 
                                 Example: ['signal peptide', 'gp120', 'gp41'] or 
                                 ['gp120', 'gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5']
        max_length (int): max indel length being plotted
    """


    if not env_specific :
        env_reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        hotspots_df['region'] = hotspots_df['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)
        only_hotspots = hotspots_df[['child','combination','region','event','position in ref']]

    else:
        only_hotspots = hotspots_df[['child','combination','region','event','position in ref']]
        regs_list_gp120 = ['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                                'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5']
        rpy_gp120 = only_hotspots.loc[only_hotspots['region'].isin(regs_list_gp120)]
        rpy_gp120.loc[only_hotspots['region'].isin(regs_list_gp120),'region'] = ['gp120']*len(rpy_gp120)
        only_hotspots = pd.concat([only_hotspots,rpy_gp120], axis=0, ignore_index=True)

        regs_list_env = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        rpy_env = only_hotspots.loc[only_hotspots['region'].isin(regs_list_env)]
        rpy_env.loc[only_hotspots['region'].isin(regs_list_env),'region'] = ['env']*len(rpy_env)
        only_hotspots = pd.concat([only_hotspots,rpy_env], axis=0, ignore_index=True)

    mask = only_hotspots['region'] == reg
    only_hotspots = only_hotspots[mask]


    only_position_indel = only_hotspots[only_hotspots['event'].isin(['insertion','deletion'])]

    only_position_indel.drop_duplicates(inplace=True)

    if not env_specific:
        start_region = 0
        if reg == 'env':
            end_region = 856
        elif reg == 'gag':
            end_region = 500
        elif reg == 'nef':
            end_region = 205
        elif reg == 'pol':
            end_region = 1003
        elif reg == 'rev':
            end_region = 116
        elif reg == 'tat':
            end_region = 100
        elif reg == 'vif':
            end_region = 192
        elif reg == 'vpr':
            end_region = 96
        elif reg == 'vpu':
            end_region = 82
    else :
        if reg == 'gp120':
            start_region = 30
            end_region = 510
        elif reg == 'gp41':
            start_region = 511
            end_region = 856
        elif reg == 'env signal peptide':
            start_region = 0
            end_region = 29
        elif reg == 'gp120 - C1':
            start_region = 30
            end_region = 128
        elif reg == 'gp120 - C2':
            start_region = 197
            end_region = 293
        elif reg == 'gp120 - C3':
            start_region = 332
            end_region = 382
        elif reg == 'gp120 - C4':
            start_region = 419
            end_region = 458
        elif reg == 'gp120 - C5':
            start_region = 470
            end_region = 510
        elif reg == 'gp120 - V1':
            start_region = 129
            end_region = 156
        elif reg == 'gp120 - V2':
            start_region = 157
            end_region = 196
        elif reg == 'gp120 - V3':
            start_region = 294
            end_region = 331
        elif reg == 'gp120 - V4':
            start_region = 383
            end_region = 418
        elif reg == 'gp120 - V5':
            start_region = 459
            end_region = 469

    print(only_position_indel)

    distance_matrix = compute_wasserstein_distance_matrix(only_position_indel)
    print(distance_matrix)

    plot_distance_matrix(distance_matrix, path_to_output_folder)

    results = compare_msa_asr_differences(distance_matrix)
    print(results)


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['Historian', 'IndelMAP','mafft','prank','ProPIP']
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP']


order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian',
            'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian',
            'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian',
            'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian',
            'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']



# Indel lengths vs region (gene), for each tool combination
hotspots_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')

        new_df['combination'] = MSA_tool+'-'+ASR_tool

        hotspots_df = pd.concat([hotspots_df, new_df], ignore_index =True)




path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/wasserstein'
alphabet = 'protein'
order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian',
            'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian',
            'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian',
            'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian',
            'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']


reg='gp120 - V1'
compute_wassertein(hotspots_df, path_to_output_folder, order, reg, env_specific=True)
