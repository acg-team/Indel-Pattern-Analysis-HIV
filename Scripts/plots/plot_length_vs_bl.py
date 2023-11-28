import pandas as pd
import seaborn as sns
import seaborn_image as isns
import matplotlib.pyplot as plt

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


def plot_indel_length_vs_branch_length_scatter(lengths_df, path_to_output_folder):
    plt.figure(figsize=(10, 6))

    sns.scatterplot(x='tMRCA', y='length', data=lengths_df, hue='combination', palette='viridis', style='combination', s=50, alpha=0.8)
    plt.xlabel('Branch Length (years)')
    plt.ylabel('Indel Length')
    # plt.yscale('log')

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(path_to_output_folder + '/indel_vs_branch_length_scatter_full.png')

def plot_indel_length_vs_branch_length_reg(lengths_df, path_to_output_folder):
    plt.figure(figsize=(10, 6))

    scatter_plot = sns.scatterplot(x='tMRCA', y='length', data=lengths_df, hue='combination', palette='viridis', style='combination', s=50, alpha=0.8)
    sns.lmplot(x='tMRCA', y='length', data=lengths_df, hue='combination', scatter=False, line_kws={'lw': 1.5}, legend=False)
    plt.xlabel('Branch Length (years)')
    plt.ylabel('Indel Length')

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(path_to_output_folder + '/indel_vs_branch_length_reg_full.png')



length_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        new_df = pd.read_csv('/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv',sep='\t')
        new_df['combination'] = MSA_tool+'-'+ASR_tool
        length_df = pd.concat([length_df, new_df], ignore_index =True)


path_to_output_folder = '/cfs/earth/scratch/seppemic/TM/results/plots/length_vs_bl'
plot_indel_length_vs_branch_length_scatter(length_df, path_to_output_folder)
plot_indel_length_vs_branch_length_reg(length_df, path_to_output_folder)



length_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        new_df = pd.read_csv('/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv',sep='\t')
        new_df['combination'] = ASR_tool
        length_df = pd.concat([length_df, new_df], ignore_index =True)

path_to_output_folder = '/cfs/earth/scratch/seppemic/TM/results/plots/length_vs_bl_vs_ASR'
plot_indel_length_vs_branch_length_scatter(length_df, path_to_output_folder)
plot_indel_length_vs_branch_length_reg(length_df, path_to_output_folder)

