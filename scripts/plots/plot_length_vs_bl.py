import pandas as pd
import seaborn as sns
import seaborn_image as isns
import matplotlib.pyplot as plt

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu'] 
MSA_tools = ['Historian', 'IndelMAP','mafft','prank','ProPIP']
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP']


def plot_indel_length_vs_branch_length_scatter(lengths_df, path_to_output_folder):

    plt.figure(figsize=(12, 8))

    g = sns.FacetGrid(lengths_df, col='MSA', hue='ASR', palette='terrain', height=4, aspect=1.2, col_wrap=2)

    # Scatter plot in each subplot
    g.map(sns.scatterplot, 'tMRCA', 'length', s=50, alpha=0.8)

    # Regression line in each subplot
    # g.map(sns.regplot, 'tMRCA', 'indel_count', scatter=False, ci=None, line_kws={'lw': 1}, lowess=True)

    # Set labels and legend
    g.set_axis_labels('Branch Length (years)', 'Indel length')

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(path_to_output_folder + '/indel_vs_branch_length_scatter.png',transparent=True)

def plot_indel_length_vs_branch_length_reg(lengths_df, path_to_output_folder):

    plt.figure(figsize=(12, 8))

    g = sns.FacetGrid(lengths_df, col='MSA', hue='ASR', palette='terrain', height=4, aspect=1.2, col_wrap=2)

    g.map(sns.regplot, 'tMRCA', 'length', scatter=False, line_kws={'lw': 1.5})

    # Set labels and legend
    g.set_axis_labels('Branch Length (years)', 'Indel length')

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(path_to_output_folder + '/indel_vs_branch_length_reg.png',transparent=True)


length_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv',sep='\t')
        new_df['MSA'] = MSA_tool
        new_df['ASR'] = ASR_tool
        length_df = pd.concat([length_df, new_df], ignore_index =True)

path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/length_vs_bl'
plot_indel_length_vs_branch_length_scatter(length_df, path_to_output_folder)
plot_indel_length_vs_branch_length_reg(length_df, path_to_output_folder)