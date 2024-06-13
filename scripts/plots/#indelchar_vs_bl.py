import pandas as pd
import seaborn as sns
import seaborn_image as isns
import matplotlib.pyplot as plt

genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu'] # 'env','gag','nef','pol','rev','tat','vif','vpr','vpu'
MSA_tools = ['Historian', 'IndelMAP','mafft','prank','ProPIP'] # 'Historian', 'IndelMAP','mafft','prank','ProPIP'
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP'] # 'ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP'


def plot_number_of_indel_char_vs_bl(indel_df, path_to_output_folder):

    indel_df = indel_df[indel_df['event'].isin(['insertion','deletion'])]
    indel_df = indel_df[['child','tMRCA','ASR','MSA']]

    indel_df['indel_count'] = 1
    grouped_df = indel_df.groupby(['child','ASR','MSA']).agg({'tMRCA': 'mean', 'indel_count': 'sum'}).reset_index()


    plt.figure(figsize=(12, 8))

    g = sns.FacetGrid(grouped_df, col='MSA', hue='ASR', palette='terrain', height=4, aspect=1.2, col_wrap=2)

    # Scatter plot in each subplot
    g.map(sns.scatterplot, 'tMRCA', 'indel_count', s=50, alpha=0.8)

    # Regression line in each subplot
    # g.map(sns.regplot, 'tMRCA', 'indel_count', scatter=False, ci=None, line_kws={'lw': 1}, lowess=True)

    # Set labels and legend
    g.set_axis_labels('Branch Length (years)', 'Inserted/Deleted Residues Count')

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(path_to_output_folder + '/indel_Number_vs_branch_length_scatter.png', transparent=True)

def plot_number_of_indel_char_vs_bl_reg(indel_df, path_to_output_folder):

    indel_df = indel_df[indel_df['event'].isin(['insertion','deletion'])]
    indel_df = indel_df[['child','tMRCA','ASR','MSA']]

    indel_df['indel_count'] = 1
    grouped_df = indel_df.groupby(['child','ASR','MSA']).agg({'tMRCA': 'mean', 'indel_count': 'sum'}).reset_index()


    plt.figure(figsize=(12, 8))

    g = sns.FacetGrid(grouped_df, col='MSA', hue='ASR', palette='terrain', height=4, aspect=1.2, col_wrap=2)

    # Regression line in each subplot
    g.map(sns.regplot, 'tMRCA', 'indel_count', scatter=False, line_kws={'lw': 1.5}) # , ci=None , lowess=True

    g.set_axis_labels('Branch Length (years)', 'Inserted/Deleted Residues Count')

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(path_to_output_folder + '/indel_Number_vs_branch_length_reg.png', transparent=True)



length_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')
        new_df['ASR'] = ASR_tool
        new_df['MSA'] = MSA_tool
        length_df = pd.concat([length_df, new_df], ignore_index =True)


path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/indel_numb_vs_bl'
plot_number_of_indel_char_vs_bl(length_df, path_to_output_folder)
plot_number_of_indel_char_vs_bl_reg(length_df, path_to_output_folder)
