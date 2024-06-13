import seaborn as sns 
import seaborn_image as isns 
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib.ticker as ticker

def plot_indel_length(length_df, path_to_output_folder,alphabet,order,max_length, tool_comb, env_specific):
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
        length_df['region'] = length_df['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)
        only_length = length_df[['region','event','length']]
        col_num = 3

    else:

        only_length = length_df[['region','event','length']]
        regs_list_gp120 = ['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                                'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5']
        rpy_gp120 = only_length.loc[only_length['region'].isin(regs_list_gp120)]
        rpy_gp120.loc[only_length['region'].isin(regs_list_gp120),'region'] = ['gp120']*len(rpy_gp120)
        only_length = pd.concat([only_length,rpy_gp120], axis=0, ignore_index=True)

        regs_list_env = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        rpy_env = only_length.loc[only_length['region'].isin(regs_list_env)]
        rpy_env.loc[only_length['region'].isin(regs_list_env),'region'] = ['env']*len(rpy_env)
        only_length = pd.concat([only_length,rpy_env], axis=0, ignore_index=True)

        col_num = 5


    mean_length = only_length.groupby(['region', 'event']).length.agg(['mean', 'median'])

    isns.set_context('paper',fontfamily='Times New Roman')
    isns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
    only_length_indel = only_length[only_length['event'].isin(['insertion','deletion'])]
    g = sns.displot(data=only_length_indel, x='length', hue='event', multiple='stack', col="region", 
                    discrete=True, facet_kws={'sharey': False, 'sharex': False, 'despine': False, 'margin_titles':False,'subplot_kws':{'title':None}}, 
                    height=1.7, aspect=1.2, col_order=(order), kind='hist',
                    palette={'insertion':'#00539C', 'deletion':'#EEA47F'}, col_wrap=col_num)

    g.tight_layout(h_pad=0.1, w_pad=0.1)
    for ax in g.axes.flat:
        isns.set_context('paper',fontfamily='Times New Roman')
        sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
        ax.set_xlim(0,max_length)
        if alphabet == 'protein':
            ax.set_xticks(range(2,max_length,2))
        else:
            ax.set_xticks(range(0,max_length,3))        
        spec = ax.get_title().split(' = ')[1]
        d = mean_length.loc[spec, :]
        if (d.index == 'insertion').any():
            ax.axvline(x=d.loc['insertion']['median'], c='#00539C', ls='-', lw=.5)
        if (d.index == 'deletion').any():
            ax.axvline(x=d.loc['deletion']['median'], c='#EEA47F', ls='--', lw=.5)
        ax.set_axisbelow(True)
        plt.setp(ax.get_yticklabels(), rotation=90, ha='center', fontsize=8)
        plt.setp(ax.get_xticklabels(), fontsize=8)
        ax.tick_params(axis='y', direction='out', pad=5)
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax.xaxis.set_label_coords(0.5, -0.2)
        ax.yaxis.set_label_coords(-0.18,0.5)
        ax.set_ylim(0, 15000)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    g.savefig(os.path.join(path_to_output_folder,'indel_length_'+tool_comb+'_'+'_'.join(order)+'.png'), dpi=400,transparent=True) # ,transparent=True
    plt.close()

def plot_indel_length_per_combination(length_df, path_to_output_folder,alphabet,order,max_length, reg, env_specific):
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
        length_df['region'] = length_df['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)
        only_length = length_df[['combination','region','event','length']]

    else:

        only_length = length_df[['combination','region','event','length']]
        regs_list_gp120 = ['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                                'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5']
        rpy_gp120 = only_length.loc[only_length['region'].isin(regs_list_gp120)]
        rpy_gp120.loc[only_length['region'].isin(regs_list_gp120),'region'] = ['gp120']*len(rpy_gp120)
        only_length = pd.concat([only_length,rpy_gp120], axis=0, ignore_index=True)

        regs_list_env = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        rpy_env = only_length.loc[only_length['region'].isin(regs_list_env)]
        rpy_env.loc[only_length['region'].isin(regs_list_env),'region'] = ['env']*len(rpy_env)
        only_length = pd.concat([only_length,rpy_env], axis=0, ignore_index=True)

    mask = only_length['region'] == reg
    only_length = only_length[mask]

    mean_length = only_length.groupby(['combination','event']).length.agg(['mean', 'median'])
    
    isns.set_context('paper',fontfamily='Times New Roman')
    isns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
    only_length_indel = only_length[only_length['event'].isin(['insertion','deletion'])]
    g = sns.displot(data=only_length_indel, x='length', hue='event', multiple='stack', col="combination", 
                    discrete=True, facet_kws={'sharey': True, 'sharex': True, 'despine': False, 'margin_titles':False,'subplot_kws':{'title':None}}, 
                    height=1.7, aspect=1.2, col_order=(order), kind='hist',
                    palette={'insertion':'#00539C', 'deletion':'#EEA47F'}, col_wrap=5, legend=False)


    g.tight_layout(h_pad=0.1, w_pad=0.1)
    for ax in g.axes.flat:
        isns.set_context('paper',fontfamily='Times New Roman')
        sns.set_context("paper", rc={"font.size":8,"axes.titlesize":10,"axes.labelsize":8})
        ax.set_xlim(0,max_length)
        if alphabet == 'protein':
            ax.set_xticks(range(2,max_length,2))
        else:
            ax.set_xticks(range(0,max_length,3))        
        spec = ax.get_title().split(' = ')[1]
        ax.set_title(spec)
        ax.set_xlabel('')
        ax.set_ylabel('')
        d = mean_length.loc[spec, :]
        if (d.index == 'insertion').any():
            ax.axvline(x=d.loc['insertion']['median'], c='#00539C', ls='-', lw=.5)
        if (d.index == 'deletion').any():
            ax.axvline(x=d.loc['deletion']['median'], c='#EEA47F', ls='--', lw=.5)
        ax.set_axisbelow(True)
        plt.setp(ax.get_yticklabels(), rotation=90, ha='center', fontsize=8)
        plt.setp(ax.get_xticklabels(), fontsize=8)
        ax.tick_params(axis='y', direction='out', pad=5)
        ax.tick_params(axis=u'both', which=u'both',length=0)

        ax.set_ylim(0, 15000)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Set global X and Y axis labels
    g.set_axis_labels("", "")  # Temporarily set them to empty strings

    # Add global X and Y axis labels using fig.text
    g.fig.text(0.5, 0, 'Indel Length', ha='center', fontsize=15)
    g.fig.text(0, 0.5, 'Count', va='center', rotation='vertical', fontsize=15)


    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    g.savefig(os.path.join(path_to_output_folder,'indel_length_'+reg+'_combination_comp.png'), dpi=10000,transparent=True) # ,transparent=True
    plt.close()


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['Historian', 'IndelMAP','mafft','prank','ProPIP']
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP']


# Indel lengths vs region (gene), for each tool combination
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        length_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv',sep='\t')
        path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/lengths_vs_region'
        alphabet = 'protein'
        order = genes
        max_length = 15
        tool_comb = MSA_tool + '_' + ASR_tool
        plot_indel_length(length_df, path_to_output_folder,alphabet,order,max_length, tool_comb, env_specific = False)


# Indel lengths vs region (specific to env), for each tool combination
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        length_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv',sep='\t')
        path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/lengths_vs_region'
        alphabet = 'protein'
        order = ['gp120','env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        tools_comb = MSA_tool + '_' + ASR_tool
        max_length = 15
        plot_indel_length(length_df, path_to_output_folder, alphabet, order,max_length, tools_comb, env_specific = True)


# Indel lengths vs tool combination, for each region
lengths_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv',sep='\t')

        new_df['combination'] = MSA_tool+'-'+ASR_tool

        lengths_df = pd.concat([lengths_df, new_df], ignore_index =True)

path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/lengths_vs_combination'
alphabet = 'protein'
order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian',
            'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian',
            'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian',
            'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian',
            'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']
max_length = 15

for reg in genes:
    plot_indel_length_per_combination(lengths_df, path_to_output_folder,alphabet,order,max_length, reg, env_specific=False)


# Indel lengths vs tool combination, for env regions only
env_reg = ['gp120','env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
gene = 'env'
for reg in env_reg:
    lengths_df = pd.DataFrame()
    for MSA_tool in MSA_tools:
        for ASR_tool in ASR_tools:
            new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df_lengths.csv',sep='\t')
            new_df['combination'] = MSA_tool+'-'+ASR_tool
            lengths_df = pd.concat([lengths_df, new_df], ignore_index =True)

    path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/lengths_vs_combination'
    alphabet = 'protein'
    order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian',
                'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian',
                'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian',
                'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian',
                'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']
    max_length = 15

    plot_indel_length_per_combination(lengths_df, path_to_output_folder,alphabet,order,max_length, reg, env_specific=True)

