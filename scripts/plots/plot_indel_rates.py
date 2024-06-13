import seaborn as sns 
import seaborn_image as isns 
import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_indel_rates_per_region(rates_df, path_to_output_folder, alphabet, order, tools_comb, env_specific):
    """ Plot indel rates for specified proteins or DNA regions.

    Args:
        path_to_csv_rates (string): path to csv file with the calculated indel rates
        path_to_output_folder (string): path to output folder for the created plots
        alphabet (string): either 'protein' or 'DNA'
        order (list of strings): list with the proteins or DNA regions that should be plotted. 
                                 Example: ['signal peptide', 'gp120', 'gp41'] or ['gp120', 'gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5']

    """
    
    if alphabet == 'protein':
            char = 'AA'
    else:
        char = 'Nt'

    if not env_specific :
        # rates_df = rates_df.where( rates_df.isin(['signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']),'env')
        env_reg = ['signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        rates_df['region'] = rates_df['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)

    rates_per_year = rates_df[['region','substitution/'+char+'/year','insertion/'+char+'/year', 'deletion/'+char+'/year']]
    rpy_melt = pd.melt(rates_per_year, id_vars=['region'], value_vars=['insertion/'+char+'/year', 'deletion/'+char+'/year'],
                    var_name='indel event', value_name='sites/'+char+'/year')

    def split(x):
        return x.split('/')[0]

    regs_list_gp120 = ['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                            'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5']
    rpy_melt['indel event'] = rpy_melt['indel event'].apply(split) 

    # add lines in the dataframe with all rates in gp120 region
    rpy_gp120 = rpy_melt.loc[rpy_melt['region'].isin(regs_list_gp120)]
    rpy_gp120.loc[rpy_melt['region'].isin(regs_list_gp120),'region'] = ['gp120']*len(rpy_gp120)
    rpy_melt = pd.concat([rpy_melt,rpy_gp120], axis=0)

    isns.set_context('paper',fontfamily='Times New Roman')
    sns.set_context("paper", rc={"font.size":5,"axes.titlesize":5,"axes.labelsize":5})  
    g = sns.FacetGrid(rpy_melt, height=1.9, aspect=1.7, despine=False)
    g.map(sns.barplot, x='sites/'+char+'/year', y="region", hue='indel event', order=(order),
        errorbar=('ci', 95), orient='h', linewidth=0.2, errwidth=0.5, errcolor='black', # err_kws={'linewidth': 0.5}, 
                capsize=.1, palette={'insertion':'#00539C', 'deletion':'#EEA47F', 'substitution':'red'},
                data=rpy_melt, width=0.8)

    for ax in g.axes.flat:
        sns.set_context("paper", rc={"font.size": 3, "axes.titlesize": 3, "axes.labelsize": 3})
        ax.set_axisbelow(True)
        ax.grid(True, which='both', axis='x', zorder=-1, linestyle='dashed', linewidth=0.55)
        ax.tick_params(axis='y', direction='out', pad=2)
        ax.tick_params(axis='x', which='both', direction='out', pad=2, length=4)  # Add more ticks to the x-axis
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6))  # Adjust the number of major ticks on the x-axis
        plt.setp(ax.get_yticklabels(), ha='right', fontsize=3)
        plt.setp(ax.get_xticklabels(), fontsize=3)  # Change the font size of x-axis labels

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    plt.legend(fontsize=6)
    g.savefig(os.path.join(path_to_output_folder,'indelrate_'+tools_comb+'_'+'_'.join(order)+'.png'), dpi=400, transparent=True) # , transparent=True
    plt.close()

def plot_indel_rates_per_combination(rates_df, path_to_output_folder, alphabet, order, reg, env_specific):
    """ Plot indel rates for specified proteins or DNA regions.

    Args:
        path_to_csv_rates (string): path to csv file with the calculated indel rates
        path_to_output_folder (string): path to output folder for the created plots
        alphabet (string): either 'protein' or 'DNA'
        order (list of strings): list with the proteins or DNA regions that should be plotted. 
                                 Example: ['signal peptide', 'gp120', 'gp41'] or ['gp120', 'gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5']

    """
    
    if alphabet == 'protein':
            char = 'AA'
    else:
        char = 'Nt'

    if not env_specific : # replace all regions specific to env gene by 'env'
        env_reg = ['signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        rates_df['region'] = rates_df['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)
    
        mask = rates_df['region'] == reg
        rates_df = rates_df[mask]
        print('Not env specific, for gene : ' + reg)
        print(rates_df.loc[1])
    else:
        regs_list_env = ['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                            'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5']
        rates_gp120 = rates_df.loc[rates_df['region'].isin(regs_list_env)]
        rates_gp120.loc[rates_df['region'].isin(regs_list_env),'region'] = ['gp120']*len(rates_gp120)
        rates_df = pd.concat([rates_df,rates_gp120], axis=0)

        mask = rates_df['region'] == reg
        rates_df = rates_df[mask]

    rates_per_year = rates_df[['region','substitution/'+char+'/year','insertion/'+char+'/year', 'deletion/'+char+'/year','combination']]
    rpy_melt = pd.melt(rates_per_year, id_vars=['combination'], value_vars=['insertion/'+char+'/year', 'deletion/'+char+'/year'],
                    var_name='indel event', value_name='sites/'+char+'/year')

    def split(x):
        return x.split('/')[0]

    rpy_melt['indel event'] = rpy_melt['indel event'].apply(split) 


    rpy_melt.loc[len(rpy_melt)] = ['blankline1', 'deletion', 0]
    rpy_melt.loc[len(rpy_melt)] = ['blankline2', 'deletion', 0]
    rpy_melt.loc[len(rpy_melt)] = ['blankline3', 'deletion', 0]
    rpy_melt.loc[len(rpy_melt)] = ['blankline4', 'deletion', 0]

    print('print rpy_melt content:')
    print(rpy_melt)

    isns.set_context('paper',fontfamily='Times New Roman')
    sns.set_context("paper", rc={"font.size":4,"axes.titlesize":4,"axes.labelsize":4,"axes.linewidth": 0.4, "xtick.major.width": 0.4, "ytick.major.width": 0.4})  
    g = sns.FacetGrid(rpy_melt, height=1.8, aspect=1.7, despine=False)
    g.map(sns.barplot, x='sites/'+char+'/year', y="combination", hue='indel event', order=order,
        errorbar=('ci', 95), orient='h', linewidth=0.2, errwidth=0.5, errcolor='black', #err_kws={'linewidth': 0.5}, 
                capsize=.1, palette={'insertion':'#00539C', 'deletion':'#EEA47F', 'substitution':'red'},
                data=rpy_melt, width=0.8)
    g.set_axis_labels('event/'+char+'/year')

    for ax in g.axes.flat:
        sns.set_context("paper", rc={"font.size": 3, "axes.titlesize": 2, "axes.labelsize": 2})
        ax.set_axisbelow(True)
        ax.grid(True, which='both', axis='x', zorder=-1, linestyle='dashed', linewidth=0.2)
        ax.tick_params(axis='y', direction='out', pad=2, length=1)
        ax.tick_params(axis='x', which='both', direction='out', pad=2, length=0)  # Add more ticks to the x-axis
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6))  # Adjust the number of major ticks on the x-axis
        plt.setp(ax.get_yticklabels(), ha='right', fontsize=3)
        plt.setp(ax.get_xticklabels(), fontsize=3)  # Change the font size of x-axis labels
        ax.set_ylabel('')

        for i, label in enumerate(ax.get_yticklabels()):
            if i % 6 == 5:
                label.set_visible(False)  # Hide the label every 5th tick
                ax.yaxis.get_major_ticks()[i].set_visible(False)  # Hide the tick every 5th tick


    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    plt.legend(fontsize=2)
    g.savefig(os.path.join(path_to_output_folder,'indelrate_'+reg+'_combination_comp.png'), dpi=2000, transparent=True) # , transparent=True
    plt.close()



genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['Historian', 'IndelMAP','mafft','prank','ProPIP']
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP']

# # # Indel rates vs region (gene), for each tool combination
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        #concatenate all rates files in one for each tool combination
        rates_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+genes[0]+'_indel_df_rates.csv',sep='\t')
        for i in range(1,len(genes)):
            new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+genes[i]+'_indel_df_rates.csv',sep='\t')
            rates_df = pd.concat([rates_df, new_df], ignore_index =True)
            print(rates_df.values[-1:])
        path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/rates_vs_region'
        alphabet = 'protein'
        order = genes
        tools_comb = MSA_tool + '_' + ASR_tool
        plot_indel_rates_per_region(rates_df, path_to_output_folder, alphabet, order, tools_comb, env_specific = False)

        print('done for combination '+MSA_tool+'-'+ASR_tool)

# # # Indel rates vs region (specific in env), for each tool combination
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        rates_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_indel_df_rates.csv',sep='\t')
        path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/rates_vs_region'
        alphabet = 'protein'
        order = ['gp120','signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        tools_comb = MSA_tool + '_' + ASR_tool
        plot_indel_rates_per_region(rates_df, path_to_output_folder, alphabet, order, tools_comb, env_specific = True)

# Indel rates vs tool combination, for each region
for gene in genes :
    rates_df = pd.DataFrame()
    for MSA_tool in MSA_tools:
        for ASR_tool in ASR_tools:
            new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_indel_df_rates.csv',sep='\t')
            new_df['combination'] = MSA_tool+'-'+ASR_tool
            rates_df = pd.concat([rates_df, new_df], ignore_index =True)

    path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/rates_vs_combination'
    alphabet = 'protein'
    order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian','blankline1',
                'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian','blankline2',
                'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian','blankline3',
                'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian','blankline4',
                'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']
    plot_indel_rates_per_combination(rates_df, path_to_output_folder, alphabet, order, gene, env_specific=False)

# Indel rates vs tool combination, for env regions only
env_reg = ['gp120','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
gene = 'env'
for reg in env_reg:
    rates_df = pd.DataFrame()
    for MSA_tool in MSA_tools:
        for ASR_tool in ASR_tools:
            new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_indel_df_rates.csv',sep='\t')

            new_df['combination'] = MSA_tool+'-'+ASR_tool

            rates_df = pd.concat([rates_df, new_df], ignore_index =True)

    path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/rates_vs_combination'
    alphabet = 'protein'
    order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian','blankline1',
                'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian','blankline2',
                'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian','blankline3',
                'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian','blankline4',
                'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']

    plot_indel_rates_per_combination(rates_df, path_to_output_folder, alphabet, order, reg, env_specific=True)
