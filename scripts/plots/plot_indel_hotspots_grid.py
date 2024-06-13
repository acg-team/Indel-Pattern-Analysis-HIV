import os
import seaborn_image as isns
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.ticker as ticker

def plot_hotspots_per_combination(hotspots_df, path_to_output_folder, order, reg, env_specific):
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

    hue_order = ['deletion','insertion']

    isns.set_context('paper',fontfamily='Times New Roman')
    isns.set_context("paper", rc={"font.size":6,"axes.titlesize":6,"axes.labelsize":6})

    bin_edges = np.arange(start_region, end_region + 0.5, 0.5)

    g = sns.displot(data=only_position_indel[only_position_indel['region']==reg], x='position in ref', hue='event', multiple='stack', col="combination", 
                    discrete=False, facet_kws={'sharey': True, 'sharex': True, 'despine': False, 'margin_titles':False,'subplot_kws':{'title':None}}, 
                    height=1.2, aspect=1.6, col_order=(order), kind='hist',
                    palette={'insertion':'#00539C', 'deletion':'#EEA47F'}, col_wrap=5, bins=bin_edges)

    g.tight_layout(h_pad=0.1, w_pad=0.1)
    for ax in g.axes.flat:
        isns.set_context('paper',fontfamily='Times New Roman')
        sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
        ax.set_xlim(start_region, end_region)
        ax.set_ylim(0, 50)
        ax.set_xlabel('')
        spec = ax.get_title().split(' = ')[1]
        ax.set_title(spec)
        ax.xaxis.set_label_coords(0.5, -0.2)
        ax.yaxis.set_label_coords(-0.1,0.5)
        ax.set_axisbelow(True)
        plt.setp(ax.get_yticklabels(), rotation=90, ha='center', fontsize=8)
        plt.setp(ax.get_xticklabels(), fontsize=8)
        ax.tick_params(axis='y', direction='out', pad=5)
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    sns.move_legend(g, "lower center", bbox_to_anchor=(.5, 1), ncol=2, title=None, frameon=False)

    g.fig.patch.set_alpha(0)
    
    g.savefig(os.path.join(path_to_output_folder, 'hotspot_'+reg+'_TEST.png'), dpi=2000,transparent=True) # ,transparent=True

    plt.close()


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['Historian', 'IndelMAP','mafft','prank','ProPIP']
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP']


# Indel lengths vs region (gene), for each tool combination
hotspots_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')

        new_df['combination'] = MSA_tool+'-'+ASR_tool

        hotspots_df = pd.concat([hotspots_df, new_df], ignore_index =True)

path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/hotspots_vs_region'
alphabet = 'protein'
order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian',
            'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian',
            'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian',
            'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian',
            'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']


for reg in genes:
    plot_hotspots_per_combination(hotspots_df, path_to_output_folder, order, reg, env_specific=False)





hotspots_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:
        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')

        new_df['combination'] = MSA_tool+'-'+ASR_tool

        hotspots_df = pd.concat([hotspots_df, new_df], ignore_index =True)

path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/hotspots_vs_region'

order = ['mafft-grasp','mafft-FastML','mafft-IndelMAP','mafft-ArPIP','mafft-Historian',
            'prank-grasp','prank-FastML','prank-IndelMAP','prank-ArPIP','prank-Historian',
            'IndelMAP-grasp','IndelMAP-FastML','IndelMAP-IndelMAP','IndelMAP-ArPIP','IndelMAP-Historian',
            'ProPIP-grasp','ProPIP-FastML','ProPIP-IndelMAP','ProPIP-ArPIP','ProPIP-Historian',
            'Historian-grasp','Historian-FastML','Historian-IndelMAP','Historian-ArPIP','Historian-Historian']
regs_list_env = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']

for reg in regs_list_env:
    plot_hotspots_per_combination(hotspots_df, path_to_output_folder, order, reg, env_specific=True)

