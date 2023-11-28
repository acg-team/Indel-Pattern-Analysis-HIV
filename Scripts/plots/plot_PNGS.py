import os
import seaborn_image as isns
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_indel_hotspots(events_df, path_to_output_folder, tools_comb, reg_to_plot, env_specific):
    """Plot the indel hotspot for each region. For gp120 each of the regions is plotted separately and combined.

    Args:
        path_to_events_csv (string): path to csv file with evolutionary events
        path_to_outputfolder (string): path to outputfolder for the created plots

    """

    only_position = events_df[['region','PNGS','position in ref']]


    if not env_specific :
        env_reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        only_position['region'] = only_position['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)

    else:
        regs_list_gp120 = ['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                                'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5']
        rpy_gp120 = only_position.loc[only_position['region'].isin(regs_list_gp120)]
        rpy_gp120.loc[only_position['region'].isin(regs_list_gp120),'region'] = ['gp120']*len(rpy_gp120)
        only_position = pd.concat([only_position,rpy_gp120], axis=0, ignore_index=True)   

    only_position_indel = only_position[only_position['PNGS'].isin(['+','-'])]

    for reg in reg_to_plot: # previously: only_position_indel.region.unique(): 

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

        isns.set_context('paper',fontfamily='Times New Roman')
        isns.set_context("paper", rc={"font.size":8,"axes.titlesize":10,"axes.labelsize":10})
        g = sns.displot(data=only_position_indel[only_position_indel['region']==reg], x='position in ref', hue= 'PNGS', multiple='stack', 
                        discrete=True, 
                        height=1.7, aspect=2,binwidth=2,
                        palette={'+':'#00539C', '-':'#EEA47F'}) #, 'PNGS':'#32a852'

        g.tight_layout()
        for ax in g.axes.flat:
            isns.set_context('paper',fontfamily='Times New Roman')
            sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})

            ax.set_axisbelow(True)
            plt.setp(ax.get_yticklabels(), rotation=90, ha='center', fontsize=8)
            plt.setp(ax.get_xticklabels(), fontsize=8)
            ax.tick_params(axis='y', direction='out', pad=5)
            ax.tick_params(axis=u'both', which=u'both',length=0)
            ax.set_xlabel('region = '+reg)
            ax.xaxis.set_label_coords(0.5, -0.2)
            ax.yaxis.set_label_coords(-0.1,0.5)
            ax.set_xlim(start_region, end_region)
        g.savefig(os.path.join(path_to_output_folder, 'hotspot_'+reg+'_'+tools_comb+'.png'), dpi=400)

        plt.close()


genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu'] # 'env','gag','nef','pol','rev','tat','vif','vpr','vpu'
MSA_tools = ['IndelMAP','mafft','prank','ProPIP'] # 'IndelMAP','mafft','prank','ProPIP'
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'IndelMAP'] # 'ArPIP', 'FastML', 'grasp', 'IndelMAP'

# Indel hotspots vs region (gene), for each tool combination
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        events_df = pd.read_csv('/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')
        path_to_output_folder = '/cfs/earth/scratch/seppemic/TM/results/plots/PNGS_vs_region'
        tools_comb = MSA_tool + '_' + ASR_tool
        reg_to_plot = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']

        plot_indel_hotspots(events_df, path_to_output_folder, tools_comb, reg_to_plot,env_specific=False)

        print('done for combination '+MSA_tool+'-'+ASR_tool)


# Indel hotspots vs region (env - specific), for each tool combination
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        events_df = pd.read_csv('/cfs/earth/scratch/seppemic/TM/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')
        path_to_output_folder = '/cfs/earth/scratch/seppemic/TM/results/plots/PNGS_vs_region'
        tools_comb = MSA_tool + '_' + ASR_tool
        reg_to_plot = ['gp120','env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']

        plot_indel_hotspots(events_df, path_to_output_folder, tools_comb, reg_to_plot,env_specific=True)

        print('done for combination '+MSA_tool+'-'+ASR_tool)

