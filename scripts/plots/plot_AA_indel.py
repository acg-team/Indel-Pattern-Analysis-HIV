import os
import matplotlib.pyplot as plt
import seaborn_image as isns
import seaborn as sns
import pandas as pd
import numpy as np

def plot_ins_del_AA(events_df,background_freq_df,path_to_output_folder,reg_to_plot,env_specific):

    if not env_specific :
        env_reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        events_df['region'] = events_df['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)

        background_freq_df = background_freq_df[background_freq_df['region'] == reg_to_plot]
    else :
        background_freq_df = background_freq_df[[reg_to_plot,'amino acid','combination']]

    events_df = events_df[events_df['region'] == reg_to_plot]
    events_df = events_df[events_df['event'].isin(['insertion', 'deletion'])]
    amino_acid_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    events_df['character'] = events_df['character'].str.upper()

    print('test event df')
    print(events_df)

    # Create a new column to represent the combination of event and character
    events_df['event_character_comb'] = events_df['event'] + '_' + events_df['character'] + '_' + events_df['combination']
    count_df = events_df['event_character_comb'].value_counts()
    count_df = count_df.to_frame()

    print('test count df')
    print(count_df)

    count_df.reset_index(inplace=True)

    count_df['amino acid'] = 0
    count_df['event'] = 0
    count_df['combination'] = 0
    count_df['total count'] = 0
    count_df['event - comb'] = 0
    count_df['bckgrnd freq'] = 0

    count_df['MSA freq'] = 0
    count_df['ASR freq'] = 0


    WAG_freq = {'A': 0.0866279, 'R': 0.043972, 'N': 0.0390894, 'D': 0.0570451, 'C': 0.0193078, 'Q': 0.0367281, 'E': 0.0580589, 'G': 0.0832518, 'H': 0.0244313, 'I': 0.048466, 'L': 0.086209, 'K': 0.0620286, 'M': 0.0195027, 'F': 0.0384319, 'P': 0.0457631, 'S': 0.0695179, 'T': 0.0610127, 'W': 0.0143859, 'Y': 0.0352742, 'V': 0.0708956}
    JTT_freq = {'A': 0.077, 'R': 0.051, 'N': 0.043, 'D': 0.052, 'C': 0.020, 'Q': 0.041, 'E': 0.062, 'G': 0.074, 'H': 0.023, 'I': 0.053, 'L': 0.091, 'K': 0.059, 'M': 0.024, 'F': 0.040, 'P': 0.051, 'S': 0.069, 'T': 0.059, 'W': 0.014, 'Y': 0.032, 'V': 0.066}
    HIVb_freq = {'A': 0.060490222, 'R': 0.066039665, 'N': 0.044127815, 'D': 0.042109048, 'C': 0.020075899, 'Q': 0.053606488, 'E': 0.071567447, 'G': 0.072308239, 'H': 0.022293943, 'I': 0.069730629, 'L': 0.098851122, 'K': 0.056968211, 'M': 0.019768318, 'F': 0.028809447, 'P': 0.046025282, 'S': 0.05060433, 'T': 0.053636813, 'W': 0.033011601, 'Y': 0.028350243, 'V': 0.061625237}


    for index, row in count_df.iterrows():
        event_character_comb = row['event_character_comb']
        amino_acid = event_character_comb.split('_')[1]
        event = event_character_comb.split('_')[0]
        combination = event_character_comb.split('_')[2]

        MSA_tool = combination.split(' - ')[0]
        ASR_tool = combination.split(' - ')[1]

        count_df.at[index, 'amino acid'] = amino_acid
        count_df.at[index, 'event'] = event
        count_df.at[index, 'combination'] = combination
        count_df.at[index, 'event - comb'] = event+'_'+combination


        condition1 = (background_freq_df['amino acid'] == amino_acid)
        condition2 = (background_freq_df['combination'] == combination)

        # Use the loc accessor to get the value at the specified position
        if not env_specific:
            result = background_freq_df.loc[condition1 & condition2, 'freq']
        else:
            result = background_freq_df.loc[condition1 & condition2, reg_to_plot]
        
        if not result.empty:
            background_freq = result.values[0]
        else:
            background_freq = 0
        count_df.at[index, 'bckgrnd freq'] = background_freq

        if MSA_tool in ['prank','ProPIP','Historian']:
            if amino_acid in WAG_freq:
                count_df.at[index, 'MSA freq'] = WAG_freq[amino_acid]
        else: # IndelMaP or mafft
            if amino_acid in HIVb_freq:
                count_df.at[index, 'MSA freq'] = HIVb_freq[amino_acid]

        if ASR_tool in ['ArPIP', 'FastML','Historian']:
            if amino_acid in WAG_freq:
                count_df.at[index, 'ASR freq'] = WAG_freq[amino_acid]
        elif ASR_tool in ['grasp']:
            if amino_acid in JTT_freq:
                count_df.at[index, 'ASR freq'] = JTT_freq[amino_acid]
        else: # IndelMAP
            if amino_acid in HIVb_freq:
                count_df.at[index, 'ASR freq'] = HIVb_freq[amino_acid]

    
    for case in count_df['event - comb'].unique():
        sub_df = count_df[count_df['event - comb'] == case]
        total_count = sub_df["count"].sum(axis=0)
        count_df.loc[count_df['event - comb'] == case, 'total count'] = total_count

    count_df['freq'] = count_df["count"] / count_df["total count"]


    combination_order = ['mafft - grasp','mafft - FastML','mafft - IndelMAP','mafft - ArPIP','mafft - Historian',
                'prank - grasp','prank - FastML','prank - IndelMAP','prank - ArPIP','prank - Historian',
                'IndelMAP - grasp','IndelMAP - FastML','IndelMAP - IndelMAP','IndelMAP - ArPIP','IndelMAP - Historian',
                'ProPIP - grasp','ProPIP - FastML','ProPIP - IndelMAP','ProPIP - ArPIP','ProPIP - Historian',
                'Historian - grasp','Historian - FastML','Historian - IndelMAP','Historian - ArPIP','Historian - Historian']

    # Create a bar plot
    plt.figure(figsize=(10, 6))
    g = sns.catplot(x='amino acid', y='freq', hue='event', col='combination', data=count_df,
                kind='bar', order=amino_acid_order,col_order=combination_order, height=4, aspect=1.2, col_wrap=5, palette={'insertion': '#00539C', 'deletion': '#EEA47F'}, alpha=1, dodge=True)

    for ax in g.axes.flatten():
        ax.set_xticks(np.arange(20))
        ax.set_xticklabels(labels=amino_acid_order, rotation=45, ha='right', fontsize=12)

    # Add background frequencies on the plot
    g.map_dataframe(
                    sns.stripplot,
                    x='amino acid',
                    y='bckgrnd freq', 
                    order=amino_acid_order,
                    alpha=0.5,
                    dodge=False,
                    size=6,
                    color = 'red'
                )

    g.map_dataframe(
                    sns.stripplot,
                    x='amino acid',
                    y='MSA freq', 
                    order=amino_acid_order,
                    alpha=0.5,
                    dodge=False,
                    size=6,
                    color = 'blue'
                )

    g.map_dataframe(
                    sns.stripplot,
                    x='amino acid',
                    y='ASR freq', 
                    order=amino_acid_order,
                    alpha=0.5,
                    dodge=False,
                    size=6,
                    color = 'green'
                )

    # Set labels / legends
    g.set_axis_labels(x_var='Amino Acid', y_var='Frequencies')

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)

    # Save the plot
    plt.savefig(path_to_output_folder + '/'+ reg_to_plot + '_AA_proportions_indels.png',transparent=True) # ,transparent=True
    plt.close()

def plot_ins_del_AA_vs_reg(events_df,background_freq_df,path_to_output_folder,env_specific):

    if not env_specific :
        env_reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        events_df['region'] = events_df['region'].apply(lambda x: 'env' if any(s in x for s in env_reg) else x)

        background_freq_df = background_freq_df[background_freq_df['combination'] == 'Historian - Historian']
        events_df = events_df[events_df['combination'] == 'Historian - Historian']
    else :
        env_reg = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        background_freq_df = background_freq_df[background_freq_df['combination'] == 'Historian - Historian']
        events_df = events_df[events_df['region'].isin(env_reg)]
        events_df = events_df[events_df['combination'] == 'Historian - Historian']



    events_df = events_df[events_df['event'].isin(['insertion', 'deletion'])]
    amino_acid_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    events_df['character'] = events_df['character'].str.upper()

    # Create a new column to represent the combination of event and character
    events_df['event_character_reg'] = events_df['event'] + '_' + events_df['character'] + '_' + events_df['region']
    count_df = events_df['event_character_reg'].value_counts()
    count_df = count_df.to_frame()

    print('background_freq_df :')
    for col in background_freq_df.columns:
        print(col)

    count_df.reset_index(inplace=True)

    count_df['amino acid'] = 0
    count_df['event'] = 0
    count_df['region'] = 0
    count_df['total count'] = 0
    count_df['event - reg'] = 0
    count_df['bckgrnd freq'] = 0

    count_df['MSA freq'] = 0
    count_df['ASR freq'] = 0


    WAG_freq = {'A': 0.0866279, 'R': 0.043972, 'N': 0.0390894, 'D': 0.0570451, 'C': 0.0193078, 'Q': 0.0367281, 'E': 0.0580589, 'G': 0.0832518, 'H': 0.0244313, 'I': 0.048466, 'L': 0.086209, 'K': 0.0620286, 'M': 0.0195027, 'F': 0.0384319, 'P': 0.0457631, 'S': 0.0695179, 'T': 0.0610127, 'W': 0.0143859, 'Y': 0.0352742, 'V': 0.0708956}
    JTT_freq = {'A': 0.077, 'R': 0.051, 'N': 0.043, 'D': 0.052, 'C': 0.020, 'Q': 0.041, 'E': 0.062, 'G': 0.074, 'H': 0.023, 'I': 0.053, 'L': 0.091, 'K': 0.059, 'M': 0.024, 'F': 0.040, 'P': 0.051, 'S': 0.069, 'T': 0.059, 'W': 0.014, 'Y': 0.032, 'V': 0.066}
    HIVb_freq = {'A': 0.060490222, 'R': 0.066039665, 'N': 0.044127815, 'D': 0.042109048, 'C': 0.020075899, 'Q': 0.053606488, 'E': 0.071567447, 'G': 0.072308239, 'H': 0.022293943, 'I': 0.069730629, 'L': 0.098851122, 'K': 0.056968211, 'M': 0.019768318, 'F': 0.028809447, 'P': 0.046025282, 'S': 0.05060433, 'T': 0.053636813, 'W': 0.033011601, 'Y': 0.028350243, 'V': 0.061625237}


    for index, row in count_df.iterrows():
        event_character_comb = row['event_character_reg']
        amino_acid = event_character_comb.split('_')[1]
        event = event_character_comb.split('_')[0]
        region = event_character_comb.split('_')[2]

        count_df.at[index, 'amino acid'] = amino_acid
        count_df.at[index, 'event'] = event
        count_df.at[index, 'region'] = region
        count_df.at[index, 'event - reg'] = event+'_'+region


        condition1 = (background_freq_df['amino acid'] == amino_acid)
        condition2 = (background_freq_df['combination'] == 'Historian - Historian')

        # Use the loc accessor to get the value at the specified position
        if not env_specific:
            condition3 = (background_freq_df['region'] == region)
            result = background_freq_df.loc[condition1 & condition2 & condition3, 'freq']
        else:
            result = background_freq_df.loc[condition1 & condition2, region]
        
        if not result.empty:
            background_freq = result.values[0]
        else:
            background_freq = 0
        count_df.at[index, 'bckgrnd freq'] = background_freq

        if amino_acid in WAG_freq:
            count_df.at[index, 'MSA freq'] = WAG_freq[amino_acid]
    
    for case in count_df['event - reg'].unique():
        sub_df = count_df[count_df['event - reg'] == case]
        total_count = sub_df["count"].sum(axis=0)
        count_df.loc[count_df['event - reg'] == case, 'total count'] = total_count

    count_df['freq'] = count_df["count"] / count_df["total count"]
    
    print(count_df)
    print(list(count_df.columns.values))

    if not env_specific:
        reg_order = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
        col_num = 3
    else:
        reg_order = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
        col_num = 4

    # Create a bar plot
    plt.figure(figsize=(10, 6))
    g = sns.catplot(x='amino acid', y='freq', hue='event', col='region', data=count_df,
                kind='bar', order=amino_acid_order,col_order=reg_order, height=4, aspect=1.2, col_wrap= col_num, palette={'insertion': '#00539C', 'deletion': '#EEA47F'}, alpha=1, dodge=True)


    # Add background frequencies on the plot
    g.map_dataframe(
                    sns.stripplot,
                    x='amino acid',
                    y='bckgrnd freq', 
                    order=amino_acid_order,
                    alpha=0.5,
                    dodge=False,
                    size=6,
                    color = 'red'
                )

    g.map_dataframe(
                    sns.stripplot,
                    x='amino acid',
                    y='MSA freq', 
                    order=amino_acid_order,
                    alpha=0.5,
                    dodge=False,
                    size=6,
                    color = 'blue'
                )


    # Set labels / legends
    g.set_axis_labels(x_var='Amino Acid', y_var='Frequencies')
    g.set_xticklabels(labels=amino_acid_order, rotation=45, ha='right', fontsize=12)

    # Set the background of the plot area to be transparent
    g.fig.patch.set_alpha(0)
    
    if env_specific :
        reg = 'env'
    else:
        reg = 'all'

    # Save the plot
    plt.savefig(path_to_output_folder + '/'+reg+'_AA_proportions_indels.png',transparent=True) # ,transparent=True
    plt.close()



genes = ['env','gag','nef','pol','rev','tat','vif','vpr','vpu']
MSA_tools = ['Historian', 'IndelMAP','mafft','prank','ProPIP']
ASR_tools = ['ArPIP', 'FastML', 'grasp', 'Historian', 'IndelMAP']


# Plot plot env specific regions (V1-V5 and gp41)
indel_df = pd.DataFrame()
background_freq_df = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')
        new_df['combination'] = MSA_tool + ' - ' + ASR_tool
        indel_df = pd.concat([indel_df, new_df], ignore_index =True)

        for gene in genes:
            new_freq_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_'+gene+'_AA_freq.csv',sep='\t')

            new_freq_df.rename(columns={"Unnamed: 0": "amino acid"}, inplace=True)
            new_freq_df['combination'] = MSA_tool + ' - ' + ASR_tool
            new_freq_df['region'] = gene

            background_freq_df = pd.concat([background_freq_df, new_freq_df], ignore_index =True)

for reg_to_plot in genes:
    path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/AA_proportions'
    plot_ins_del_AA(indel_df, background_freq_df, path_to_output_folder, reg_to_plot, False)


path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/AA_proportions_vs_reg'
plot_ins_del_AA_vs_reg(indel_df, background_freq_df, path_to_output_folder, False)



# Plot plot env specific regions (V1-V5 and gp41)
indel_df = pd.DataFrame()
env_spec_freq = pd.DataFrame()
for MSA_tool in MSA_tools:
    for ASR_tool in ASR_tools:

        new_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_indel_df.csv',sep='\t')
        new_df['combination'] = MSA_tool + ' - ' + ASR_tool
        indel_df = pd.concat([indel_df, new_df], ignore_index =True)

        new_env_freq_df = pd.read_csv('/cfs/earth/scratch/sepe/TM/results/results/indel_analysis/'+MSA_tool+'_'+ASR_tool+'/'+MSA_tool+'_'+ASR_tool+'_env_Vreg_AA_freq.csv',sep='\t',index_col=0)
        new_env_freq_df.reset_index(inplace=True)
        new_env_freq_df.rename(columns={"index": "amino acid"}, inplace=True)

        new_env_freq_df = new_env_freq_df[['amino acid','env signal peptide freq','gp120 - C1 freq','gp120 - V1 freq','gp120 - V2 freq','gp120 - C2 freq','gp120 - V3 freq','gp120 - C3 freq','gp120 - V4 freq','gp120 - C4 freq','gp120 - V5 freq','gp120 - C5 freq','gp41 freq']]
        new_env_freq_df.rename(columns={'gp120 - V1 freq': 'gp120 - V1', 
                                    'gp120 - V2 freq': 'gp120 - V2',
                                    'gp120 - V3 freq': 'gp120 - V3',
                                    'gp120 - V4 freq': 'gp120 - V4',
                                    'gp120 - V5 freq': 'gp120 - V5',
                                    'gp120 - C1 freq': 'gp120 - C1',
                                    'gp120 - C2 freq': 'gp120 - C2',
                                    'gp120 - C3 freq': 'gp120 - C3',
                                    'gp120 - C4 freq': 'gp120 - C4',
                                    'gp120 - C5 freq': 'gp120 - C5',
                                    'env signal peptide freq': 'env signal peptide',
                                    'gp41 freq': 'gp41',
                                    }, inplace=True)

        new_env_freq_df['combination'] = MSA_tool + ' - ' + ASR_tool

        env_spec_freq = pd.concat([env_spec_freq, new_env_freq_df], ignore_index =True)

for reg_to_plot in ['gp120 - V1','gp120 - V2','gp120 - V3','gp120 - V4','gp120 - V5','gp41']:
    path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/AA_proportions'
    plot_ins_del_AA(indel_df, env_spec_freq, path_to_output_folder, reg_to_plot, True)

path_to_output_folder = '/cfs/earth/scratch/sepe/TM/results/results/plots/AA_proportions_vs_reg'
plot_ins_del_AA_vs_reg(indel_df, env_spec_freq, path_to_output_folder, True)