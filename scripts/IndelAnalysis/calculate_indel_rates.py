import re
import pandas as pd
from ete3 import PhyloNode
from helper_functions import *
from Bio import AlignIO


def calculate_indel_rates(events_df, path_to_events_all, path_to_nwk, protein):
    """

    Args:
        events_df (DataFrame): output from write_evolutionary_events
        path_to_events_all (string): path to multiple sequence alignment of internal and leave sequences in indelMaP format
        path_to_nwk (string): path to guide tree in newick format
        protein (string): protein which is analyzed

    Returns:
        DataFrame: holds all inferred rates for insertion, deletions and substitutions
    """
    sub_rates_dic = {'env': 0.0057326, 
                'gag': 0.00226134,
                'nef': 0.00256465,
                'pol': 0.000902068,
                'rev': 0.00546282,
                'tat': 0.00557396,
                'vif': 0.00125458,
                'vpr': 0.00122074,
                'vpu': 0.00208927}

    # infer minimum branch length
    tree_string = [line for line in open(path_to_nwk).readlines()]
    bl = re.findall(r":[-+]?(?:\d*\.*\d+)",tree_string[0])
    bl_float = [float(item.split(':')[1]) for item in bl]
    min_bl = min(bl_float)

    msa = AlignIO.read(path_to_events_all, format='fasta')
    tree = PhyloNode(newick=path_to_nwk, alignment=path_to_events_all, format=1)
    if protein == "env":
        start_sites, end_sites = extract_regions(path_to_events_all)

    rates_df = pd.DataFrame(columns=['parent', 'child', 'tMRCA', 'region',
                                    'number of residues in region',
                                    'number of substitutions', 'substitution/AA', 'substitution/AA'+'/year',
                                    'number of deletions', 'deletion/AA', 'deletion/AA'+'/year',
                                    'number of insertions', 'insertion/AA', 'insertion/AA'+'/year'])
    if protein == 'env':
        regs = ['signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
    else:
        regs = [protein]
    
    for node in events_df['child'].unique():
        sub_df = events_df[events_df['child'] == node]
        parent = sub_df['parent'].unique()[0]
        tMRCA =  sub_df['tMRCA'].unique()[0]

        if min_bl >= 1.0:
            min_bl = 0.0

        if tMRCA > (min_bl / sub_rates_dic[protein]) and tMRCA >= 0.001 : # second condition added to avoid inflating rates with really short bl          
            
            seq = tree.search_nodes(name=node)[0].sequence

            # take the length of the parent sequence for rates computation
            parent_seq = tree.search_nodes(name=parent)[0].sequence
            

            for j in range(len(regs)): # 0 to number of regions -1
                if protein == 'env':
                    start_reg = start_sites[j*2]
                    end_reg = end_sites[j*2]
                else:
                    start_reg = 0
                    end_reg = len(seq)-1 #previously end_reg = len(msa)-1
                length_reg = 0

                for site in range(start_reg,end_reg+1):
                    if parent_seq[site] not in ['*','-']:       # take the length of the parent sequence for rates computation
                        length_reg += 1

                if not re.search(r'K03455', node) and length_reg != 0:
                    subs = 0
                    ins = 0
                    dels = 0

                    for row in range(len(sub_df)):
                        if regs[j] in sub_df['region'].iloc[row] and row == 0:
                            if sub_df['event'].iloc[row] == 'substitution':
                                subs += 1
                            if sub_df['event'].iloc[row] == 'insertion':
                                ins += 1
                            if sub_df['event'].iloc[row] == 'deletion':
                                dels += 1
                        elif regs[j] in sub_df['region'].iloc[row]:
                            if sub_df['event'].iloc[row] == 'substitution' and sub_df['event'].iloc[row] != sub_df['event'].iloc[row-1] and sub_df['position w/o placeholders'].iloc[row] != sub_df['position w/o placeholders'].iloc[row-1] + 1:
                                subs += 1
                            if sub_df['event'].iloc[row] == 'insertion' and sub_df['event'].iloc[row] != sub_df['event'].iloc[row-1] and sub_df['position w/o placeholders'].iloc[row] != sub_df['position w/o placeholders'].iloc[row-1] + 1:
                                ins += 1
                            if sub_df['event'].iloc[row] == 'deletion' and sub_df['event'].iloc[row] != sub_df['event'].iloc[row-1] and sub_df['position w/o placeholders'].iloc[row] != sub_df['position w/o placeholders'].iloc[row-1] + 1:
                                dels += 1

                    if protein == 'env':
                        region = regs[j]
                    else:
                        region = protein

                    rates_df.loc[len(rates_df)] = [parent, node, tMRCA, region, length_reg, 
                                                subs, subs/length_reg, (subs/length_reg)/tMRCA,
                                                dels, dels/length_reg, (dels/length_reg)/tMRCA,
                                                ins, ins/length_reg, (ins/length_reg)/tMRCA]
    return rates_df

    