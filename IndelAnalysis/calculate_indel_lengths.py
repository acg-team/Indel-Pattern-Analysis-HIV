import re
import os
import pandas as pd

def calculate_indel_lengths(events_df):
    """

    Args:
        events_df (DataFrame): output from write_evolutionary_events

    Returns:
        DataFrame: holds indel lengths per branch
    """
    min_bl = min(events_df['tMRCA']) # time to most recent common ancestor (time = node.dist/sub_rate -> in years)
    length_df = pd.DataFrame(columns=['parent', 'child', 'tMRCA', 'region','event', 'length'])
    length = 1

    for row in range(1,len(events_df)):
        node = events_df['child'].iloc[row]
        prev_node = events_df['child'].iloc[row-1]
        region = events_df['region'].iloc[row]
        prev_region =  events_df['region'].iloc[row-1]
        prev_parent = events_df['parent'].iloc[row-1]
        prev_tMRCA = events_df['tMRCA'].iloc[row-1] # previously events_df['time in years'].iloc[row-1]
        event = events_df['event'].iloc[row]
        prev_event = events_df['event'].iloc[row-1]
        pos = events_df['position w/o placeholders'].iloc[row]
        prev_pos = events_df['position w/o placeholders'].iloc[row-1]

        if event == prev_event and pos == prev_pos+1 and node == prev_node and region == prev_region:
            length += 1
            pass
        else:
            if prev_tMRCA > min_bl: 
                length_df.loc[len(length_df)] = [prev_parent, prev_node, prev_tMRCA, prev_region, prev_event, length]
                length = 1
            else:
                length = 1
                
    return length_df