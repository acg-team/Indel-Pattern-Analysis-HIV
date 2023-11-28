from ete3 import PhyloNode
import pandas as pd
import re
from helper_functions import *


def write_evolutionary_events(path_to_events_all, path_to_nwk, sub_rate, protein):
    """
    Args:
        path_to_events_all (string): path to multiple sequence alignment of internal and leave sequences in indelMaP format
        path_to_nwk (string): path to guide tree in newick format
        sub_rate (float): substitution rate inferred by the tree searching algorithm for the specific protein
        protein (string): protein which is analyzed

    Returns:
        DataFrame: holding all inferred evolutionary events
    """

    map_dict = map_alignment_to_ref(path_to_events_all)

    tree = PhyloNode(newick=path_to_nwk, alignment=path_to_events_all, format=1)
    if protein=='env':
        start_sites, end_sites = extract_regions(path_to_events_all)
    events_df = pd.DataFrame(columns=['parent', 'child', 'tMRCA',
                                'event', 'character', 'region', 
                                'site number in alignment', 
                                'position in alignment',
                                'position w/o placeholders','position in ref',
                                'PNGS'])

    for node in tree.traverse('preorder'):
        if not node.is_root() and not re.search(r'K03455', node.name):
            pos_wo_placeholders = 0

            for site in range(len(node.sequence)):
                char = node.sequence[site]
                parent = node.up.name
                char_up = node.up.sequence[site]
                time = node.dist/sub_rate
                if protein=='env':
                    region = define_region_env(start_sites, end_sites, site)
                else:
                    region = protein
                
                site_no = site+1

                if char not in ['*','-']:
                    pos_wo_placeholders += 1

                if char.islower() and char_up == '*':
                    event = 'insertion'

                    # get position in the reference
                    pos_in_ref = map_dict[site]

                    # check if Potential N-linked Glycosylation site (PNGS)
                    PNGS_status = check_if_PNGS(node.sequence,node.up.sequence,site)


                    events_df.loc[len(events_df)] = [parent, node.name, time, 
                                                event, char, region, 
                                                site_no, site, pos_wo_placeholders, pos_in_ref, PNGS_status]

                elif char == '-' and char_up != '-':
                    event = 'deletion' 

                    pos_wo_placeholders += 1

                    # get position in the reference
                    pos_in_ref = map_dict[site]

                    # check if Potential N-linked Glycosylation site (PNGS)
                    PNGS_status = check_if_PNGS(node.sequence,node.up.sequence,site)

                    events_df.loc[len(events_df)] = [parent, node.name, time, 
                                                event, char_up, region, 
                                                site_no, site, pos_wo_placeholders, pos_in_ref, PNGS_status]

                elif char not in ['*', '-'] and char.upper() != char_up.upper() and char_up not in ['*', '-']:
                    event = 'substitution'

                    # get position in the reference
                    pos_in_ref = map_dict[site]

                    # check if Potential N-linked Glycosylation site (PNGS)
                    PNGS_status = check_if_PNGS(node.sequence,node.up.sequence,site)

                    char = char_up+' -> '+char
                    events_df.loc[len(events_df)] = [parent, node.name, time, 
                                                event, char, region, 
                                                site_no, site, pos_wo_placeholders, pos_in_ref, PNGS_status]
    
    return events_df
