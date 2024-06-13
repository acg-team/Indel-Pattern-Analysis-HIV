from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
import re

def extract_regions(path_to_msa):
    msa = AlignIO.read(path_to_msa, format='fasta')
    for record in msa._records:
        if re.search(r'K03455', record.id):
            count = 0
            for i in range(len(record.seq)):
                if i == 0:
                        start_signal_peptide=i
                if record.seq[i] != '-' and record.seq[i] != '*':
                    if count == 29:
                        end_signal_peptide=i
                        print('End of env signal peptide', i)
                    elif count == 30:
                        start_c1=i
                    elif count == 128:
                        end_c1=i
                        print('C1 from', start_c1, 'to', end_c1)
                    elif count == 129:
                        start_v1=i
                    elif count == 156:
                        end_v1=i
                        print('V1 from', start_v1, 'to', end_v1)
                    elif count == 157:
                        start_v2=i
                    elif count == 196:
                        end_v2=i
                        print('V2 from', start_v2, 'to', end_v2)
                    elif count == 197:
                        start_c2=i
                    elif count == 293:
                        end_c2=i
                        print('C2 from', start_c2, 'to', end_c2)
                    elif count == 294:
                        start_v3=i
                    elif count == 331:
                        end_v3=i
                        print('V3 from', start_v3, 'to', end_v3)
                    elif count == 332:
                        start_c3=i
                    elif count == 382:
                        end_c3=i
                        print('C3 from', start_c3, 'to', end_c3)
                    elif count == 383:
                        start_v4=i
                    elif count == 418:
                        end_v4=i
                        print('V4 from', start_v4, 'to', end_v4)
                    elif count == 419:
                        start_c4=i
                    elif count == 458:
                        end_c4=i
                        print('C4 from', start_c4, 'to', end_c4)
                    elif count == 459:
                        start_v5=i
                    elif count == 469:
                        end_v5=i
                        print('V5 from', start_v5, 'to', end_v5)
                    elif count == 470:
                        start_c5=i
                    elif count == 510:
                        end_c5=i
                        print('C5 from', start_c5, 'to', end_c5) 
                    elif count == 510:
                        print('End of gp120', i)
                    elif count == 511:
                        start_gp41=i    
                    elif count == 855 or i == len(record.seq)-1:
                        end_gp41=i
                        print('gp41 from', start_gp41, 'to', end_gp41,
                        'or end of alignment', len(record.seq))
                    count += 1
                if i == len(record.seq)-1:
                        end_gp41=i
                        print('gp41 from', start_gp41, 'to', end_gp41,
                        'or end of alignment', len(record.seq))
   
    start_signal_peptide_c1 = end_signal_peptide+1
    end_signal_peptide_c1 = start_c1-1
    start_c1_v1=end_c1+1
    end_c1_v1=start_v1-1
    start_v1_v2=end_v1+1
    end_v1_v2=start_v2-1
    start_v2_c2=end_v2+1
    end_v2_c2=start_c2-1
    start_c2_v3=end_c2+1
    end_c2_v3=start_v3-1
    start_v3_c3=end_v3+1
    end_v3_c3=start_c3-1
    start_c3_v4=end_c3+1
    end_c3_v4=start_v4-1
    start_v4_c4=end_v4+1
    end_v4_c4=start_c4-1
    start_c4_v5=end_c4+1
    end_c4_v5=start_v5-1
    start_v5_c5=end_v5+1
    end_v5_c5=start_c5-1
    start_c5_gp41=end_c5+1
    end_c5_gp41=start_gp41-1
    start_gp41_end=end_gp41+1
    end_gp41_end=len(record.seq)-1

    start_sites = [start_signal_peptide, start_signal_peptide_c1, start_c1, start_c1_v1, start_v1, start_v1_v2, 
                start_v2, start_v2_c2, start_c2, start_c2_v3, start_v3, start_v3_c3, start_c3, 
                start_c3_v4, start_v4, start_v4_c4, start_c4, start_c4_v5, start_v5, start_v5_c5,
                start_c5, start_c5_gp41, start_gp41, start_gp41_end]

    end_sites =  [end_signal_peptide, end_signal_peptide_c1, end_c1, end_c1_v1, end_v1, end_v1_v2, 
                end_v2, end_v2_c2, end_c2, end_c2_v3, end_v3, end_v3_c3, end_c3, 
                end_c3_v4, end_v4, end_v4_c4, end_c4, end_c4_v5, end_v5, end_v5_c5, end_c5, end_c5_gp41, end_gp41, end_gp41_end]  
                 
    return start_sites, end_sites

def define_region_env(start_sites, end_sites, site):
    regions = ['env signal peptide', 'between signal peptide and C1', 'gp120 - C1', 'gp120 - between C1 and V1', 
               'gp120 - V1', 'gp120 - between V1 and V2', 'gp120 - V2', 'gp120 - between V2 and C2', 'gp120 - C2',
               'gp120 - between C2 and V3', 'gp120 - V3', 'gp120 - between V3 and C3', 'gp120 - C3', 'gp120 - between C3 and V4', 
               'gp120 - V4', 'gp120 - between V4 and C4', 
               'gp120 - C4', 'gp120 - between C4 and V5', 'gp120 - V5', 'gp120 - between V5 and C5', 'gp120 - C5', 
               'between C5 and gp41', 'gp41', 'between gp41 and end']
    for i in range(len(regions)):
        if site >= start_sites[i] and site <= end_sites[i]:
            region = regions[i]
    return region

def extract_non_placeholder_chars_from_position(seq, start_position):
    non_placeholder_chars = set('*-')
    chars_after = []
    chars_before = []

    for char in seq[:]:
        if char not in non_placeholder_chars:
            chars_after.append(char)
            if len(chars_after) == 4:
                break

    for char in reversed(seq[:start_position]):
        if char not in non_placeholder_chars:
            chars_before.append(char)
            if len(chars_before) == 3:
                break

    chars_before.reverse()
    seq_of_interest = ''.join(chars_before) + seq[start_position] + ''.join(chars_after)

    return seq_of_interest

def check_if_PNGS(child_seq,parent_seq,site):

    child_extracted_sequence = extract_non_placeholder_chars_from_position(child_seq, site)
    parent_extracted_sequence = extract_non_placeholder_chars_from_position(parent_seq, site)

    pattern = re.compile(r'(?=(N[^P][ST][^P]))')
    match_child = [(match.group(1), match.start()) for match in pattern.finditer(child_extracted_sequence)]
    match_parent = [(match.group(1), match.start()) for match in pattern.finditer(parent_extracted_sequence)]

    if match_child and not match_parent:
        return ['inserted PNGS', match_child[0][1]]
    elif match_parent and not match_child:
        return ['deleted PNGS', match_parent[0][1]]
    elif match_child and match_parent:
        return ['conserved PNGS', match_child[0][1]]
    else:
        return ['NA', 0]


def map_alignment_to_ref(path_to_msa):
    msa = AlignIO.read(path_to_msa, format='fasta')
    map_dict = {}
    for record in msa._records:
        if re.search(r'K03455', record.id):

            pos_in_ref = 0

            for pos_in_msa in range(len(record.seq)):

                if record.seq[pos_in_msa] in ['*','-']:
                    map_dict[pos_in_msa] = pos_in_ref + 0.5
                else:
                    pos_in_ref += 1
                    map_dict[pos_in_msa] = pos_in_ref

    return map_dict

def extract_region_from_msa(msa_filename, start_position, end_position, desired_sequence_ids):
    # Load the MSA from file
    alignment = AlignIO.read(msa_filename, "fasta")

    # Check if the positions are within the alignment length
    if start_position < 0 or end_position >= alignment.get_alignment_length() or start_position > end_position:
        print("Invalid start or end positions.")
        print('alignement length : '+ str(alignment.get_alignment_length()))
        print(start_position)
        print(end_position)
        return None

    # Extract the specified region from the alignment
    region_alignment = alignment[:, start_position:end_position + 1]

    # Filter sequences based on the desired sequence IDs
    filtered_alignment = MultipleSeqAlignment(
        seq for seq in region_alignment if seq.id in desired_sequence_ids
    )

    return filtered_alignment

def compute_background_frequencies(sequences):

    # Create a DataFrame to store the counts of each amino acid
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    counts_df = pd.DataFrame(index=amino_acids, columns=["Count"])
    counts_df["Count"] = 0
    print('start of freq computation')

    # Count the occurrences of each amino acid across all positions
    for sequence in sequences:
        sequence = sequence.upper()
        for aa in amino_acids:
            counts_df.loc[aa, "Count"] = counts_df.loc[aa, "Count"] + int(sequence.seq.count(aa))

    total_count = counts_df["Count"].sum()

    if total_count !=0 :
        counts_df['freq'] = counts_df["Count"] / total_count
    else:
        print('div by 0 !!!')
        print(counts_df)

    return counts_df

def append_count_aa(sequences, previous_count_df,col_name):

    # Create a DataFrame to store the counts of each amino acid
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    counts_df = pd.DataFrame(index=amino_acids, columns=[col_name])
    counts_df[col_name] = 0
    print('start of freq computation')

    # Count the occurrences of each amino acid across all positions
    for sequence in sequences:
        sequence = sequence.upper()
        for aa in amino_acids:

            counts_df.loc[aa, col_name] = counts_df.loc[aa, col_name] + int(sequence.seq.count(aa))

            previous_count_df.loc[aa, col_name] = int(counts_df.loc[aa, col_name]) + int(previous_count_df.loc[aa, col_name])


    return previous_count_df

def read_fasta(file_path):
    sequences = []
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    return sequences