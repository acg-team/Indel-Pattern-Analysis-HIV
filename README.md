# Comparative Analysis of Phylogenetic Inference Tools in HIV-1 Subtype B

The repository includes the data and the code used for the manuscript Impact of Methodological Choices on Indel Detection: A Comparative Analysis of Phylogenetic Inference Tools in HIV-1 Subtype B, currently under review.

## Manuscript Abstract

Insertions and deletions (indels) play a critical role in the evolutionary dynamics of genomes, yet their accurate detection and interpretation in phylogenetic studies remain challenging. Our study investigates the influence of different multiple sequence alignment (MSA) and ancestral sequence reconstruction (ASR) tools on indel pattern reconstruction, focusing on HIV-1 subtype B. We aim to understand how methodological choices affect the detection of indels, thereby emphasizing the importance of selecting appropriate tools for evolutionary analyses to improve phylogenetic accuracy. We conducted a comparative analysis using five MSA tools (MAFFT, PRANK+F, IndelMaP, ProPIP, Historian) and five ASR tools (GRASP, FastML, IndelMaP, ARPIP, Historian). By examining inferred indel events across all tool combinations, we evaluated their rates, lengths, and positions within the genome, specifically analyzing the env gene and its V1 variable loop. Even though each method tested was able to reconstruct known variable regions in the env gene, our results highlight that the choice of MSA tool significantly impacts indel conservation and interpretation, more so than the choice of ASR tool. This finding underscores the necessity of context-specific MSA tool selection in phylogenetic studies and provides crucial insights for improving the accuracy of indel detection and evolutionary inferences in phylogenetic studies of HIV-1 and other genomes.

## Data
The data folder contains the protein sequences of all nine HIV-1 subtype B genes used in our study, collected from the Los Alamos National Laboratory HIV database. Only sequences originating from Europe were used. Subsets used for certain method combination (as described in the manuscript) are also available here. Newick trees were obtained using IQ-TREE on a mafft alignement of the sequences. 

## Code
The code used in the manuscript is available in the scripts folder, classified between subsets (code used for subsetting), exec (code used to call and run the different phylogenetic tools), IndelAnalysis (code used to extract indel information from the result and compute indel position, rates, lengths) and plots. Each python script was executed using a slurm script, also available.
