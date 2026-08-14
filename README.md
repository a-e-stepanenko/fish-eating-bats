# fish-eating-bats
This repository contains code and scripts used in the paper on the comparative genomics of fish-eating bats (Stepanenko et al. link-to-add).

**make_gene_loss_status_table_with_genenames.py** & **make_gene_dupl_status_table_with_genenames.py** 
Extract gene states and orthology classifications based on TOGA2 (https://github.com/hillerlab/TOGA2) annotations. The scripts were run as: 
make_gene_*_status_table_with_genenames.py full_dirs_sp_list.txt > genes.txt

**geme_losses_duplications.R** 
Using the extracted tables from the above scripts and the phylogenetic tree, this script filters for genes lost or duplicated in the fishing bats whereas intact or not duplicated in the sister species and most of the non-fishing bats in the tree. See details in the methods of our paper.

**filter_for_transcripts.R**
This script takes the output table from transcript_selection_table_TOGA2.py and filters for the longest and most complete transcripts for each gene. Then the selected transcripts were  processed by sleasy (https://github.com/casparbein/sleasy) to run Hyphy aBSREL (and RELAX for some of the genes, see details in the paper). Detailed settings of sleasy runs are in the script **sleasy_commands**.

**selection_screen.R**
This script?

**mol_convergence.R** 
These are the codes for running RERconverge (https://github.com/nclark-lab/RERconverge) using the same sleasy-generated transcript alignments ($transcript.manual.fa) for HyPhy runs and the phylogenetic tree. 

