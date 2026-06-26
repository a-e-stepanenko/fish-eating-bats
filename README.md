# fish-eating-bats
This repository contains the code I used in the paper on the comparative genomics of fish-eating bats. 

Tables with the gene states and orthology classifications were generated using the following commands:

./make_gene_loss_status_table_with_genenames.py full_dirs_list.txt > lost_genes.txt

./make_gene_dupl_status_table_with_genenames.py full_dirs_list.txt > dupl_genes.txt

full_dirs_list.txt contains full paths to directories for TOGA2 outputs for all species in the tree.

The outputs contain tables with gene and their states in all genomes from the list (lost_genes.txt), and gene orthology states relative to humann genes (dupl_genes.txt).  

geme_losses_duplications.R takes the tables from the previous step, the table with bat assembly names, scientific names, and diets, and the tree, and filters the lost and duplicated genes according to the methods.

filter_for_transcripts.R filters the longest and most complete transcripts from all transcripts to then run the selection screen pipeline.

mol_convergence.R needs a tree and a folder with the alignments of the genes across all species that have these genes, and performs molecular convergence analysis.
