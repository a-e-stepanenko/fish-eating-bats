library(dplyr)
library(readr)
library(stringr)

transcripts = read.csv("Documents/loss/selection_toga2/selected_transcripts.txt", sep = "\t", header = TRUE)
#View(transcripts)

filtered_transcripts <- transcripts %>%
  group_by(Gene) %>%
  filter(FI_Count >= 0.8 * max(FI_Count), CDS_Length >= 0.9 * max(CDS_Length)) %>%
  filter(FI_Count == max(FI_Count)) %>%
  filter(Mean_Across == max(Mean_Across))
#View(filtered_transcripts)
nrow(filtered_transcripts)
nrow(filtered_transcripts[!duplicated(filtered_transcripts$Gene),])


#dupl_transcripts <- filtered_transcripts[filtered_transcripts$Gene %in% filtered_transcripts$Gene[duplicated(filtered_transcripts$Gene)],]
#View(dupl_transcripts)

filtered_transcripts <- filtered_transcripts[!duplicated(filtered_transcripts$Gene),]
nrow(filtered_transcripts)

gene_loss_full = read.csv("Documents/loss/toga2_loss/lost_genes_new_wo_anc_pl_fin.txt", sep="\t", row.names = 1, header = T) %>% as.data.frame()
genes_for_selection_screen = gene_loss_full[gene_loss_full$HLmyoViv5 %in% c("FI", "I") | gene_loss_full$HLnocLep2 %in% c("FI", "I"),]
genes_for_selection_screen = strsplit(rownames(genes_for_selection_screen),"#") %>% lapply(function(x) x[1]) %>% unlist

filtered_transcripts <- filtered_transcripts[filtered_transcripts$Gene %in% genes_for_selection_screen,]
nrow(filtered_transcripts)

write_delim(filtered_transcripts , "Documents/loss/selection_toga2/filtered_transcripts_full.txt", delim = '\t') 

# FI, I, UL
# one2one, many2one?