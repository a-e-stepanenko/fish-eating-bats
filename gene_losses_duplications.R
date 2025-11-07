library(ggplot2); theme_set(theme_bw())
library(ggVennDiagram)
library(gprofiler2)
library(ggpubr)
library(ape)
library(phytools)
library(ggtree)
library(tidytree)
library(cowplot)
library(tidyverse)
library(reshape)
library(phangorn)
library(stringr)
library(kableExtra)
library(stringr)
library(egg)
library(org.Hs.eg.db)
library(ggnewscale)
library(treeio)
library(ggvenn)
library(wrapr)

## tree with diets
setwd("/home/ana/Documents/loss")
tree <- read.tree("tree21fam0.tree") # Bat1K tree with all 103 bats in Newick format
table_bats <- read.csv("Yangochiroptera77.csv") # table with assembly names, scientific bat names, and other information about assemblies
synonyms_df <- data.frame(
  old = c(table_bats$Assembly),
  new = c(table_bats$ScientificName)
)
rownames(table_bats) = table_bats$Assembly
tree <- drop.tip(tree, tree$tip.label[tree$tip.label %notin% c(rownames(table_bats), "hg38")]) # filtering out non-Yangochiropteran bats
write.tree(tree, file='tree_for_selection_screen.tree') # saving tree for the selection screen
tree <- read.tree("Ariadna_21Famout_Yango77_noNA_nodelab.tre") # reading in a subset tree
tree$edge.length <- rep(1, length(tree$edge.length)) # setting all edge lengths to 1
tree$label = rep("background", nrow(tree$edge)+1) # labeling
tree$label[MRCA(tree, "HLmyoDau2")] = tree$label[MRCA(tree, "HLmyoRic2")] = "midground"
tree$label[MRCA(tree, "HLnocLep2")] = tree$label[MRCA(tree, "HLmyoViv5")] = "focus\nspecies"
table_bats = table_bats[, c("Arthropods", "Blood", "Terrestrial.vertebrates", "Fish", "Leaves.and.flower.pieces", "Pollen.and.nectar", "Fruit")]
table_bats = as.data.frame(table_bats)
diets = c("Arthropods", "Blood", "Terrestrial.vertebrates", 
          "Fish", "Leaves.and.flower.pieces", "Pollen.and.nectar", "Fruit", "Seed")
diet.states = c("Absent" = 0, "Complementary" = 1, "Predominant" = 2, "Strict" = 3)
state.color = c("#d2dae1", "#6baed6", "#08519c", "black")

renamed_tree <- rename_taxa(tree, data = synonyms_df, key = "old", value = "new")
rownames(table_bats) = synonyms_df[synonyms_df$old %in% rownames(table_bats),]$new

tree_for_plot = ggplot(renamed_tree, aes(color=tree$label)) + scale_color_manual(values=c("focus\nspecies"="#d1182a", "midground"="#d15500", "background"="#555555")) + 
  geom_tree() + theme_tree() + geom_tiplab(align=TRUE, fontface = "italic", size=7) # plotting a tree amd a heatmap with diets of fish-eating bats
gheatmap(tree_for_plot, table_bats, offset=18, width = 0.5, color = "black",, colnames_angle=90, hjust = 1) +
  scale_fill_continuous(name="Diet code", breaks = diet.states, high = "#132B43", low = "#56B1F7") +
  theme(axis.text.x=element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=15)) + ggtree::vexpand(0, -1) + vexpand(-1)  + vexpand(-1)

## Gene losses TOGA2

gene_loss_full = read.csv("toga2_loss/lost_genes_16072025.txt", sep="\t", row.names = 1, header = T) %>% as.data.frame()
gene_loss_full <- gene_loss_full[, !(colnames(gene_loss_full) %in% c("HLmyoViv2", "HLmyoViv3"))]
gene_loss_full = gene_loss_full[rowMeans(gene_loss_full == "N") != 1,]
genes_for_selection_screen = gene_loss_full[gene_loss_full$HLmyoViv5 %in% c("FI", "I", "UL") & gene_loss_full$HLnocLep2 %in% c("FI", "I", "UL"),]
genes_for_selection_screen = strsplit(rownames(genes_for_selection_screen),"#") %>% lapply(function(x) x[1]) %>% unlist
myoviv_L_intersection <- list(HLmyoViv1 = rownames(gene_loss_full[gene_loss_full$HLmyoViv1 == "L",]),
                              HLmyoViv5 = rownames(gene_loss_full[gene_loss_full$HLmyoViv5 == "L",]))

myoviv1_L_uniq_status_myoviv5 <- gene_loss_full[(gene_loss_full$HLmyoViv1 == "L") & (gene_loss_full$HLmyoViv5 != "L"),]
myoviv5_L_uniq_status_myoviv1 <- gene_loss_full[(gene_loss_full$HLmyoViv5 == "L") & (gene_loss_full$HLmyoViv1 != "L"),]

myoviv_L_UL_intersection <- list(HLmyoViv1 = rownames(gene_loss_full[gene_loss_full$HLmyoViv1 %in% c("L", "UL"),]),
                                 HLmyoViv5 = rownames(gene_loss_full[gene_loss_full$HLmyoViv5 %in% c("L", "UL"),]))
myoviv1_L_UL_uniq_status_myoviv5 <- gene_loss_full[(gene_loss_full$HLmyoViv1 %in% c("L", "UL")) & (gene_loss_full$HLmyoViv5 %notin% c("L", "UL")),]
myoviv5_L_UL_uniq_status_myoviv1 <- gene_loss_full[(gene_loss_full$HLmyoViv5 %in% c("L", "UL")) & (gene_loss_full$HLmyoViv1 %notin% c("L", "UL")),]

noclep_L_intersection <- list(HLnocLep1 = rownames(gene_loss_full[gene_loss_full$HLnocLep1 == "L",]),
                              HLnocLep2 = rownames(gene_loss_full[gene_loss_full$HLnocLep2 == "L",]))
noclep1_L_uniq_status_noclep2 <- gene_loss_full[(gene_loss_full$HLnocLep1 == "L") & (gene_loss_full$HLnocLep2 != "L"),]
noclep2_L_uniq_status_noclep1 <- gene_loss_full[(gene_loss_full$HLnocLep2 == "L") & (gene_loss_full$HLnocLep1 != "L"),]

noclep_L_UL_intersection <- list(HLnocLep1 = rownames(gene_loss_full[gene_loss_full$HLnocLep1 %in% c("L", "UL"),]),
                                 HLnocLep2 = rownames(gene_loss_full[gene_loss_full$HLnocLep2 %in% c("L", "UL"),]))
noclep1_L_UL_uniq_status_noclep2 <- gene_loss_full[(gene_loss_full$HLnocLep1 %in% c("L", "UL")) & (gene_loss_full$HLnocLep2 %notin% c("L", "UL")),]
noclep2_L_UL_uniq_status_noclep1 <- gene_loss_full[(gene_loss_full$HLnocLep2 %in% c("L", "UL")) & (gene_loss_full$HLnocLep1 %notin% c("L", "UL")),]

# 2D Venn diagram
p1 <- ggVennDiagram(myoviv_L_intersection) + ggtitle("Myotis vivesi L") + coord_flip() + scale_x_continuous(expand = expansion(mult = .3))
p1_1 <- ggplot(myoviv1_L_uniq_status_myoviv5, aes(x=factor(HLmyoViv5)))+geom_bar(stat="count")+theme_minimal()+ggtitle("L HLmyoViv1 only")+xlab("state in HLmyoViv5")
p1_2 <- ggplot(myoviv5_L_uniq_status_myoviv1, aes(x=factor(HLmyoViv1)))+geom_bar(stat="count")+theme_minimal()+ggtitle("L HLmyoViv5 only")+xlab("state in HLmyoViv1")
p2 <- ggvenn(myoviv_L_UL_intersection, fill_color=c("darkblue", "darkgreen"), show_percentage = F, set_name_size = 10, text_size = 10)# + ggtitle("Myotis vivesi L and UL") + theme(plot.title = element_text(size = 20))
# + coord_flip() + scale_x_continuous(expand = expansion(mult = .3)) 
p2_1 <- ggplot(myoviv1_L_UL_uniq_status_myoviv5, aes(x=factor(HLmyoViv5)))+geom_bar(stat="count")+theme_minimal()+ggtitle("L+UL HLmyoViv1 only")+xlab("state in HLmyoViv5")+theme(axis.text=element_text(size=14), axis.title=element_text(size=14), plot.title = element_text(size = 14))
p2_2 <- ggplot(myoviv5_L_UL_uniq_status_myoviv1, aes(x=factor(HLmyoViv1)))+geom_bar(stat="count",)+theme_minimal()+ggtitle("L+UL HLmyoViv5 only")+xlab("state in HLmyoViv1")+theme(axis.text=element_text(size=14), axis.title=element_text(size=14), plot.title = element_text(size = 14))
p3 <- ggVennDiagram(noclep_L_intersection) + ggtitle("Noctilio leporinus L") + coord_flip() + scale_x_continuous(expand = expansion(mult = .3))
p3_1 <- ggplot(noclep1_L_uniq_status_noclep2, aes(x=factor(HLnocLep2)))+geom_bar(stat="count")+theme_minimal()+ggtitle("L HLnocLep1 only")+xlab("state in HLnocLep2")
p3_2 <- ggplot(noclep2_L_uniq_status_noclep1, aes(x=factor(HLnocLep1)))+geom_bar(stat="count")+theme_minimal()+ggtitle("L HLnocLep2 only")+xlab("state in HLnocLep1")
p4 <- ggvenn(noclep_L_UL_intersection, fill_color=c("darkblue", "darkgreen"), show_percentage = F, set_name_size = 10, text_size = 10) #+ ggtitle("Noctilio leporinus L and UL") + theme(plot.title = element_text(size = 20))
p4_1 <- ggplot(noclep1_L_UL_uniq_status_noclep2, aes(x=factor(HLnocLep2)))+geom_bar(stat="count")+theme_minimal()+ggtitle("L+UL HLnocLep1 only")+xlab("state in HLnocLep2")+theme(axis.text=element_text(size=14), axis.title=element_text(size=14), plot.title = element_text(size = 14))
p4_2 <- ggplot(noclep2_L_UL_uniq_status_noclep1, aes(x=factor(HLnocLep1)))+geom_bar(stat="count")+theme_minimal()+ggtitle("L+UL HLnocLep2 only")+xlab("state in HLnocLep1")+theme(axis.text=element_text(size=14), axis.title=element_text(size=14), plot.title = element_text(size = 14))

ggarrange(p4, p2, nrow = 1)
ggarrange(p4_1, p2_1, p4_2, p2_2, nrow = 2)

### Myotis vivesi
all_species = colnames(gene_loss_full)
gene_loss_viv_table = gene_loss_full[gene_loss_full$HLmyoViv1 %notin% c("I", "FI") & gene_loss_full$HLmyoViv5 == "L",]# & gene_loss_full$HLmyoViv3 == "L" & gene_loss_full$HLmyoViv5 == "L",]
gene_loss_lep_table = gene_loss_full[gene_loss_full$HLnocLep1 %notin% c("I", "FI") & gene_loss_full$HLnocLep2 == "L",]#, "UL"),]

gene_L_UL_viv_table = gene_loss_full[gene_loss_full$HLmyoViv1 %notin% c("I", "FI") & gene_loss_full$HLmyoViv5 %in% c("L", "UL"),]# & gene_loss_full$HLmyoViv3 %in% c("L", "UL") & gene_loss_full$HLmyoViv5 %in% c("L", "UL"),]
gene_L_UL_lep_table = gene_loss_full[gene_loss_full$HLnocLep1 %notin% c("I", "FI") & gene_loss_full$HLnocLep2 %in% c("L", "UL"),]#, "UL"),]
gene_common_L_UL_viv_lep = gene_loss_full[gene_loss_full$HLnocLep1 %notin% c("I", "FI") & gene_loss_full$HLnocLep2 %in% c("L", "UL") & gene_loss_full$HLmyoViv1 %notin% c("I", "FI") & gene_loss_full$HLmyoViv5 %in% c("L", "UL"),]

losses_by_background_viv = data.frame(row.names = c("node_number", rownames(gene_loss_viv_table)))
losses_by_background_lep = data.frame(row.names = c("node_number", rownames(gene_loss_lep_table)))

L_UL_by_background_viv = data.frame(row.names = c("node_number", rownames(gene_L_UL_viv_table)))
L_UL_by_background_lep = data.frame(row.names = c("node_number", rownames(gene_L_UL_lep_table)))

common_L_UL_viv_lep = data.frame(row.names = c("node_number", rownames(gene_common_L_UL_viv_lep)))

foreground = c("HLmyoViv5", "HLnocLep2")
midground = c("HLmyoDau2", "HLmyoRic2")

for (label1 in all_species){
  for (label2 in all_species) {
    label_in_tree = paste0(label1, "_", label2)
    if (label_in_tree %in% c(tree$node.label)) {
      node_number = MRCA(tree, label_in_tree)
      tips = tree$tip.label[unlist(Descendants(tree, node_number, type = c("tips")))]
      if ("HLmyoViv5" %in% tips) {
        if (length(tips) < 4) {
          losses_by_background_viv[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_loss_viv_table) %in% 
                                                            c(rownames(gene_loss_viv_table %>% 
                                                                       filter(if_all(tips[tips != "HLmyoViv5"], 
                                                                                     ~ .x %in% c("FI", "I")))))))
          L_UL_by_background_viv[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_L_UL_viv_table) %in% 
                                                            c(rownames(gene_L_UL_viv_table %>% 
                                                                       filter(if_all(tips[tips != "HLmyoViv5"], 
                                                                                     ~ .x %in% c("FI", "I")))))))
        }
        else {
          losses_by_background_viv[,label_in_tree] = c(node_number, 
                                                        c(rownames(gene_loss_viv_table) %in% 
                                                            rownames(gene_loss_viv_table[(rowMeans(gene_loss_viv_table[,tips[tips %notin% c(foreground, midground)]] == "I") +
                                                                                            rowMeans(gene_loss_viv_table[,tips[tips %notin% c(foreground, midground)]] == "FI")) >= 0.9,])))
          L_UL_by_background_viv[,label_in_tree] = c(node_number, 
                                                        c(rownames(gene_L_UL_viv_table) %in% 
                                                            rownames(gene_L_UL_viv_table[(rowMeans(gene_L_UL_viv_table[,tips[tips %notin% c(foreground, midground)]] == "I") +
                                                                                            rowMeans(gene_L_UL_viv_table[,tips[tips %notin% c(foreground, midground)]] == "FI")) >= 0.9,])))
          common_L_UL_viv_lep[,label_in_tree] = c(node_number, 
                                                        c(rownames(gene_common_L_UL_viv_lep) %in% 
                                                            rownames(gene_common_L_UL_viv_lep[(rowMeans(gene_common_L_UL_viv_lep[,tips[tips %notin% c(foreground, midground)]] == "I") +
                                                                                            rowMeans(gene_common_L_UL_viv_lep[,tips[tips %notin% c(foreground, midground)]] == "FI")) >= 0.5,])))
        }
      }
      if ("HLnocLep2" %in% tips) {
        if (length(tips) < 3) {
          losses_by_background_lep[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_loss_lep_table) %in% 
                                                            c(rownames(gene_loss_lep_table %>% 
                                                                       filter(if_all(tips[tips != "HLnocLep2"], 
                                                                                     ~ .x %in% c("FI", "I")))))))
          L_UL_by_background_lep[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_L_UL_lep_table) %in% 
                                                            c(rownames(gene_L_UL_lep_table %>% 
                                                                       filter(if_all(tips[tips != "HLnocLep2"], 
                                                                                     ~ .x %in% c("FI", "I")))))))
        }
        else {
          losses_by_background_lep[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_loss_lep_table) %in% 
                                                            rownames(gene_loss_lep_table[(rowMeans(gene_loss_lep_table[,tips[tips %notin% c(foreground, midground)]] == "I") +
                                                                                            rowMeans(gene_loss_lep_table[,tips[tips %notin% c(foreground, midground)]] == "FI")) >= 0.9,])))
          L_UL_by_background_lep[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_L_UL_lep_table) %in% 
                                                            rownames(gene_L_UL_lep_table[(rowMeans(gene_L_UL_lep_table[,tips[tips %notin% c(foreground, midground)]] == "I") +
                                                                                            rowMeans(gene_L_UL_lep_table[,tips[tips %notin% c(foreground, midground)]] == "FI")) >= 0.9,])))
        }
      }
    }
  }
}

losses_by_background_viv = losses_by_background_viv[-1,order(as.numeric(losses_by_background_viv["node_number",]), decreasing = TRUE)]
losses_by_background_lep = losses_by_background_lep[-c(1),order(as.numeric(losses_by_background_lep["node_number",]), decreasing = TRUE)]

L_UL_by_background_viv = L_UL_by_background_viv[-c(1),order(as.numeric(L_UL_by_background_viv["node_number",]), decreasing = TRUE)]
L_UL_by_background_lep = L_UL_by_background_lep[-c(1),order(as.numeric(L_UL_by_background_lep["node_number",]), decreasing = TRUE)]
common_L_UL_viv_lep = common_L_UL_viv_lep[-c(1),order(as.numeric(common_L_UL_viv_lep["node_number",]), decreasing = TRUE)]
strsplit(rownames(common_L_UL_viv_lep[common_L_UL_viv_lep$HLmyzAur1_HLnycThe1A == 1,]), "#") %>% lapply(function(x) x[2]) %>% unlist
common_L_UL_viv_lep = common_L_UL_viv_lep[common_L_UL_viv_lep$HLmyzAur1_HLnycThe1A == 1,]
common_L_UL_all_states = gene_common_L_UL_viv_lep[rownames(common_L_UL_viv_lep), colnames(gene_common_L_UL_viv_lep) %notin% c("HLmyoViv1", "HLmyoViv2", "HLmyoViv3", "HLnocLep1")]
species = synonyms_df[order(match(synonyms_df$old, colnames(common_L_UL_all_states))),]
colnames(common_L_UL_all_states) = species[species$old %in% colnames(common_L_UL_all_states),]$new
tree$label[MRCA(tree, "HLnocLep2")] = "focus\nspecies"
tree$label[MRCA(tree, "HLmyoViv5")] = "focus\nspecies"
renamed_tree <- rename_taxa(tree, data = synonyms_df, key = "old", value = "new")
tree_common = ggplot(renamed_tree, aes(color=tree$label)) + scale_color_manual(values=c("focus\nspecies"="#d1182a", "background"="#555555")) + 
  geom_tree() + theme_tree() + geom_tiplab(align=TRUE, fontface = "italic", size=7)
rownames(common_L_UL_all_states) = strsplit(rownames(common_L_UL_all_states),"#") %>% lapply(function(x) x[2]) %>% unlist
common_L_UL_all_states = common_L_UL_all_states[rownames(common_L_UL_all_states) != "RAET1E",]
p1 = gheatmap(tree_common, as.data.frame(t(common_L_UL_all_states)), offset=20, width = 3, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_manual(values=setNames(gene.colors, gene.states), name="gene state") + 
  theme(axis.text=element_text(size=25), axis.text.x=element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=15)) + ggtree::vexpand(0, -1) + vexpand(-1)#+
p2 = p1 + new_scale_fill()
gheatmap(p2, table_bats, offset=75, width = .6, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_continuous(name="Diet code", breaks = diet.states, high = "#132B43", low = "#56B1F7") +
  theme(axis.text=element_text(size=20), axis.text.x=element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=15)) + ggtree::vexpand(0, -1)

losses_by_background_viv = losses_by_background_viv[losses_by_background_viv$HLmyoViv5_HLmyoLuc2 != 0,]
losses_by_background_viv = losses_by_background_viv[order(losses_by_background_viv[,2]),]
losses_by_background_viv = losses_by_background_viv[losses_by_background_viv$HLmyzAur1_HLnycThe1A == 1,]
losses_viv_all_states = gene_loss_viv_table[rownames(losses_by_background_viv), colnames(gene_loss_viv_table) %notin% c("HLmyoViv1", "HLmyoViv2", "HLmyoViv3", "HLnocLep1")]
colnames(losses_viv_all_states) = species[species$old %in% colnames(losses_viv_all_states),]$new
tree$label[MRCA(tree, "HLnocLep2")] = "background"
tree$label[MRCA(tree, "HLmyoViv5")] = "focus\nspecies"
renamed_tree <- rename_taxa(tree, data = synonyms_df, key = "old", value = "new")
tree_viv = ggplot(renamed_tree, aes(color=tree$label)) + scale_color_manual(values=c("focus\nspecies"="#d1182a", "background"="#555555")) + 
  geom_tree() + theme_tree() + geom_tiplab(align=TRUE, fontface = "italic", size=7)
rownames(losses_viv_all_states) = strsplit(rownames(losses_viv_all_states),"#") %>% lapply(function(x) x[2]) %>% unlist
gheatmap(tree_viv, as.data.frame(t(losses_viv_all_states)), offset=6, width = 0.2, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_manual(values=setNames(gene.colors, gene.states), name="gene state") + 
  theme(axis.text.x=element_blank()) + ggtree::vexpand(0, -1) + vexpand(-1)

L_UL_by_background_viv = L_UL_by_background_viv[L_UL_by_background_viv$HLmyoViv5_HLmyoLuc2 != 0,]
L_UL_by_background_viv = L_UL_by_background_viv[order(L_UL_by_background_viv[,2]),]
L_UL_by_background_viv = L_UL_by_background_viv[L_UL_by_background_viv$HLmyzAur1_HLnycThe1A == 1,]
L_UL_viv_all_states = gene_L_UL_viv_table[rownames(L_UL_by_background_viv), colnames(gene_L_UL_viv_table) %notin% c("HLmyoViv1", "HLmyoViv2", "HLmyoViv3", "HLnocLep1")]
colnames(L_UL_viv_all_states) = species[species$old %in% colnames(L_UL_viv_all_states),]$new
L_UL_viv_all_states_for_heat <- melt(as.matrix(L_UL_viv_all_states))
colnames(L_UL_viv_all_states_for_heat) <- c("y", "x", "value")
rownames(L_UL_viv_all_states) = strsplit(rownames(L_UL_viv_all_states),"#") %>% lapply(function(x) x[2]) %>% unlist
L_UL_viv_all_states = L_UL_viv_all_states[rownames(L_UL_viv_all_states) %notin% c("LSS", "CD81", "ODAD2", "DDX4", "FLT4", "CPZ"),]
p1 = gheatmap(tree_viv, as.data.frame(t(L_UL_viv_all_states)), offset=15, width = 1.5, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_manual(values=setNames(gene.colors, gene.states), name="gene state") + 
  theme(axis.text.x=element_blank(), axis.text=element_text(size=20)) + ggtree::vexpand(0, -1) + vexpand(-1)#+
p2 = p1 + new_scale_fill()
gheatmap(p2, table_bats, offset=43, width = .3, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_continuous(name="Diet code", breaks = diet.states, high = "#132B43", low = "#56B1F7") +
  theme(axis.text=element_text(size=20), axis.text.x=element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=15)) + ggtree::vexpand(0, -1)

losses_by_background_lep = losses_by_background_lep[losses_by_background_lep[,1] != 0,]
losses_by_background_lep = losses_by_background_lep[losses_by_background_lep$HLmyzAur1_HLnycThe1A == 1,]
losses_lep_all_states = gene_loss_lep_table[rownames(losses_by_background_lep), colnames(gene_loss_lep_table) %notin% c("HLmyoViv1", "HLmyoViv2", "HLmyoViv3", "HLnocLep1")]
colnames(losses_lep_all_states) = species[species$old %in% colnames(losses_lep_all_states),]$new
losses_lep_all_states_for_heat <- melt(as.matrix(losses_lep_all_states))
colnames(losses_lep_all_states_for_heat) <- c("y", "x", "value")
tree$label[MRCA(tree, "HLnocLep2")] = "focus\nspecies"
tree$label[MRCA(tree, "HLmyoViv5")] = "background"
renamed_tree <- rename_taxa(tree, data = synonyms_df, key = "old", value = "new")
tree_lep = ggplot(renamed_tree, aes(color=tree$label)) + scale_color_manual(values=c("focus\nspecies"="#d1182a", "background"="#555555")) + 
  geom_tree() + theme_tree() + geom_tiplab(align=TRUE, fontface = "italic", size=7)
rownames(losses_lep_all_states) = strsplit(rownames(losses_lep_all_states),"#") %>% lapply(function(x) x[2]) %>% unlist
gheatmap(tree_lep, as.data.frame(t(losses_lep_all_states)), offset=6, width = 1.4, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_manual(values=setNames(gene.colors, gene.states), name="gene state") + 
  theme(axis.text.x=element_blank()) + ggtree::vexpand(0, -1) + vexpand(-1)
                                                                                           
L_UL_by_background_lep = L_UL_by_background_lep[L_UL_by_background_lep[,1] != 0,]
L_UL_by_background_lep = L_UL_by_background_lep[L_UL_by_background_lep$HLmyzAur1_HLnycThe1A == 1,]
L_UL_lep_all_states = gene_L_UL_lep_table[rownames(L_UL_by_background_lep), colnames(gene_L_UL_lep_table) %notin% c("HLmyoViv1", "HLmyoViv2", "HLmyoViv3", "HLnocLep1")]
colnames(L_UL_lep_all_states) = species[species$old %in% colnames(L_UL_lep_all_states),]$new
rownames(L_UL_lep_all_states) = strsplit(rownames(L_UL_lep_all_states),"#") %>% lapply(function(x) x[2]) %>% unlist
L_UL_lep_all_states = L_UL_lep_all_states[rownames(L_UL_lep_all_states) %notin% c("ZNF446", "CCDC158", "EIF4E1B", "CA13", "ATN1", "DNAJB13", "SCML1", "SCML2",
                                                                                  "TRBJ2-2", "TRAJ25", "TRAJ32", "ENSG00000288644", "PPP1R13L"),]
p1 = gheatmap(tree_lep, as.data.frame(t(L_UL_lep_all_states)), offset=35, width = 9, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_manual(values=setNames(gene.colors, gene.states), name="gene state") + 
  theme(axis.text.x=element_blank()) + ggtree::vexpand(0, -1) + vexpand(-1)#+
p2 = p1 + new_scale_fill()
gheatmap(p2, table_bats, offset=195, width = .6, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_continuous(name="Diet code", breaks = diet.states, high = "#132B43", low = "#56B1F7") +
  theme(axis.text.x=element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=15)) + ggtree::vexpand(0, -1)

## Duplicated genes TOGA2

gene_dupl_full = read.csv("toga2_loss/dupl_genes_16072025.txt", sep="\t", row.names = 1, header = T) %>% as.data.frame()
names(gene_dupl_full)[names(gene_dupl_full) == 'HLnocLep2_upd'] <- 'HLnocLep2'
gene_dupl_full = gene_dupl_full[,colnames(gene_dupl_full) %notin% c("HLmyoViv2", "HLmyoViv3")]
all_species = colnames(gene_dupl_full)
gene_dupl_viv_table = gene_dupl_full[gene_dupl_full$HLmyoViv1 == "one2many" & gene_dupl_full$HLmyoViv5 == "one2many",]# & gene_dupl_full$HLmyoViv3 == "one2many" & gene_dupl_full$HLmyoViv5 == "one2many",]
gene_dupl_lep_table = gene_dupl_full[gene_dupl_full$HLnocLep1 == "one2many" & gene_dupl_full$HLnocLep2 == "one2many",]

dupl_by_background_viv = data.frame(row.names = c("node_number", rownames(gene_dupl_viv_table)))
dupl_by_background_lep = data.frame(row.names = c("node_number", rownames(gene_dupl_lep_table)))

for (label1 in all_species){
  for (label2 in all_species) {
    label_in_tree = paste0(label1, "_", label2)
    if (label_in_tree %in% c(tree$node.label)) {
      node_number = MRCA(tree, label_in_tree)
      tips = tree$tip.label[unlist(Descendants(tree, node_number, type = c("tips")))]
      if ("HLmyoViv5" %in% tips) {
        if (length(tips) < 4) {
          dupl_by_background_viv[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_dupl_viv_table) %in% 
                                                            c(rownames(gene_dupl_viv_table %>% 
                                                                       filter(if_all(tips[tips != "HLmyoViv5"], 
                                                                                     ~ .x == "one2one"))))))
        }
        else {
          dupl_by_background_viv[,label_in_tree] = c(node_number, 
                                                        c(rownames(gene_dupl_viv_table) %in% 
                                                            rownames(gene_dupl_viv_table[rowMeans(gene_dupl_viv_table[,tips[tips %notin% c("HLmyoViv5")]] == "one2one") >= 0.9 ,])))
        }
      }
      if ("HLnocLep2" %in% tips) {
        if (length(tips) < 3) {
          dupl_by_background_lep[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_dupl_lep_table) %in% 
                                                            c(rownames(gene_dupl_lep_table %>% 
                                                                       filter(if_all(tips[tips != "HLnocLep2"], 
                                                                                     ~ .x == "one2one"))))))
        }
        else {
          dupl_by_background_lep[,label_in_tree] = c(node_number,
                                                        c(rownames(gene_dupl_lep_table) %in% 
                                                            rownames(gene_dupl_lep_table[rowMeans(gene_dupl_lep_table[,tips[tips != "HLnocLep2"]] == "one2one") >= 0.9,])))
        }
      }
    }
  }
}
dupl_by_background_viv = dupl_by_background_viv[-c(1),order(as.numeric(dupl_by_background_viv["node_number",]), decreasing = TRUE)]
dupl_by_background_lep = dupl_by_background_lep[-c(1),order(as.numeric(dupl_by_background_lep["node_number",]), decreasing = TRUE)]
dupl.colors = c("blue", "lightblue", "orange", "gray50", "yellow")
dupl.states = c("many2many", "many2one", "one2many", "one2one", "one2zero")

dupl_by_background_viv = dupl_by_background_viv[dupl_by_background_viv$HLmyoViv5_HLmyoLuc2 != 0,]
dupl_by_background_viv = dupl_by_background_viv[order(dupl_by_background_viv[,2]),]
dupl_by_background_viv = dupl_by_background_viv[dupl_by_background_viv$HLmyzAur1_HLnycThe1A == 1,]
dupl_viv_all_states = gene_dupl_viv_table[rownames(dupl_by_background_viv), colnames(gene_dupl_viv_table) %notin% c("HLmyoViv1", "HLmyoViv2", "HLmyoViv3", "HLnocLep1")]
dupl_viv_all_states = dupl_viv_all_states[,colnames(dupl_viv_all_states) %in% species$old]
colnames(dupl_viv_all_states) = species[species$old %in% colnames(dupl_viv_all_states),]$new
#dupl_viv_all_states = dupl_viv_all_states[c("FABP3", "ORMDL2", "UQCRQ", "PEA15", "NDUFB1", "RPL22"),]
rownames(dupl_viv_all_states) = strsplit(rownames(dupl_viv_all_states),"#") %>% lapply(function(x) x[2]) %>% unlist
dupl_viv_all_states = dupl_viv_all_states[rownames(dupl_viv_all_states) %notin% c("SCT", "MED9"),]
tree$label[MRCA(tree, "HLmyoViv5")] = "focus\nspecies"
tree$label[MRCA(tree, "HLnocLep2")] = "background"
renamed_tree <- rename_taxa(tree, data = synonyms_df, key = "old", value = "new")
p1 = gheatmap(tree_viv, as.data.frame(t(dupl_viv_all_states)), offset=10, width = .2, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_manual(values=setNames(dupl.colors, dupl.states), name="dupl state") + 
  theme(axis.text = element_text(size=15), axis.text.x=element_blank()) + ggtree::vexpand(0, -1) + vexpand(-1)
p2 = p1 + new_scale_fill()
gheatmap(p2, table_bats, offset=15, width = .2, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_continuous(name="Diet code", breaks = diet.states, high = "#132B43", low = "#56B1F7") +
  theme(axis.text.x=element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=15)) + ggtree::vexpand(0, -1)

dupl_by_background_lep = dupl_by_background_lep[dupl_by_background_lep[,1] != 0,]
dupl_by_background_lep = dupl_by_background_lep[order(dupl_by_background_lep[,2]),]
dupl_by_background_lep = dupl_by_background_lep[dupl_by_background_lep$HLmyzAur1_HLnycThe1A == 1,]
dupl_lep_all_states = gene_dupl_lep_table[rownames(dupl_by_background_lep), colnames(gene_dupl_lep_table) %notin% c("HLmyoViv1", "HLmyoViv2", "HLmyoViv3", "HLnocLep1")]
dupl_lep_all_states = dupl_lep_all_states[,colnames(dupl_lep_all_states) %in% species$old]
colnames(dupl_lep_all_states) = species[species$old %in% colnames(dupl_lep_all_states),]$new
rownames(dupl_lep_all_states) = strsplit(rownames(dupl_lep_all_states),"#") %>% lapply(function(x) x[2]) %>% unlist
dupl_lep_all_states = dupl_lep_all_states[rownames(dupl_lep_all_states) %notin% c("ALG13", "TSPAN31", "CTXN1"),]
tree$label[MRCA(tree, "HLnocLep2")] = "focus\nspecies"
tree$label[MRCA(tree, "HLmyoViv5")] = "background"
renamed_tree <- rename_taxa(tree, data = synonyms_df, key = "old", value = "new")
p1 = gheatmap(tree_lep, as.data.frame(t(dupl_lep_all_states)), offset=10, width = .2, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_manual(values=setNames(dupl.colors, dupl.states), name="dupl state") + 
  theme(axis.text = element_text(size=15), axis.text.x=element_blank()) + ggtree::vexpand(0, -1) + vexpand(-1)
p2 = p1 + new_scale_fill()
gheatmap(p2, table_bats, offset=10, width = .2, color = "black", colnames_angle=90, hjust = 1) +
  scale_fill_continuous(name="Diet code", breaks = diet.states, high = "#132B43", low = "#56B1F7") +
  theme(axis.text.x=element_blank(), legend.text = element_text(size=15), legend.title = element_text(size=15)) + ggtree::vexpand(0, -1)

                                                                                           
