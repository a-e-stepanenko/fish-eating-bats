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
library(egg)
library(org.Hs.eg.db)
library(ggnewscale)
library(treeio)
library(ggvenn)
library(wrapr)

`%notin%` <- Negate(`%in%`)
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

absrel_corr_1_3000 = cbind(read.csv("corr_1_3000.txt", sep = "\t", header = FALSE), V4 = "1_3000")
absrel_corr_3001_6000 = cbind(read.csv("corr_3001_6000.txt", sep = "\t", header = FALSE), V4 = "3001_6000")
absrel_corr_6001_9000 = cbind(read.csv("corr_6001_9000.txt", sep = "\t", header = FALSE), V4 = "6001_9000")
absrel_corr_9001_12000 = cbind(read.csv("corr_9001_12000.txt", sep = "\t", header = FALSE), V4 = "9001_12000")
absrel_corr_12001_15000 = cbind(read.csv("corr_12001_15000.txt", sep = "\t", header = FALSE), V4 = "12001_15000")
absrel_corr_15001_18075 = cbind(read.csv("corr_15001_18075.txt", sep = "\t", header = FALSE), V4 = "15001_18075")

absrel_corr = rbind(absrel_corr_1_3000, absrel_corr_3001_6000, absrel_corr_6001_9000, absrel_corr_9001_12000, absrel_corr_12001_15000, absrel_corr_15001_18075)  %>% data.frame
colnames(absrel_corr) = c("transcript", "label", "pval", "batch")

absrel_corr_sig = absrel_corr[absrel_corr$pval <= 0.05,]
absrel_corr_sig[c("transctipt2", "gene")] = str_split_fixed(absrel_corr_sig$transcript, '[#]', 2)
absrel_corr_sig = filter(absrel_corr_sig, !grepl("_",label))
myoviv_strict = absrel_corr_sig[absrel_corr_sig$label == "HLmyoViv5",]
index_df = table(absrel_corr_sig[(absrel_corr_sig$transcript %in% myoviv_strict$transcript) & (absrel_corr_sig$label %notin% c("HLnocLep2", "HLmyoDau2", "HLmyoRic2")),]$transcript) %>% t %>% data.frame
myoviv_strict$index = index_df[(index_df$Var2 %in% myoviv_strict$transcript) & order(index_df$Var2),]$Freq
myoviv_strict = myoviv_strict[myoviv_strict$index == 1,]
myoviv_strict

noclep_strict = absrel_corr_sig[absrel_corr_sig$label == "HLnocLep2",]
index_df = table(absrel_corr_sig[(absrel_corr_sig$transcript %in% noclep_strict$transcript) & (absrel_corr_sig$label %notin% c("HLmyoViv5", "HLmyoDau2", "HLmyoRic2")),]$transcript) %>% t %>% data.frame
noclep_strict$index = index_df[(index_df$Var2 %in% noclep_strict$transcript) & order(index_df$Var2),]$Freq
noclep_strict = noclep_strict[noclep_strict$index == 1,]
noclep_strict

absrel_corr_lax = absrel_corr[absrel_corr$pval <= 0.2,]
absrel_corr_lax[c("transctipt2", "gene")] = str_split_fixed(absrel_corr_lax$transcript, '[#]', 2)
absrel_corr_lax = filter(absrel_corr_lax, !grepl("_",label))
myoviv_lax = absrel_corr_lax[absrel_corr_lax$label == "HLmyoViv5",]
index_df = table(absrel_corr_lax[(absrel_corr_lax$transcript %in% myoviv_lax$transcript) & (absrel_corr_lax$label %notin% c("HLnocLep2", "HLmyoDau2", "HLmyoRic2")),]$transcript) %>% t %>% data.frame
myoviv_lax$index = index_df[(index_df$Var2 %in% myoviv_lax$transcript) & order(index_df$Var2),]$Freq
myoviv_lax

noclep_lax = absrel_corr_lax[absrel_corr_lax$label == "HLnocLep2",]
index_df = table(absrel_corr_lax[(absrel_corr_lax$transcript %in% noclep_lax$transcript) & (absrel_corr_lax$label %notin% c("HLmyoViv5", "HLmyoDau2", "HLmyoRic2")),]$transcript) %>% t %>% data.frame
noclep_lax$index = index_df[(index_df$Var2 %in% noclep_lax$transcript) & order(index_df$Var2),]$Freq
noclep_lax
