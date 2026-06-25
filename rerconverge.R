library(RERconverge)
library(pathview)
library(ggplot2)

tree_rer = "/home/ana/Documents/tree_for_selection_screen_t2_BL.tree"
ali_files = "fg_srv_prank_ali"
output = "est_phangorn_260426.txt"
estimatePhangornTreeAll(alndir = ali_files, treefile = tree_rer, output.file = output, type = "DNA", submodel = "GTR")

table_bats <- read.csv("Documents/loss/Yangochiroptera77.csv")

treefile = readTrees("est_phangorn_260426.txt")
names(treefile$lengths) <- gsub(".manual", "", names(treefile$lengths))
names(treefile$trees) <- gsub(".manual", "", names(treefile$trees))
batRERw = getAllResiduals(treefile, useSpecies = table_bats$Assembly)

fish_eaters <- c("HLnocLep2", "HLmyoViv5")
par(mfrow=c(1,2))
avgtree=plotTreeHighlightBranches(treefile$masterTree, outgroup="hg38",
                                  hlspecies=fish_eaters, hlcols=c("red"),
                                  main="Average tree") #plot average tree

par(mfrow=c(1,1))
phenvFish <- foreground2Paths(fish_eaters,treefile,clade="terminal", useSpecies = table_bats$Assembly)

fish2d = foreground2Tree(fish_eaters, treefile, clade="all", weighted = TRUE)

corFish=correlateWithBinaryPhenotype(batRERw, phenvFish, min.sp=10, min.pos=2,
                                     weighted="auto")
head(corFish[order(corFish$P),], 10)

stats = getStat(corFish)
names(stats) = sub(".*#", "", names(stats))
annots = read.gmt("Documents/c2.cp.v2025.1.Hs.symbols.gmt")
annotslist = list(annots)
names(annotslist) = "MSigDBpathways"
enrichment = fastwilcoxGMTall(stats, annotslist,alternative = "two.sided", outputGeneVals = T, num.g = 10)
head(enrichment)

hist(corFish$P, breaks=15, xlab="Kendall P-value")

enrich_df <- enrichment$MSigDBpathways
top_enrich <- enrich_df[order(enrich_df$p.adj), ][1:20, ]

ggplot(top_enrich, aes(x=reorder(rownames(top_enrich), log10(p.adj)), y=-log10(p.adj))) +
  geom_bar(stat="identity", fill="dodgerblue") +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(x="Pathway", y="-log10(p.adj)", title="Top Pathway Enrichments")
