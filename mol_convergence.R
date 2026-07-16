library(RERconverge)
library(pathview)
library(ggplot2)

## only using bats
table_bats <- read.csv("Yangochiroptera77.csv")

###################################################################
################### generate gene trees ############################
tree_rer = "tree_for_selection_screen_t2_BL.tree"
ali_files = "fg_srv_prank_ali"

library(parallel)
# Run the tree estimation using 4 cores
estimatePhangornTreeAll.function = function(x){
  estimatePhangornTreeAll(alnfiles = paste0(ali_files, "/", list.files(ali_files)[x]),
                          treefile = tree_rer, 
                          output.file = paste0("outputs/output.", x), 
                          type = "DNA", submodel = "GTR")
}
length(list.files(ali_files))
mclapply(1:17018, estimatePhangornTreeAll.function, mc.cores = 4)
# cat outputs/* > estimatePhangornTreeAll.txt

###################################################################
########################## run RERconverge ###################################
treefile = readTrees("estimatePhangornTreeAll.txt", useSpecies = table_bats$Assembly)

#Estimating relative evolutionary rates (RER) with getAllResiduals
batRERw = getAllResiduals(treefile, useSpecies = table_bats$Assembly)
saveRDS(batRERw, file="batRERw.rds")
# batRERw = readRDS("batRERw.rds")

# Import Pathway Annotations
# downloaded from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/
## using the canonical pathways
annots = read.gmt("c2.cp.v2025.1.Hs.symbols.gmt")
# reformat
annotslist = list(annots)
names(annotslist) = "MSigDBpathways"

fish_eaters_2 <- c("HLnocLep2", "HLmyoViv5")
fish_eaters_4 <- c("HLnocLep2", "HLmyoViv5", "HLmyoRic2", "HLmyoDau2")

for(i in c(2,4)){
  #Generating paths using foreground2Paths which does not require a binary tree as input
  phenvFish <- foreground2Paths(get(paste0("fish_eaters_", i)),
                                useSpecies = table_bats$Assembly,
                                treefile, clade="terminal")
  assign(paste0("phenvFish_", i), phenvFish)
  
  #Correlating gene evolution with binary trait evolution, default no bootstrap
  # default Kendall rank correlation coefficient, or Tau
  corFish = correlateWithBinaryPhenotype(batRERw, phenvFish, min.sp=10, min.pos=2,  weighted="auto")
  
  # Calculates Rho-signed negative log-base-ten p-value for use in enrichment functions (NA removed)
  stats = getStat(corFish)
  corFish$stats = sapply(row.names(corFish), FUN=function(x){ stats[x] })
  assign(paste0("corFish_", i), corFish)
  
  # enrichment 
  names(stats) = sub(".*#", "", names(stats))
  enrichment = fastwilcoxGMTall(stats, annotslist, alternative = "two.sided", outputGeneVals = T, num.g = 10)
  assign(paste0("enrich_", i), enrichment$MSigDBpathways)
  
  rm(phenvFish, corFish, stats, enrichment)
}

