#!/usr/bin/env Rscript

#####################################
## name: plot_trees.R
## authors: Edyth Parker <ep484@cam.ac.uk>,
## Dylan Morris <dhmorris@princeton.edu>
##
## plots phylogenetic trees for
## Morris et al paper on asynchrony
## between diversity and selection
####################################


###############################
## package loading
###############################
library(ape, quietly = TRUE)
library(ggtree, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(colorspace, quietly = TRUE)



###############################
## styling pre-definitions / setup
###############################
col_vector <- c(sequential_hcl(5, "viridis")[5], 
                sequential_hcl(5, "blues")[2],
                sequential_hcl(5, "Light Grays")[3])

###############################
## arg parsing
###############################

args <- commandArgs(trailingOnly=TRUE)

H1_nwk <- args[1]
H3_08_nwk <- args[2]
H3_12_nwk <- args[3]

H1_sub <- args[4]
H3_08_sub <- args[5]
H3_12_sub <- args[6]


H1_path <- args[7]
H3_08_path <- args[8]
H3_12_path <- args[9]

###############################
## H1N1 tree
###############################

cat("Plotting H1 tree...\n\n")
cat("Reading newick file...\n")
tree <- read.tree(H1_nwk)
tree <- groupClade(tree, .node=c(3564, 5028, 5547))
annodf=read.table(H1_sub, sep="\t", stringsAsFactors = F)
names(annodf)<-c("TAXA", "LINEAGE")
annodf$LINEAGE[!(annodf$LINEAGE=="K") & !(annodf$LINEAGE=="E")]<-"Other"
annodf$LINEAGE = factor(annodf$LINEAGE, levels = c("E", "K", "Other"))
cat("Plotting...\n")
f1<-ggtree(tree, size=0.3, ladderize = TRUE) %<+% annodf+
  geom_hilight(node=3564, fill="steelblue", alpha=.2) +
  geom_hilight(node=5028, fill="darkgreen", alpha=.2) +
  geom_hilight(node=5547, fill="firebrick", alpha=.2)+
  geom_treescale(fontsize=3,offset=-60) +
  geom_cladelabel(node=5547, label="3", fontsize=5) +
  geom_cladelabel(node=5028, label="2", fontsize=5) +
  geom_cladelabel(node=3564, label="1", fontsize=5) +
  geom_tippoint(aes(fill=LINEAGE),size=2,pch = I(21), color="black")+
  scale_fill_manual(values=c(col_vector))


edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge))
colnames(edge)=c("parent", "node", "edge_num")
edge$edge_num<-NA
edge$edge_num[42]<-"Y94H"
edge$edge_num[4022]<-"K73R"
edge$edge_num[4030]<-"V128A"
edge$edge_num[4039]<-"A128T"
edge$edge_num[4043]<-"K140E"
edge$edge_num[4145]<-"P270S"
edge$edge_num[2985]<-"S36N"
edge$edge_num[2988]<-"A189T"
edge$edge_num[2996]<-"T193K"
edge$edge_num[2999]<-"K140E"
edge$edge_num[3583]<-"N244S"
edge$edge_num[3005]<-"K82R"
edge$edge_num[3020]<-"I47K"
edge$edge_num[3028]<-"E68G"
edge$edge_num[50]<-"R188K"
edge$edge_num[58]<-"E273K"
edge$edge_num[65]<-"K140E"
edge$edge_num[70]<-"D35N"
edge$edge_num[71]<-"K145R"



f1%<+% edge + geom_text2(aes(x=branch, label=edge_num), size=1.6, color="firebrick", vjust=-0.5, fontface = "bold")+
  theme(legend.position = c(0, 1), legend.title = element_text(size=12),
       legend.justification = c(0, 1), legend.text=element_text(size=10), legend.background = element_blank(),
       axis.text=element_text(size=3))+
  labs(col="Amino acid at position 140")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  annotate("text", x=0.067, y=1650, label="R188M", size=1.6, color="firebrick", fontface = "bold")+
  labs(fill = "Amino acid at position 140")

cat("Saving to ", H1_path, "...\n\n")
ggsave(H1_path, width = 210, height = 297, units = "mm")


###############################
## H3N2 tree, 2008
###############################

cat("Plotting H3 2008 tree...\n\n")
cat("Reading newick file...\n")
tree<-read.tree(H3_08_nwk)
tree <- groupClade(tree, .node=c(8380, 6626))
annodf=read.table(H3_08_sub, sep="\t", header=T, stringsAsFactors = F)
names(annodf)<-c("TAXA", "LINEAGE")
annodf$LINEAGE[!(annodf$LINEAGE=="158K/189N") & !(annodf$LINEAGE=="158N/189K")]<-"Other"
annodf$LINEAGE = factor(annodf$LINEAGE, levels = c("158N/189K",
                                                   "158K/189N",
                                                   "Other"))

cat("Plotting...\n")
f1<-ggtree(tree, size=0.3, ladderize = TRUE) %<+% annodf+
  geom_hilight(node=8380, fill="steelblue", alpha=.2) +
  geom_hilight(node=6626, fill="darkgreen", alpha=.2) +
  geom_treescale(fontsize=3,offset=-80) +
  geom_cladelabel(node=8380, label="1", fontsize=5)+
  geom_cladelabel(node=6626, label="2", fontsize=5)+
  geom_tippoint(aes(fill=LINEAGE),size=1.8,pch = I(21), color="black")+
  scale_fill_manual(values=c(col_vector))


edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge))
colnames(edge)=c("parent", "node", "edge_num")
edge$edge_num<-NA
edge$edge_num[1769]<-"K189N"
edge$edge_num[1771]<-"E62K"
edge$edge_num[1773]<-"K158N, N144K"
edge$edge_num[2020]<-"R261Q"
edge$edge_num[2090]<-"I260M"
edge$edge_num[2196]<-"P162S"
edge$edge_num[2266]<-"E50K"
edge$edge_num[2655]<-"V213A"
edge$edge_num[3211]<-"N133D"
edge$edge_num[3299]<-"T212A"
edge$edge_num[3305]<-"R142G"
edge$edge_num[5282]<-"K158N, N189K"
edge$edge_num[5284]<-"T212A"
edge$edge_num[5701]<-"S45N"
edge$edge_num[6542]<-"T48A"
edge$edge_num[6546]<-"R42R"
edge$edge_num[6774]<-"Q57H"
edge$edge_num[7005]<-"A198S"
edge$edge_num[7017]<-"V223I"
edge$edge_num[7017]<-"N312S"
edge$edge_num[7326]<-"N278K"
edge$edge_num[7338]<-"Q33R"
edge$edge_num[7791]<-"N145S"
edge$edge_num[8208]<-"G5E"
edge$edge_num[8208]<-"E62V"
edge$edge_num[5311]<-"V223I"
edge$edge_num[5313]<-"N145S"
edge$edge_num[8816]<-"D53N"
edge$edge_num[8816]<-"E280A"
edge$edge_num[8818]<-"I230V"
edge$edge_num[8820]<-"Y94H"
edge$edge_num[10400]<-"S199A"
edge$edge_num[7026]<-"S45N"
edge$edge_num[7026]<-"T46I"
edge$edge_num[6161]<-"N312S"
edge$edge_num[5331]<-"N144D"


f1%<+% edge +
  geom_text2(aes(x=branch, label=edge_num), size=1.6, color="firebrick", vjust=-0.5, fontface = "bold")+
  theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        axis.text=element_text(size=3),legend.background=element_blank())+
  labs(col="Amino acid at positions 158/189")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
    labs(fill = "Amino acid at positions 158 and 189")

cat("Saving to ", H3_08_path, "...\n\n")
ggsave(H3_08_path, width = 210, height = 297, units = "mm")


###############################
## H3N2 tree, 2012
############################### 

cat("Plotting H3 2012 tree...\n\n")
cat("Reading newick file...\n")
tree<-read.tree(H3_12_nwk)
annodf=read.table(H3_12_sub, sep="\t", header=T, stringsAsFactors = F)
annodf$COMBO = factor(annodf$COMBO, levels = c("Y",
                                               "F",
                                               "S"))

tree <- groupClade(tree, .node=c(11957, 14555))


cat("Plotting...\n")
f1<-ggtree(tree, size=0.3, ladderize = TRUE) %<+% annodf+
  geom_hilight(node=11957, fill="steelblue", alpha=.2) +
  geom_hilight(node=14555, fill="darkgreen", alpha=.2)+
  geom_cladelabel(node=11957, label="1", fontsize=5)+
  geom_cladelabel(node=14555, label="2", fontsize=5)+
    geom_tippoint(aes(fill=COMBO),size=1.8,pch = I(21), color="black")+
    scale_fill_manual(values=c(col_vector)) +
  labs(fill = "Amino acid at position 159")


edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge))
colnames(edge)=c("parent", "node", "edge_num")
edge$edge_num<-NA
edge$edge_num[8948]<-"L3I, Q311H, N144S"
edge$edge_num[8959]<-"F159Y"
edge$edge_num[11150]<-"R142K"
edge$edge_num[10474]<-"R261L"
edge$edge_num[68]<-"T128A"
edge$edge_num[3684]<-"A138S"
edge$edge_num[3686]<-"N225D"
edge$edge_num[3691]<-"F159S"


plot3<-f1%<+% edge + geom_text2(aes(x=branch, label=edge_num), size=1.6, color="firebrick", vjust=-0.5, fontface = "bold")+
  theme(legend.position = c(0, 1), legend.title = element_text(size=7),
        legend.justification = c(0, 1), legend.text=element_text(size=7),
        axis.text=element_text(size=3))+
  labs(col="Amino acid at position 159")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
    annotate("text", x=0.033, y=3600, label="N225D", size=1.6, color="firebrick", fontface = "bold")+
  annotate("text", x=0.0335, y=6050, label="R142G", size=1.6, color="firebrick", fontface = "bold")+
  annotate("text", x=0.0407, y=3605, label="K160T", size=1.6, color="firebrick", fontface = "bold")

cat("Saving to ", H3_12_path, "...\n\n")
ggsave(H3_12_path, width = 210, height = 297, units = "mm")


cat("All phylogenetic tree figures generated and saved successfully!\n")
