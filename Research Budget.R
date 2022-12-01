## Torpor phylogeny script
## Database originally from Ruf and Geiser 2015
## Additional data compiled by Santiago Tabares Erices and Grace Li Guo
## Supervised by Anusha Shankar
## Code by Anusha Shankar, github/nushiamme
## Started Oct 3, 2022

### Contents
## Setup, read files in, format data
## Figures

#library(MCMCglmm)
#library(nlme)
library(ape)
library(geiger) # for treedata() function
library(caper)
library(phytools)
#library(RColorBrewer)
library(ggplot2)
library(here)
library(readxl)
library(here)
library(tidyverse)
library(phyloch)  #remotes::install_github("fmichonneau/phyloch")
## For putting trait data on phylogeny
library(ggimage)
library(ggtree)
library(TDbook)
## To have two color scales
library(ggnewscale)

#### Setup ####
torpor_dat <- read.csv(here("..", "TorporExpandedPhylogeny_AGS.csv"))
torpor_dat$Use[torpor_dat$Use==""] <- torpor_dat$Species_Geiser[torpor_dat$Use==""]
#aviandiet <- read.csv(here())
#mammaldiet <- read.table(here())

## General plotting functions ####
## Generic theme
my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) 

# Axis labels
Tbmin.lab <- expression(atop(paste("Minimum torpid body temperature (", degree,"C)")))
bodymass.lab <- paste("log(Body Mass (g))")
# SurfTemp_short.lab <- expression(atop(paste("Ts (", degree,"C)")))
# SurfTemp.lab <- expression(atop(paste("Maximum Surface Temperature (", degree,"C)")))

##Fixed color scale for categories ####
my_colors <- c("#23988aff", "#440558ff", "#9ed93aff", "darkgoldenrod2") #NTDA
my_colors_five <- c("#23988aff",  "#F38BA8", "#440558ff",  "#9ed93aff", "darkgoldenrod2") #NSTDA
my_colors7 <- c("#ffe74c", "#508aa8", "#242f40", "#c60f7b", "#bbc7a4")
#names(my_colors) <- c("Normothermic", "ShallowTorpor", "Transition", "DeepTorpor", "Arousal")
colScale <- scale_colour_manual(name = "Category", values = my_colors)

## Function for plotting trees ####
f.treplot <-function(tree,...){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  plotTree(tree,fsize=fsize,lwd=1,...)
}

## Read phylogenetic tree for birds ####
avTree <- read.tree(here("..", "Stage2_Hackett_MCC_no_neg_fromEliotMiller.tre")) ## From Eliot, based on Jetz

f.treplot(avTree) 
avTree1 <- ladderize(avTree)

## IGNORE: Plot bird tree to PDF ####
pdf(file=here("..", "BirdPhylogeny",
              "AvianPhylo.pdf"), width=8.5, height=188, onefile=TRUE)

#quartz(width=8.5, height=144)
plot(avTree1, cex=0.3, label.offset=0.01, y.lim=c(130,3950))

node.support(avTree1$node.label, mode="dots", col = "red", cex=0.3)
node.support(avTree1$node.label, mode="dots", cutoff=75, col = "black", cex=0.35)
node.support(avTree1$node.label, mode="numbers", digits=2,pos="pretty", col = "black", font=2, cex=0.3)
add.scale.bar(cex=0.6)

dev.off()


## IGNORE: Source circular phylogeny functions ####
source(here("..", "ELE", "ELEcode",'DrawCircularPhylogeny.R'))

## IGNORE: Eliot suggested for Quick plotting and checking out data ####
phytools::phylosig(tree = avTree,)
contmap()

## Match tip labels from tree with our database for birds ####
torpor_av <- torpor_dat[torpor_dat$Class=="AVES",]
torpor_av$species_av <- str_replace(torpor_av$Species_Geiser, " ", "_")
species_av <- torpor_av$Species_Geiser
species_av <- str_replace(species_av, " ", "_")
torpor_av$newick_label <- species_av

## Making a data frame with all the bird species, columns without data will be empty
torpor_all_birds <- data.frame(matrix(nrow=length(avTree$tip.label),data = avTree$tip.label))
colnames(torpor_all_birds) <-  "newick_label"
torpor_all_birds$newick_label <- row.names(torpor_all_birds)


torpor_abirds <- torpor_all_birds %>%
  left_join(torpor_av, by = 'newick_label')

## Setting Type as "UK" for birds not in our database
torpor_abirds$Type <- as.character(torpor_abirds$Type)
torpor_abirds$Type[torpor_abirds$newick_label %in% setdiff(avTree$tip.label, species_av)] <- "UK" 
torpor_abirds$Type <- as.factor(torpor_abirds$Type)

torpor_abirds$newick_label[is.na(torpor_abirds$Species_Geiser)] <- as.numeric(sub("_", "", torpor_abirds$newick_label[is.na(torpor_abirds$Species_Geiser)]))

#check for mismatches in tree family tips labels and family names in tables
setdiff(species_av, avTree$tip.label) # appears in x but not y
intersect(species_av, avTree$tip.label) # appears in both x and y

## Find number of times a string occurs in the phylogeny
sum(grepl("Selasphorus", avTree$tip.label, fixed = TRUE), na.rm=T)

## Find what those names are in the phylogeny
avTree$tip.label[grepl("scandiaca", avTree$tip.label, fixed = TRUE)]

## Pruning full avian tree to match torpor tree, to test for phylogeny errors ####
torpor_av$species_av <- as.factor(torpor_av$species_av)
tips_av<-data.frame(levels(torpor_av$species_av))
colnames(tips_av) <- "tips"
rownames(tips_av)<-tips_av$tips

#match tree to data, prune tree, species names should be in rownnames of "data" 
tre_av<-treedata(avTree, tips)$phy

plot(tre_av, cex=0.5, edge.width = 0.5)

# cols<-setNames(c("black","blue","red", "green"),c("DT","HIB","ST", "NO"))
# plot(tre1,colors=cols,lwd=4)

f.treplot(tre_av) 


torpor_dat$Type <- as.factor(torpor_dat$Type)
torpor_av$Type <- as.factor(torpor_av$Type)
torpor_mam_pruned$Type <- droplevels(as.factor(torpor_mam_pruned$Type))



## Tree with tips colored by factor, for birds ####
## THIS WORKS
# From https://yulab-smu.top/treedata-book/chapter7.html
torpor_av <- torpor_av %>% relocate(newick_label, .before = Source)
ggtree(tre_av, layout = "circular") %<+% torpor_av + xlim(-.1, 150) +
  geom_tiplab(offset = .4, hjust=-0.2) +
  geom_tippoint(aes(col=Type, size = BM)) + 
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20)) +
  theme(legend.position = "right") + 
  scale_size_continuous(range = c(3, 8), name = "Body Mass (g)") +
  scale_colour_manual(values = my_colors) +
  guides(colour = guide_legend(override.aes = list(size=3)))

## Separate colors for Type and Tbmin
ggtree(tre_av, layout = "circular") %<+% torpor_av + xlim(-.1, 150) +
  geom_tiplab(aes(label=Species_Geiser, col=Type), align = T, size=3.5, hjust=-0.2) +
  scale_color_manual(values = c("black", "red", "green", "orange")) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  new_scale_color() +
  geom_tippoint(aes(col=Tbmin, size = BM)) + 
  scale_colour_viridis_c() +
  scale_size_continuous(range = c(3, 8), name = "Body Mass (g)") +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20)) +
  theme(legend.position = "right")
  
## Trying to plot all birds in the phylogeny
## Pruning full avian tree to match torpor tree, to test for phylogeny errors ####
torpor_abirds$newick_label <- as.factor(torpor_abirds$newick_label)
tips_av<-data.frame(levels(torpor_abirds$newick_label))
colnames(tips_av) <- "tips"
rownames(tips_av)<-tips_av$tips

#match tree to data, prune tree, species names should be in rownnames of "data" 
tre_av<-treedata(avTree, tips_av)$phy

plot(tre_av, cex=0.5, edge.width = 0.5)

torpor_abirds <- torpor_abirds %>% relocate(newick_label, .before = Source)
ggtree(tre_av, layout = "circular") %<+% torpor_abirds + xlim(-.1, 150) +
  geom_tiplab(aes(label=newick_label, col=Type), align = T, size=1, hjust=-0.2) +
  scale_color_manual(values = c("black", "red", "green", "orange", "grey")) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  new_scale_color() +
  geom_tippoint(aes(col=Tbmin, size = BM)) + 
  scale_colour_viridis_c() +
  scale_size_continuous(range = c(3, 8), name = "Body Mass (g)") +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20)) +
  theme(legend.position = "right")


## Mammal tree ####
global<-read.tree(here("..", "Upham2019_Data_S3_globalRAxML_files",
                       "RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick"))
global2<-drop.tip(global,"_Anolis_carolinensis") # just don't use drop.tip2 on an ML tree !!
global3<-ladderize(global2)

## IGNORE: Plot mammal tree to PDF ####
pdf(file=here("..", "MammalPhylogeny",
              "MamPhy_globalRAxML_FIN4_108bootstraps_4098mam_linear_agesBP.pdf"), width=8.5, height=188, onefile=TRUE)

## Or double-wide
# pdf(file=here("..", "MammalPhylogeny","MamPhy_globalRAxML_FIN4_108bootstraps_4098mam_linear_agesBP_2xWide.pdf"), 
#     width=22, height=238, onefile=TRUE)

#quartz(width=8.5, height=144)
plot(global3, cex=0.3, label.offset=0.01, y.lim=c(130,3950))

node.support(global3$node.label, mode="dots", col = "red", cex=0.3)
node.support(global3$node.label, mode="dots", cutoff=75, col = "black", cex=0.35)
node.support(global3$node.label, mode="numbers", digits=2,pos="pretty", col = "black", font=2, cex=0.3)
add.scale.bar(cex=0.6)

dev.off()

## Match tip labels from tree with our database, for mammals ####
torpor_mam <- torpor_dat[torpor_dat$Class=="MAMMALIA",]
species_mam <- torpor_mam$Species_Geiser
## Replace space between genus and species with underscore
species_mam <- str_replace(species_mam, " ", "_")
species_fam_order_mam <- paste(species_mam, toupper(torpor_mam$Family), torpor_mam$Order, sep = "_")
torpor_mam$newick_label <- species_fam_order_mam

#check for mismatches in tree family tips labels and family names in tables ####
setdiff(species_fam_order_mam, global$tip.label) # appears in x but not y
intersect(species_fam_order_mam, global$tip.label) # appears in both x and y

## For now, 11/29, removing mammal rows that don't match any names on phylogeny
unmatched_mam <- setdiff(species_fam_order_mam, global$tip.label) # appears in x but not y
torpor_mam_pruned <- torpor_mam %>% 
  filter(!newick_label %in% unmatched_mam)
setdiff(torpor_dat_pruned$newick_label, global$tip.label)

## Find number of times a string occurs in the phylogeny
sum(grepl("Spermophilus", global$tip.label, fixed = TRUE), na.rm=T)

## FOR SANTI: Can use the following line of code to test whether a particular string (i.e. a sequence of chracters) 
## appears in the phylogeny's tips
## Then if you figure out how to resolve it, go to the WorkingCopy Excel and change it in the Species_Geiser column
## So that it matches the spelling in the phylogeny. Once you're happy with the changes, 
## You can re-export the TorporTable2013 tab of the excel sheet
## To a csv named "TorporExpandedPhylogeny_AGS.csv" and can re-run the code to check.

## Find what those names are in the phylogeny
global$tip.label[grepl("Vespadelus", global$tip.label, fixed = TRUE)]

## Missing species
# Pseudomys_albocinereus missing from phylogeny?


## Tree with tips colored by factor, for mammals ####
## Pruning full mammal tree to match torpor tree, to test for phylogeny errors ####
torpor_mam_pruned$newick_label <- as.factor(torpor_mam_pruned$newick_label)
tips<-data.frame(levels(torpor_mam_pruned$newick_label))
colnames(tips) <- "tips"
rownames(tips)<-tips$tips

#match tree to data, prune tree, species names should be in rownnames of "data" 
tre_mam<-treedata(global, tips)$phy

plot(tre_mam, cex=0.5, edge.width = 0.5)

#cols_mam <-setNames(c("black","red"),c("DT","HIB"))
#plot(tre_mam,colors=cols_mam,lwd=4)
## THIS WORKS
# From https://yulab-smu.top/treedata-book/chapter7.html
torpor_mam_pruned$Type <- as.factor(torpor_mam_pruned$Type)
torpor_mam_pruned <- torpor_mam_pruned %>% relocate(newick_label, .before = Source)
torpor_mam_pruned$Type <- droplevels(torpor_mam_pruned$Type)
ggtree(tre_mam, layout = "circular", aes(col=Type)) %<+% torpor_mam_pruned + xlim(-.1, 2.5) +
  geom_tiplab(aes(label=Species_Geiser) ,align = T, size=2.75, hjust=-0.1) +
  geom_tippoint(aes(size = BM)) + 
  scale_colour_manual(values = my_colors) +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20)) +
  theme(legend.position = "right") + 
  scale_size_continuous(range = c(3, 8), name = "Body Mass (g)") +
  guides(color = guide_legend(override.aes = list(size=3)))

## Separate colors for Type and Tbmin
ggtree(tre_mam, layout = "circular") %<+% torpor_mam_pruned + xlim(-.1, 2.5) +
  geom_tiplab(aes(label=Species_Geiser, col=Type),align = T, size=2.75, hjust=-0.1) +
  scale_color_manual(values = c("black", "red", "gray")) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  new_scale_color() +
  geom_tippoint(aes(col=Tbmin, size = BM)) + 
  scale_colour_viridis_c() +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=20)) +
  theme(legend.position = "right") + 
  scale_size_continuous(range = c(3, 8), name = "Body Mass (g)")
  


df_ght <- torpor_mam_pruned[,c("Type", "BM", "Tbmin")]
row.names(df_ght) <- torpor_mam_pruned$newick_label
circ <- ggtree(tre_mam, layout = "circular") 
p1 <- gheatmap(circ, df_ght, offset=.8, width=.2) + #,
  #       colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")
p2 <- p1 + new_scale_fill()
gheatmap(p2, df_ght, offset=15, width=.3,
         colnames_angle=90, colnames_offset_y = .25) +
  scale_fill_viridis_c(option="A", name="continuous\nvalue")





## From eggshell wettability code ####
torpor_av$Type_categ <- torpor_av$Type # create variable you want to use
myvars <- c("Type_categ") #create list of variables
newdata <- torpor_av[myvars] #create a new dataframe with desired variables
Type_categ <- as.matrix(newdata)[,1] #your new matrix just has that variable.
names(Type_categ) <- species_av

setdiff(tre1$tip.label, species_av)

#This function performs fast estimation of the ML ancestral states for a continuous state
fit <- phytools::fastAnc(tre1,Type_categ,vars=TRUE,CI=TRUE) #Ancestral character estimate
td <- data.frame(node = nodeid(tree, names(response_continuous)), trait = response_continuous) #Ancestral character estimate
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree1 <- full_join(tree, d, by = 'node')

legend_title <- "response name" #put the name your response name here

tree_plot  =ggtree(tree1, aes(color=trait),
                   ladderize = FALSE, continuous = TRUE, size=0.5, layout='circular', right=T) +
  scale_color_gradientn(legend_title, colours=c("#820817", '#fc8d59', '#fee090', '#91bfdb', '#1a5199'))+ ##b2182b #4575b4
  theme(legend.position = c(0.2, .85)) 
tree_plot

#Put tip labels on - they will be scientific names with no underscore, and in italics
p1  =ggtree(tree1, aes(color=trait),
            ladderize = FALSE, continuous = TRUE, size=0.5, layout="circular") +
  scale_color_gradientn(legend_title, colours=c("#b2182b", '#fc8d59', '#fee090', '#91bfdb', '#4575b4'))+
  theme(legend.position = c(.01, .85)) 
p1

#Add labels to plot that are spaced and italics
lb = get.tree(tree1)$tip.label
lc= sub("_", " ", lb)
d = data.frame(label=lb, label2 = lc)
a = p1 %<+% d + geom_tiplab(aes(label=label2), size=0.000000000000005, fontface="italic", colour="black") 
tree_plot1 = a + xlim_tree(62)
tree_plot1










## Non-phylogeny plots ####
## Body mass vs. Tb_min
ggplot(torpor_dat, aes(log(BM), Tbmin)) + geom_point(aes(col=Family, shape=Type), size=3) + 
  geom_smooth(aes(col=Type), method = 'lm') +
  facet_grid(.~Class, scales = "free") + my_theme +
  scale_color_manual(values=my_colors7) + xlab(bodymass.lab) + ylab(Tbmin.lab) 

## Body mass vs. Tb_min with family
ggplot(torpor_dat[torpor_dat$Class=="AVES",], aes(log(BM), Tbmin)) + geom_point(aes(col=Order), size=3) + 
  #geom_smooth(aes(col=Type), method = 'lm') +
  facet_grid(.~Class, scales = "free") + my_theme +
  scale_color_viridis_d() +
  #scale_color_manual(values=my_colors7) + 
  xlab(bodymass.lab) + ylab(Tbmin.lab) 

## Body mass vs. MR_min
ggplot(torpor_dat, aes(log(BM), MRmin)) + geom_point(aes(col=Type), size=3) + 
  geom_smooth(aes(col=Type), method = 'lm') +
  facet_wrap(Class~., scales = "free") + my_theme +
  scale_color_manual(values=my_colors7) + xlab(bodymass.lab) + ylab("Min metabolic rate in torpor (ml O2/g/h)") 

## Body mass vs. MR_min colored by tempdiff for birds
ggplot(torpor_dat[torpor_dat$Class=="AVES",], aes(log(BM), MRmin)) + geom_point(aes(col=Tbdiff), size=3) + 
  #geom_smooth(aes(col=Type), method = 'lm') +
  #facet_wrap(Class~., scales = "free") + 
  my_theme +
  #scale_color_manual(values=my_colors7) + 
  xlab(bodymass.lab) + ylab("Min metabolic rate in torpor (ml O2/g/h)") 


## Latitude vs. Tb_min
ggplot(torpor_dat, aes(lat_study_loc, Tbmin)) + geom_point(aes(col=Type), size=3) + 
  geom_smooth(aes(col=Type), method = 'lm') +
  facet_wrap(Class~., scales = "free") + my_theme +
  scale_color_manual(values=my_colors7) + xlab("Latitude") + ylab(Tbmin.lab) 

## Latitude vs. Body mass
ggplot(torpor_dat, aes(lat_study_loc, log(BM))) + geom_point(aes(col=Type), size=3) + 
  geom_smooth(aes(col=Type), method = 'lm') +
  facet_wrap(Class~., scales = "free") + my_theme +
  scale_color_manual(values=my_colors7) + xlab("Latitude") + ylab(bodymass.lab) 


## Test ####
## Traits to add on phylogeny
traits1<-familyData$Portion_of_family_with_at_least_one_full_specimen
names(traits1)<-familyData$Birdlife.Family
traits2<-familyData$Portion_of_family_with_at_least_four_specimens
names(traits2)<-familyData$Birdlife.Family
traits3<-familyData$Portion_of_family_with_at_least_four_main_traits_measured_per_species
names(traits3)<-familyData$Birdlife.Family

## Trying out phylogeny
plotTree.wMultipleBars(familyPhyTipsL_av, #xT1=traits1, xT2=traits2, xT3=traits3,   
                       scale=0.2, widthT1=2.8,widthT2=3.7,widthT3=4.25,
                       ThresholdT1=75,ThresholdT2=75,ThresholdT3=75,
                       method="plotTree",type="fan", tip.labels=TRUE,tip.label.Offset=65,buffer=44,
                       alphaVal=100,colT1="gold",colT2="red4",colT3="dodgerblue4",border=FALSE)


