#Load packages

packageVersion("phytools")

####################
####################
####################
####################Plot species_iridescence, male_iridescence, and female_iridescence
####################
####################
####################
#Should I use phylosignal
install.packages("phylosignal")
library(adephylo)
library(phylobase)
library(phylosignal)
p4d <- phylo4d(ir_tree,color$species_iridescence)

#should I use diversitree?
library(diversitree)
trait.plot(ir_tree, color, cols =list(color$Species_Iridescence ("pink","black")))
trait.plot(ir_tree, color, cols = list(color$species_iridescence = c("white", "black"), color$female_iridescence = c("white", "red"), color$male_iridescence = c("white","blue")))

#should I use phytools?
fmode<- as.factor(setNames(ir_dataframe[,1],rownames(ir_dataframe)))
dotTree(ir_tree,fmode,colors=setNames(c("blue","red"),
                                      c("1","0")),ftype="i",fsize=0.7)
