
#Iridescence workflow

####Remove below section - used to edit askpass error
use_git_config(user.name = "map469", user.email = "map469@cornell.edu")
install.packages("gitcreds")
library(gitcreds) # install.packages("gitcreds")
gitcreds_set()
####Remove above section for final data analyses - used to edit askpass error

#Read in csv file
color <- read.csv("12_01_2022_irid.csv")

#Check what the database structure is 
str(color)
#check to make sure that you have BirdTree_Species appears first in the structure
#check to be sure all NA species are removed from your dataset
#


#Convert relevant columns to factors. 
##NOTE: This is only the first three of many columns. Others still need adding.
##NOTE: Currently most columns have text or missing values that need changing. 

color$species_iridescence <- factor(color$Species_Iridescence)
color$male_iridescence <- factor(color$Male_Iridescence)
color$female_iridescence <- factor(color$Female_Iridescence)
color$body_feathers <- factor(color$Body_Feathers)
color$body_male <- factor(color$Body_Male)
color$body_female <-factor(color$Body_Female)
color$flight_feathers <- factor(color$Flight_Feathers)
color$flight_male <- factor(color$Flight_Male)
color$flight_female <-factor(color$Flight_Female)
color$nonflight_feathers <- factor(color$NonFlight_Feathers)
color$nonflight_male <- factor(color$NonFlight_Male)
color$nonflight_female <-factor(color$NonFlight_Female)
color$dorsal_feathers <- factor(color$Dorsal_Feathers)
color$dorsal_male <- factor(color$Dorsal_Male)
color$dorsal_female <-factor(color$Dorsal_Female)
color$primary.lifestyle <- factor(color$Primary.Lifestyle)
color$habitat <- factor(color$Habitat)

#Load packages for analysis
library(ape)
library(MCMCglmm)
library(geiger)
library(phytools)

#Read in the tree
ir_names <- read.csv("ir_names.csv") #This is a file with JETZ Tree names from the database 
ir_dataframe <- data.frame(ir_names)
rownames(ir_dataframe) <- ir_names [,1] 
ir_tree <- read.tree ("BirdTree_Jetz.tre") #This is the tree file from Eliot
name.check (ir_tree, ir_dataframe) #Identifies mistmatched species
#$tree_not_data are species that are in the tree and not your dataset
#$data_not_tree are species that are in your data and don't appear in tree
#We can reconcile this by dropping the unneccessary tips


#check for overlap between the tree and dataframe, and drop tips on the tree
ir_overlap <- name.check(ir_tree, ir_dataframe)
ir_tree_drop <- drop.tip(ir_tree,ir_overlap$tree_not_data)
ir_tree_drop
name.check (ir_tree_drop, ir_dataframe)
#Should output "OK" if all names match up. 
#If problem persists, name check the species of concern


#Compute branch lengths of the tree
tree <- compute.brlen(ir_tree_drop, method = "Grafen")
inv.tree <- inverseA (tree, nodes = "ALL", scale = TRUE)

#Set priors for the MCMCglmm model. 
##Note: these are relatively uninformative priors. I suggest keeping them for your initial analysis. 
##Before finalizing the analysis it would be good to check how much changing the priors affects results (e.g., change to v = 1, nu = 2) 
prior1 <- list (G=list (G1=list (V=1,nu=0.002)), 
              R=list(V=1,nu=0.002))
#In this section of code the number of "G-lists" you have equates with the number of random effects in your model. 


#Run the model 
##Note: this is code copied from a previous analysis. You'll need to change the names to your column and tree names. 
##Note: I find that for trial runs using nitt=50000, burnin=5000 and thin=200 works well, and doesn't take too long to run. 
##Note: For your final analyses you'll want to change nitt to 1,000,000.

#Sample model looking at the centroid of latitude as a predictor of iridescence in flight feathers. 
model_lat <- MCMCglmm(flight_feathers~1 + 
                        Centroid.Latitude,
                random = ~species,
                family= "categorical",
                ginverse=list(species=inv.tree$Ainv), 
                prior=prior1,
                data=color,
                nitt=50000,#specifies how many times the simulation is performed (ORIGINALLY 100k but changed for practice run)
                burnin=5000,#specifies how many simulations are discarded at the beginning
                thin=200) #specifies that every nth iteration is saved. 
#"gaussian", "poisson", "categorical", "multinomial", "ordinal", "threshold", 
# "exponential", "geometric", "cengaussian", "cenpoisson", "cenexponential", "zipoisson", 
#"zapoisson", "ztpoisson", "hupoisson", "zibinomial", "threshold", "nzbinom" , "ncst", "msst" , 
#"hubinomial", "ztmb" and "ztmultinomial"
?MCMCglmm
summary(model_lat)

#Make sure to rerun with more informative priors (e.g., v=1, nu=1)


##Here are some ways you can (and should!) explore your model:

# 1) Calculate posterior probability of phylogenetic signal (lambda) 
    # Units = residual variance, built into code already.
lambda <- model_lat$VCV[,'species']/(model_lat$VCV[,'species']+model_lat$VCV[,'units'])

# calculate the posterior mean (mean of the posterior distribution), 
# posterior mode (most likely value regarding the posterior distribution) and # the 95% credible interval of lambda
mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)

# 2) Create Density and trace plots so that you can visually inspect them for independence and consistency of the posterior distribution. 
  ##Note: I had some trouble with this when re-running code in May 2022, error message "figure margins too large"
  ##Update: error message continues "figure margins too large"
plot(model_lat$VCV)

# 3) Check Autocorrelation:
  ##Having everything below 0.05 is good. Ignore the first row (Lag 0) when examining your results. 
  ##Note: trial models run with a low number of iterations will probably have high autocorrelation. So this doesn't really make sense to do until you run "final" models with close to 1M iterations. 
autocorr(model_lat$VCV)

####################
####################  Water Repellency
####################  Hypothesis 1: Iridescence is predicted by primary.lifestyle - aquatic
####################
####################

#Model examines primary.lifestyle as predictor of iridescence in flight feathers. 
??MCMCglmm

prior1 <- list (G=list (G1=list (V=1,nu=0.002)), 
                R=list(V=1,nu=0.002))

prior2 <- list (G=list (G1=list (V=1,nu=2)), 
                R=list(V=1,nu=0.002))

model_wet <- MCMCglmm(Species_Iridescence~1 + 
                        Primary.Lifestyle,
                      random = ~species,
                      family="binomial",
                      ginverse=list(species=inv.tree$Ainv), 
                      prior=prior1,
                      data=color,
                      nitt=50000,#specifies how many times the simulation is performed (ORIGINALLY 100k but changed for practice run)
                      burnin=5000,#specifies how many simulations are discarded at the beginning
                      thin=200) #specifies that every nth iteration is saved. 

summary(model_wet)

model_wet2 <- MCMCglmm(Species_Iridescence~1 + 
                        Primary.Lifestyle,
                      random = ~Family,
                      family="gaussian",
                      ginverse=list(species=inv.tree$Ainv), 
                      prior=prior1,
                      data=color,
                      nitt=50000,#specifies how many times the simulation is performed (ORIGINALLY 100k but changed for practice run)
                      burnin=5000,#specifies how many simulations are discarded at the beginning
                      thin=200) #specifies that every nth iteration is saved. 

summary(model_wet)

