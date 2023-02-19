# Aplicações Avançadas em Biologia

# Masters' Degree in Bioinformatics and Computational Biology
# FCUL - 2020/2021jass
# May 2021
# author: Miguel Casanova, fc24475

############################################################################################################################
########################## Individual exercises - Fst and permutation tests ################################################
############################################################################################################################
# FST and permutations

# Fst is a summary statistic that is used to measure allele frequency similarity, between populations. 
# The Fst index can give us an idea about the population structure.
# It is a function of allele frequencies within each pair of populations and the overal allele frequencies across all populations.
# There are several ways to estimate Fst. Here we use the Hudson estimator (^FSt = 1 - [^He(W) within/ ^He(B) between])

# The FSt can vary between 0 and 1. A FSt of 0, means that there is a high flux of genes between the two populations.
# An FSt of 1, means that there is a very reduced genetic flux between the populations. 
# As such, this means that the populations are very "differentiated".

############################################################################################################################
####################################### Loading the genotypic data #########################################################
############################################################################################################################

# Get and set the working directory.
getwd()
setwd("C:/Cloud/GoogleDrive/Mestrado_BCG/4th_Semester/AplicacoesAvancadasBiologia")

# We next will load into memory, the functions contained within the file 'functions.R'.
# These functions, include a function to calculate Fst between population (from the classes of AAB),
# and a function I have made to implement permutation tests and test whether Fst between populations is different from zero.
# It also includes functions to make several plots for the project.
source("Individual_Homework/functions.R")

# We need to put all the data from a pair of populations in the same vector 
# as we are assuming they are from the same population.

# Let's load the genotype matrix from Henn et al.
genotypes <- as.matrix(read.table("Day2/Practical_HumanData/Henn_et_al_data/Hennetal_genotypeMatrix.geno", 
                                  header=T,
                                  stringsAsFactors = F, 
                                  na.strings = "NA"))
str(genotypes)

# Let's now create a metadata file, with the population IDs of each individual
metadata <- data.frame(Individuals = colnames(genotypes), PopID = colnames(genotypes))
# To create the population metadata information, we will profit from the fact that all individuals 
# are identified with the population name, followed by an underscore and 9 digits. 
# As such, we will remove the 10 last digits from each individual ID, to get their populations.
metadata$PopID = substr(metadata$PopID, 1, nchar(metadata$PopID)-10)

# What are all the unique populations in the matrix?
populationNames <- unique(metadata$PopID)
populationNames

# Using the above metadata information, we can easily load the list of individuals belonging to a given population and get their number.
# Let's pick two populations.
# For Maya:
mayaGenotypes <- genotypes[, metadata$PopID == "Maya"]
numberMaya <- ncol(mayaGenotypes)
numberMaya

# For San:
sanGenotypes <- genotypes[, metadata$PopID == "San"]
numberSan <- ncol(sanGenotypes)
numberSan

# We can now fuse the two genotype matrixes, for performing our permutations, down the line.
mayaSanGenotypes <- cbind(mayaGenotypes, sanGenotypes)
totalNumberIndividuals <- ncol(mayaSanGenotypes)
totalNumberIndividuals

############################################################################################################################
########################### Exercise 1 - Creating and Running the function #################################################
############################################################################################################################
# The function for exercise 1 is contained within the 'functions.R' script (where it is properly commented).
# Succinctly, the function 'PermutFST' takes two matrices with genotypes of individuals belonging to two samples populations,
# as well as a number of permutations (it assumes 1000 permutations as a default parameter). 
# I've changed the requested output for the function, to have a list of three elements: 
# 1- The values for all of the simulated pairwise Fst.
# 2- The value for the observed Fst between the two sampled populations provided.
# 3- The p-value for the permutation test

# As such, using the above genotype matrices produced above, we can easily run the function by doing:
test <- PermutFST(mayaGenotypes, sanGenotypes, 10000)
test

############################################################################################################################
######################## Exercise 2 - Calculating pvalues for pairwise Fst values ##########################################
############################################################################################################################
# To assess the significance of the pairwise Fst values between all pairwise comparison for the dataset of Henn et al. (2015),
# I got inspiration from the loop used in the script from class 2 to calculate the pairwise Fst between all populations.
# I repurposed the loop to calculate the significance (p-value) of the pairwise Fst values.

# Let's first create a matrix that will be used to fill with the pairwise p-values for the permutation test between each pair of populations
# I will start by creating a matrix inputing the population names for each row and columns
pairPvalues <- matrix(NA, ncol=length(populationNames), nrow=length(populationNames), dimnames = list(populationNames, populationNames))

# I'll next create a loop to go through each pair of populations (using the metadata file to help us slice our genotype matrices) 
# and run the 'PermutFST' function, which will calculate the significance of the pairwise Fst values.
for(i in 1:(length(populationNames)-1)) {
  for(j in (i+1):length(populationNames)) {
    # call the function to compute the pairwise p-value for the Fst permutation tests between each pair of populations
    pval <- PermutFST(genotypes[, metadata$PopID == populationNames[[i]]], genotypes[, metadata$PopID == populationNames[[j]]], 1000)
    pairPvalues[i,j] <-  pval$pvalue       
  }
}

# print the matrix with the pairwise FST p-values
View(as.matrix(pairPvalues))

############################################################################################################################
# Let's visualize the pairwise FST p-values graphically
# I will use ggplot2 to plot the pairwise matrix for the pairwise FST p-value.
library(reshape2)

longData<-melt(pairPvalues)

ggplot(longData, aes(Var2, Var1, col = value, fill = value, label = value)) +
  geom_tile() +
  geom_text(col = "black") +
  theme_minimal() +
  scale_fill_gradient2(na.value = "gray95", low = "white", mid = "seashell", high = "indianred1") +
  scale_color_gradient2(na.value = "gray95", low = "white", mid = "seashell", high = "indianred1") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "p-values for pairwise FST") +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5),
        axis.text = element_text(size=10, face="italic"))


############################################################################################################################
################# Exercise 3 - Histogram of the sampling distribution of simulated Fst #####################################
############################################################################################################################
# To perform the plots, I'll use the 'ggplot2' library.

# ggplot2, works with dataframes. As such, our data needs to be  transformed into a dataframe first.This is all done under the hood,
# by a plotting function.
# As an output for the 'PermutFST' function I've created above, we have as 1st and 2nd elements, the values for all simulated Fst, 
# and the observed Fst between the two sampled populations. As such, it is a breeze to create the histograms required.  
sampling_Fst(PredFst = test$simFstPerm,
             ObsFst = test$obsFst,
             title = "Distribution of simulated Fst\n10000 permutations\nMaya Vs San")

############################################################################################################################
######################################## Individual exercises - ABC methods ################################################
############################################################################################################################
# Approximate Bayesian Computation (ABC), is a computational method based on Bayesian statistics, and that provides an 
# approximation of the posterior probability, when a likelihood function is not available by using simulations.
# It replaces data, with summary statistics, and allows inferring about complex models. This is possible by making several 
# assumptions and approximations, that might not always reflect the reality of the data. As such, careful consideration
# has to be taken when using ABC.

# ABC broadly follows the following steps:
# - Define a model to describe your data and sample parameter values from a specified prior distributions;
# - Compute the summary statistics for the simulated data;
# - Compute distance between observed and simulated data.
# - Accept parameters of simulations close to the observed data;

############################################################################################################################
########################################## Loading the data ################################################################
############################################################################################################################
getwd()
setwd("C:/Cloud/GoogleDrive/Mestrado_BCG/4th_Semester/AplicacoesAvancadasBiologia")

# load the package that simulated genetic data according to a given model
# load the functions in coalfunctions.r
source("Day4/ABC/coalfunctions.r")
# load functions in 'functions.r'.
source("Individual_Homework/functions.R")

# Load and visualize the table with the information about the two sampled populations of chimpanzees.
ChimpsSegSites <- read.table("Individual_Homework/SegregatingSitesChimp.txt", 
                             header=T,
                             stringsAsFactors = F,
                             na.strings = "NA")
ChimpsSegSites
############################################################################################################################
########################################### Exercise 4 - ABC methods #######################################################
############################################################################################################################
# Let's start by storing all the information about our data, about our samples populations.

# Observed sumstat - The summary statistics. The summary statistics represents the data.
# In this case, the summary statistics represents the segregating sites that were observed in the 
# two different populations of Central and Western Chimpanzees.
obs_central <- ChimpsSegSites$SegregatingSites[1]
obs_central

obs_western <- ChimpsSegSites$SegregatingSites[2]
obs_western

# All of the remaining parameters are similar, between the two populations.
# Let's define the number of simulations that will be used to sample the posterior.
# In this case, I will use 10000 simulations.
nsim <- 10000

# Next, I'll define the tolerance (closest 10% of simulations) to accept a summary statistic.
# In this case, I'll accept the 10% simulated summary statistics that are the closest to our observed statistics.
tol <- 0.1 

# I'll next define the sample size, taking into consideration that chimps are diploids and, as such, the 
# sample size is multiplied by 2 (two sets of chromosomes)
n <- ChimpsSegSites$NumIndividuals[1]*2 
n

# I'll next define the total mutation per site,
mutrate <- ChimpsSegSites$MutRatePerSite[1]  

# The number of sites in the locus (loci) under study,
nsites <- ChimpsSegSites$NumSites[1] 

# The number of loci under study. In this case, I'm unsure about the number of loci studies, 
# as there isn't any relevant information provided. Nevertheless, we are assuming this number is similar between populations, 
# and as such, shouldn't contribute to differences between the groups.
nloci <- 1 

############################################################################################################################
############################################################################################################################

# 1. Define the prior for Ne, and sample nsim values at random from the prior 
# We start by defining a model (uniform distribution) and sample parameter values 
# from this distribution
param <- runif(nsim, 10, 100000) 

# In this case, the sampled parameter is the effective size of a population.
# The effective size is, roughly, the number of individuals that reproduce. 

# We can plot the above, with ggplot2, using a function defined in 'functions.R'
prior_plot(parameters = param, 
           title = "Prior Distribution of Effective Population Number\nUniform Distribution")

############################################################################################################################
############################################################################################################################
# 2. Simulate the segregating sites - Compute the summary statistics for the simulated data
# Next, we use the data we simulated above and, using the function 'sim.tree.mut', 
# we use the information about the population and the simulated Ne, to calculate a simulated
# summary statistics.

# I start by creating an empty array to save the simulated segregating sites.
simS <- numeric(nsim) 

# I then call a function to simulate genetic data within the for loop
# and for each simulation compute the summary statistic (number of segregating sites)
for(i in 1:nsim) {
  # 2.1. call function to simulate data (coalescent simulator)
  muttree <- sim.tree.mut(sample = n, 
                          current = param[i], 
                          ancestral = param[i],
                          time = 0,
                          nrep = nloci, 
                          mu = mutrate, 
                          L = nsites)
  
  # 2.2. get the number of segregating sites from simulated data
  # this is given by the number of the positions in the resulting matrix (assuming an infinite sites model)
  simS[i] <- ncol(muttree$seg_sites[[1]])
}
############################################################################################################################
############################################################################################################################
# 3. Compute the distance between simulated and observed summary statistics
dist_sumstat_central <- abs(simS - obs_central)
dist_sumstat_western <- abs(simS - obs_western)

# 3.1. Define the tolerance and reject all params with distance larger than tol_dst
tol_dst_central <- quantile(dist_sumstat_central, tol)
tol_dst_western <- quantile(dist_sumstat_western, tol)

# get the index of simulations that are closer to observed data
closest_sims_central <- which(dist_sumstat_central<tol_dst_central)
closest_sims_western <- which(dist_sumstat_western<tol_dst_western)

############################################################################################################################
############################################################################################################################
# 4. Visualize the joint distribution of the summary statistics and parameters
#    and the area of parameters points accepted

# To plot the joint distributions and the posteriors, I used a combination of ggplot, tidyverse and cowplot.
# I've written a couple of functions for plotting, which can be found in 'functions.R'.

abc_plot_V2(SimStatistics = simS,
            obsStatistics = obs_central,
            tol_distance =  tol_dst_central,
            parameters = param,
            dist_sumstat = dist_sumstat_central,
            title = "Joint distribution of Segregating Sites\nand Effective Population\n-Central chimpanzees-")

abc_plot_V2(SimStatistics = simS,
            obsStatistics = obs_western,
            tol_distance = tol_dst_western,
            parameters = param,
            dist_sumstat = dist_sumstat_western,
            title = "Joint distribution of Segregating Sites\nand Effective Population\n-Western chimpanzees-")

############################################################################################################################
############################################################################################################################
# 5. Plot the posterior
# Similarly, I have written a function to plot the posterior probability, using ggplot.
posterior_plot(parameters = param, 
               closest_sims = closest_sims_central,
               title = "Posterior Distribution of Effective Population Number\nCentral chimpanzees")

posterior_plot(parameters = param, 
               closest_sims = closest_sims_western,
               title = "Posterior Distribution of Effective Population Number\nWestern chimpanzees")

############################################################################################################################
############################################################################################################################
# 6. Get the summary of the posterior distribution
# Several post-processing steps can be taken, before computing the posterior distribution.
# In our case, this will be skipped and, as such, we can directly explore our posterior distribution for our 
# parameter: Effective Population Size (Ne).

# Let's get some information about our posterior distribution.
mean(param[closest_sims_central])
quantile(param[closest_sims_central], c(0.025,0.975))

# To calculate the most likely Ne, I'll calculate the density estimation (using function 'density') for the posterior
# and calculate the mode of the density distribution. This will give me the most likely Ne, given the posterior distribution. 

# Let's start by calculating the density estimation with the "density" function
# The density estimates for our posterior include the x (in this case, accepted posteriors) and y values (densities).
posterior_central <- density(param[closest_sims_central])
posterior_central

posterior_western <- density(param[closest_sims_western])
posterior_western

# Using the above, it's very easy to calculate the mode and, as such, the most likely Ne.
posterior_central$x[posterior_central$y==max(posterior_central$y)]

posterior_western$x[posterior_western$y==max(posterior_western$y)]

# To make this more graphical and pretty, we can plot the distribution function and it's mode.
# For this, I use 'ggplot' and 'ggtext', and a function that I wrote and can be found in 'functions.R'.
posterior_density_plot(posterior_density = posterior_central,
                       title = "Density Curve of Posterior Effective Sizes\nCentral chimpanzees")

posterior_density_plot(posterior_density = posterior_western,
                       title = "Density Curve of Posterior Effective Sizes\nWestern chimpanzees")

# We can take this a step further, and plot the combined density distributions. For this, I created a function that
# takes the information about the density distributions of both populations, plotting them together and adding
# the value for the mode, for each of these.
combined_density_plot(posterior_density_population_1 = posterior_central, 
                      posterior_density_population_2 = posterior_western,
                      name_population_1 = "Central",
                      name_population_2 = "Western",
                      title = "Density Curves of Posterior Effective Sizes\n Central and Western chimpanzees")

############################################################################################################################
################################################ Effect of tolerance #######################################################
############################################################################################################################
# To test the effect of tolerance on the posterior, we can easily change the tolerance level.
# Let's change the tolerance level to 1% and 20%, and see how this affects the posterior distribution of the effective number.
# I'll first store the density distribution we have just obtained, with a explicit name.
posterior_central_tol10 <- posterior_central
posterior_western_tol10 <- posterior_western

############################################################################################################################
############################################################################################################################
# We now change the tolerance, and run all the code to calculate the density distribution of posterior Ne.
tol <- 0.01 

# We will start at step 3.1: Define the tolerance and reject all params with distance larger than tol_dst
tol_dst_central <- quantile(dist_sumstat_central, tol)
tol_dst_western <- quantile(dist_sumstat_western, tol)

# get the index of simulations that are closer to observed data
closest_sims_central <- which(dist_sumstat_central<tol_dst_central)
closest_sims_western <- which(dist_sumstat_western<tol_dst_western)

# And finally, we calculate the density functions:
posterior_central_tol1 <- density(param[closest_sims_central])
posterior_central_tol1

posterior_western_tol1 <- density(param[closest_sims_western])
posterior_western_tol1

############################################################################################################################
############################################################################################################################
# Let's do the same, for tolerance of 20%
tol <- 0.2 

# We will start at step 3.1: Define the tolerance and reject all params with distance larger than tol_dst
tol_dst_central <- quantile(dist_sumstat_central, tol)
tol_dst_western <- quantile(dist_sumstat_western, tol)

# get the index of simulations that are closer to observed data
closest_sims_central <- which(dist_sumstat_central<tol_dst_central)
closest_sims_western <- which(dist_sumstat_western<tol_dst_western)

# And finally, we calculate the density functions:
posterior_central_tol20 <- density(param[closest_sims_central])
posterior_central_tol20

posterior_western_tol20 <- density(param[closest_sims_western])
posterior_western_tol20

############################################################################################################################
############################################################################################################################
# Now we only have to compare what happens to the posterior, when different tolerance levels are considered for the rejection step.
# For this, I wrote a function using 'ggplot' and 'ggextra', to plot the density distributions of the posteriors for the three
# tolerance levels. This plot will allow me to visualize the overlaid density curves for the posteriors, as well as their more likely
# effective size. 

# Let's do this for the Central chimpanzees' population
combined_density_plot_X3(posterior_density_1 = posterior_central_tol1,
                         posterior_density_2 = posterior_central_tol10,
                         posterior_density_3 = posterior_central_tol20,
                         name_1 = "Tol 1%",
                         name_2 = "Tol 10%",
                         name_3 = "Tol 20%",
                         title = "Density Curves of Posterior Effective Sizes\nfor different tolerance levels\nCentral chimpanzees")

# And finally, let's do it for the Western chimpanzees' population
combined_density_plot_X3(posterior_density_1 = posterior_western_tol1,
                         posterior_density_2 = posterior_western_tol10,
                         posterior_density_3 = posterior_western_tol20,
                         name_1 = "Tol 1%",
                         name_2 = "Tol 10%",
                         name_3 = "Tol 20%",
                         title = "Density Curves of Posterior Effective Sizes\nfor different tolerance levels\nWestern chimpanzees")

