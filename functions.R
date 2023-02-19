#' @title "functions for the individual project of AAB"
#' @author Miguel Casanova, fc24475
#' @date May 2021
#' @description  The following script, contains a list of functions to be used with the main script for the "Individual Homework" of AAB.

############################################################################################################################
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggtext)

############################################################################################################################
################################################## Fst Functions ###########################################################
############################################################################################################################
#' getFst
#' @author Vitor Sousa
#' 
#' This function was taken from a script provided by Vitor Sousa, for the AAB course.
#' It computes the FST according to the Hudson's estimator following Bathia.
#' 
#' INPUT:
#'  @param geno1  matrix with nsites x nind1 with the genotype data for population 1
#'  @param geno2  matrix with nsites x nind2 with the genotype data for population 2
#' OUTPUT
#'  @return FST between the two populations

getFst <- function(geno1, geno2) {
  
  # compute the sample size for each site
  ss1 <- 2*rowSums(!is.na(geno1))
  ss2 <- 2*rowSums(!is.na(geno2))
  
  # compute allele frequency for each site
  freq1 <- rowSums(geno1, na.rm=T)/ss1
  freq2 <- rowSums(geno2, na.rm=T)/ss2
  
  # compute the terms p1(1-p1) and p2(1-p2)
  p1 <- freq1*(1-freq1)
  p2 <- freq2*(1-freq2)
  
  
  # compute the square of the difference among the allele frequenies
  pdiffsquare <- (freq1-freq2)^2
  
  # compute the numerator
  numerator <- pdiffsquare-(p1/(ss1-1))-(p2/(ss2-1))
  
  # compute the denominator
  denominator <- (freq1*(1-freq2))+(freq2*(1-freq1))
  
  
  # output FST estimators
  sum(numerator, na.rm=T)/sum(denominator, na.rm=T)
}

############################################################################################################################
#' PermutFST
#' 
#' The function 'PermutFST' implements a permutation test to assess if the population Fst from
#' which the samples were taken, is different from zero.
#' The function takes two matrices with genotypes of individuals belonging to two samples populations,
#' as well as a number of permutations (it assumes 1000 permutations as a default parameter). 
#' I've changed the requested output for the function, to have a list of three elements: 
#' 1- The values for all of the simulated pairwise Fst.
#' 2- The value for the observed Fst between the two sampled populations provided.
#' 3- The p-value for the permutation test
#' 
#' INPUT:
#'  @param geno1  matrix with nsites x nind1 with the genotype data for population 1
#'  @param geno2  matrix with nsites x nind2 with the genotype data for population 2
#'  @param perm   number of permutations
#' OUTPUT
#'  @return list with the following elements:
#'  - simulated Fst from n permutations between two populations
#'  - observed Fst between two populations
#'  - p-value for the permutation test
 
PermutFST <- function(geno1, geno2, perm = 1000) {
  
  # Compute the observed FSt. This is our test statistic.
  observedFst <- getFst(geno1, geno2) 
  
  # Fuse the genotype matrices
  fusedGeno <- cbind(geno1, geno2)
  
  # Perform the permutations. 
  permutationNumber <- perm
  simulatedpairwiseFst_permutations <- c()
  
  for (i in 1:permutationNumber) {
    
    fusedGeno_Permutation <- fusedGeno[, sample(ncol(fusedGeno), size = ncol(fusedGeno), replace = FALSE), drop = FALSE]
    simulatedFst <- getFst(fusedGeno_Permutation[, 1:ncol(geno1)], 
                           fusedGeno_Permutation[, (ncol(geno1)+1):ncol(fusedGeno)])
    simulatedpairwiseFst_permutations <- append(simulatedpairwiseFst_permutations, simulatedFst)}
  
  # Output the p-value
  # p-value is given by the proportion of permutations where we obtain a value more extreme than the observed sample statistic
  pvalue <- sum(simulatedpairwiseFst_permutations >= observedFst)/permutationNumber
  
  outputs <- list(simulatedpairwiseFst_permutations, observedFst, pvalue)
  names(outputs) <- c("simFstPerm", "obsFst", "pvalue")
  
  return(outputs)
}  

############################################################################################################################
################################################ Plotting Functions ########################################################
############################################################################################################################

#' sampling_Fst
#' This function takes a list of simulated Fst between two populations, from a given number of permutation tests.
#' It will then plot a histogram with the frequencies of simulated Fst values. The function also takes the observed Fst, 
#' and plots it, in order to compare the simulated Fst, with its observed value.
#' 
#' INPUT:
#'  @param PredFst list with n simulated Fst between two populations (permutation test)
#'  @param ObsFst observed Fst value, for the two sampled populations
#'  @param title title of the plot  
#' OUTPUT
#'  @return plot with the joint distribution and histograms of prior and posterior distribution

sampling_Fst <- function(PredFst, ObsFst, title) {
  
  df = data.frame(PredFst)
  
  ggplot(df, aes(x=PredFst)) +
    theme_bw() +
    geom_histogram(aes(y=..density..), bins=100, color="black", fill = "blue", alpha = 0.3) +
    theme(plot.title=element_text(size=16, face="bold", hjust=0.5)) +
    #geom_density(alpha=.1, fill="red") +
    scale_x_continuous(limits = c(-0.01,0.3)) +
    geom_vline(aes(xintercept = ObsFst), color = "red", linetype = 2, size = 1.01) +
    labs(title = title,
         x = "Simulated pairwise Fst",
         y = "Frequency")
}

############################################################################################################################
#' abc_plot
#' Function to plot the joint distribution of parameters and summary statistics, and show the accepted points after rejection step.
#' it also plots the prior and posterior distributions, as marginal plots.
#' 
#' INPUT:
#'  @param SimStatistics vector of size n, corresponding to n simulation, with the summary statistics for each of this.
#'  @param obsStatistics observed summary statistic of the data. 
#'  @param tol_distance value of the distance threshold, defined by the quantile of distances defined in tol (e.g. if tol is 0.10, this value 
#'  corresponds to the maximum distance of the 10% closest simulations).
#'  @param parameters vector of size n, corresponding to n simulation, with the parameter value used for each simulation.
#'  @param dist_sumstat vector of size n, corresponding to n simulations, with the absolute distance between simulated and observed data for each of this.   
#'  @param title  title of the plot.
#' OUTPUT
#'  @return plot with the joint distribution and histograms of prior and posterior distribution


abc_plot <- function(SimStatistics, obsStatistics, tol_distance, parameters, dist_sumstat, title) {
  
  df <- data.frame(simS = SimStatistics, param = parameters, thresh = (dist_sumstat<tol_distance))
  
  # With this, we can do lots of pretty graphs.
  plt <- ggplot(df, aes(x = simS, y = param, color = thresh)) +
    geom_point(alpha = 0.4, size = 1, aes(color = thresh)) +
    scale_color_manual(values=c("black", "red"),
                       labels = c("Prior", "Posterior\nAccepted")) +
    geom_vline(xintercept = obsStatistics, col = "blue", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c((obsStatistics+tol_distance), (obsStatistics-tol_distance)), col = "darkred", linetype = "dashed", size = 0.51) +
    
    labs(title = title, 
         x = "Simulated Summary Statistics\n(Segregating Sites)",
         y = "Model Parameter\n(Population Effective Size)", 
         color = "Parameters") +
    theme_bw() +
    theme(plot.title=element_text(size=16, 
                                  face="bold", 
                                  hjust=0.5),
          legend.position = c(0.85,0.15),
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) 
  
  plt <-  plt + annotate("rect", xmin = (obsStatistics-tol_distance), xmax = (obsStatistics+tol_distance), ymin = -Inf, ymax = +Inf,
                         alpha = .1,fill = "red")
  
  ggExtra::ggMarginal(plt, 
                      type = "density", 
                      groupColour = T,
                      groupFill = T,
                      alpha = 0.4)
}

############################################################################################################################
#' abc_plot_V2
#' This function was implemented to make ABC plots for joint distribution of parameters and summary statistics, 
#' as well as marginal distributions. It is an alternative version, that uses cowplot and tidyverse,
#' to allow a more accurate representation of marginal plots.
#' 
#' INPUT:
#'  @param SimStatistics vector of size n, corresponding to n simulation, with the summary statistics for each of this.
#'  @param obsStatistics observed summary statistic of the data. 
#'  @param tol_distance value of the distance threshold, defined by the quantile of distances defined in tol (e.g. if tol is 0.10, this value 
#'  corresponds to the maximum distance of the 10% closest simulations).
#'  @param parameters vector of size n, corresponding to n simulation, with the parameter value used for each simulation.
#'  @param dist_sumstat vector of size n, corresponding to n simulations, with the absolute distance between simulated and observed data for each of this.   
#'  @param title  title of the plot.
#' OUTPUT
#'  @return plot with the joint distribution and histograms of prior and posterior distribution

abc_plot_V2 <- function(SimStatistics, obsStatistics, tol_distance, parameters, dist_sumstat, title) {
  
  df <- data.frame(simS = SimStatistics, param = parameters, thresh = (dist_sumstat<tol_distance))
  
  
  g <- ggplot(df, aes(x = simS, y = param, color = thresh)) +
    geom_point(alpha = 0.3, size = 1, aes(color = thresh)) +
    scale_color_manual(values=c("black", "red"),
                       labels = c("Prior", "Posterior\nAccepted")) +
    geom_vline(xintercept = obsStatistics, col = "blue", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c((obsStatistics+tol_distance), (obsStatistics-tol_distance)), col = "darkred", linetype = "dashed", size = 0.51) +
    
    labs(title = title,
         x = "Simulated Summary Statistics\n(Segregating Sites)",
         y = "Model Parameter\n(Population Effective Size)",
         color = "Parameters") +
    theme_bw() +
    theme(plot.title=element_text(size=16, 
                                  face="bold", 
                                  hjust=0.5),
          legend.position = c(0.85,0.15),
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) +
    annotate("rect", xmin = (obsStatistics-tol_distance), xmax = (obsStatistics+tol_distance), ymin = -Inf, ymax = +Inf,
             alpha = .1,fill = "red")
  
  xhist <- ggplot(df, aes(x = simS)) +
    geom_histogram(aes(y=..density..), bins=100, color="black", fill = "black", alpha = 0.3) +
    geom_histogram(data = df[df$thresh==TRUE,], aes(y=..density..), bins=100, color="black", fill = "red", alpha = 0.3) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  yhist <- ggplot(df, aes(x = param)) +
    geom_histogram(aes(y=..density..), bins=100, color="black", fill = "black", alpha = 0.3) +
    geom_histogram(data = df[df$thresh==TRUE,], aes(y=..density..), bins=100, color="black", fill = "red", alpha = 0.3) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    coord_flip()
  
  
  g %>%
    insert_xaxis_grob(xhist, grid::unit(1.4, "cm"), position = "top") %>%
    insert_yaxis_grob(yhist, grid::unit(1.4, "cm"), position = "right") %>%
    ggdraw()
}

############################################################################################################################
#' posterior_plot
#' This function was implemented to make ABC plots for posterior distribution of parameters, after the rejection step.

#' INPUT:
#'  @param parameters vector of size n, corresponding to n simulation, with the parameter value used for each simulation.
#'  @param closets_sims list with the indexes of the simulated statistics that are under a distance threshold. Used to compute the value of
#'  parameter, from the points with a summary statistics closer to the observed statistics.
#'  @param title the title of the plot.  
#' OUTPUT
#'  @return plot the posterior distribution of the accepted parameters as an overlayed histogram and density plot.

posterior_plot <- function(parameters, closest_sims, title) {
  
  dfPosterior = data.frame(posterior = parameters[closest_sims])
  
  ggplot(dfPosterior, aes(x=posterior)) +
    geom_histogram(aes(y=..density..), bins=50, color="black", fill = "red", alpha = 0.3) +
    geom_density(alpha=.4, fill="#FF6666") +
    theme_bw() +
    theme(plot.title=element_text(size=16, face="bold", hjust=0.5)) +
    scale_x_continuous(limits = c(0.0,100000)) +
    labs(title = title,
         x = "Effective Size")
}

############################################################################################################################
#' prior_plot
#' This function was implemented to make ABC plots for the prior distribution of parameters.
#' 
#' INPUT:
#'  @param parameters vector of size n, corresponding to n simulation, with the parameter value used for each simulation.
#'  @param title the title of the plot.  
#' OUTPUT
#'  @return plot the prior distribution of parameters as an histogram and density plots.

prior_plot <- function(parameters, title) {
  
  dfPrior = data.frame(priors = parameters)
  
  ggplot(dfPrior, aes(x=priors)) +
    geom_histogram(aes(y=..density..), bins=50, color="black", fill = "blue", alpha = 0.3) +
    theme_bw() +
    theme(plot.title=element_text(size=16, face="bold", hjust=0.5)) +
    geom_density(alpha=.2, fill="blue") +
    scale_x_continuous(limits = c(0.0,100000)) +
    labs(title = title,
         x = "Effective Size")
}


############################################################################################################################
#' posterior_density_plot
#' Function to plot the density estimation of the posterior distribution of the parameters. It also calculates the mode, plotting it
#' as a line. It also displays the value for the mode, as a text-box.
#' 
#' INPUT:
#'  @param posterior_density object of class "density". Contains information about the values for the parameter, and their respective density.
#'  @param title the title of the plot.  
#' OUTPUT
#'  @return plot with density curve for the posterior parameters, and value for most likely value of the parameter (as the density's mode)

posterior_density_plot <- function(posterior_density, title) {
  
  dfDensity <- data.frame(Posterior = posterior_density$x, Density = posterior_density$y)
  mode <- posterior_density$x[which.max(posterior_density$y)]
  
  annotations <- data.frame(
    xpos = c(Inf),
    ypos =  c(Inf),
    annotateText = paste0("**Mode:**", round(mode, 0)),
    hjustvar = c(1.5),
    vjustvar = c(3))
  
  ggplot(dfDensity, aes(x = Posterior, y = Density)) +
    theme_bw() +
    geom_density(stat = "identity", linetype = "dotted", fill = "mistyrose", alpha = 0.5, size = .5) +
    geom_vline(xintercept = mode, col = "darkred", linetype = "twodash", size = 1) +
    xlim(0, 100000) +
    theme_bw() +
    theme(plot.title=element_text(size=16, face="bold", hjust=0.5))  +
    labs(title = title, 
         x = "Posterior Effective Size",
         y = "Density") +
    geom_richtext(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), fill = "mistyrose")
  
}

############################################################################################################################
#' combined_density_plot
#' Function to plot the combined density estimation of the posterior distribution of the parameters for two populations. 
#' It also calculates the modes, plotting them as dotted lines. Finally, it also displays the value for the modes, as text-boxes.
#' 
#' INPUT:
#'  @param posterior_density_population_1 object of class "density" for population 1. Contains information about the values for the parameter, and their respective density.
#'  @param posterior_density_population_2 object of class "density" for population 2. Contains information about the values for the parameter, and their respective density.
#'  @param name_population_1 string with the name of population 1.
#'  @param name_population_2 string with the name of population 2.
#'  @param title the title of the plot.  
#' OUTPUT
#'  @return plot with density curves for the posterior parameters for two populations, 
#'  and value for most likely value of the parameter (as the density's mode).

combined_density_plot <- function(posterior_density_population_1, posterior_density_population_2, name_population_1, name_population_2, title) {
  
  dfDensityPop1 <- data.frame(Posterior = posterior_density_population_1$x, Density = posterior_density_population_1$y, Population = name_population_1)
  dfDensityPop2 <- data.frame(Posterior = posterior_density_population_2$x, Density = posterior_density_population_2$y, Population = name_population_2)
  dfDensityMerged <- rbind(dfDensityPop1, dfDensityPop2)
  
  modePop1 <- posterior_density_population_1$x[which.max(posterior_density_population_1$y)]
  modePop2 <- posterior_density_population_2$x[which.max(posterior_density_population_2$y)]
  
  annotations <- data.frame(
    xpos = c(Inf),
    ypos =  c(Inf),
    annotateText = c(paste0(name_population_1, ":", round(modePop1, 0)), paste0(name_population_2, ":", round(modePop2, 0))),
    hjustvar = c(1.5),
    vjustvar = c(3,5),
    Population = c(name_population_1, name_population_2))
  
  
  ggplot(dfDensityMerged, aes(x = Posterior, y = Density, fill = Population)) +
    theme_bw() +
    geom_density(stat = "identity", linetype = "dotted", alpha = 0.4, size = .5, color = "black") +
    scale_fill_manual(values=c("mistyrose", "lightskyblue")) +
    geom_vline(xintercept = modePop1, col = "darkred", linetype = "twodash", size = .8) +
    geom_vline(xintercept = modePop2, col = "darkblue", linetype = "twodash", size = .8) +
    xlim(0, 100000) +
    theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
          legend.position = "none")  +
    labs(title = title, 
         x = "Posterior Effective Size",
         y = "Density") +
    geom_richtext(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black")
  
}

############################################################################################################################
#' combined_density_plot_X3
#' Function to plot the combined density estimation of the posterior distribution of the parameters for three posterior distributions. 
#' It also calculates the modes, plotting them as dotted lines. Finally, it also displays the value for the modes, as text-boxes.
#' 
#' INPUT:
#'  @param posterior_density_1 object of class "density" for population 1. Contains information about the values for the parameter, and their respective density.
#'  @param posterior_density_2 object of class "density" for population 2. Contains information about the values for the parameter, and their respective density.
#'  @param posterior_density_3 object of class "density" for population 3. Contains information about the values for the parameter, and their respective density.
#'  @param name_1 string with the name of population 1.
#'  @param name_2 string with the name of population 2.
#'  @param name_3 string with the name of population 3.
#'  @param title the title of the plot.  
#' OUTPUT
#'  @return plot with density curves for the posterior parameters for three posterior density distributions, 
#'  and values for most likely values of the parameter (as the density's mode).

combined_density_plot_X3 <- function(posterior_density_1, posterior_density_2, posterior_density_3, name_1, name_2, name_3, title) {
  
  dfDensity_1 <- data.frame(Posterior = posterior_density_1$x, Density = posterior_density_1$y, Tolerance = name_1)
  dfDensity_2 <- data.frame(Posterior = posterior_density_2$x, Density = posterior_density_2$y, Tolerance = name_2)
  dfDensity_3 <- data.frame(Posterior = posterior_density_3$x, Density = posterior_density_3$y, Tolerance = name_3)
  
  dfDensityMerged <- rbind(dfDensity_1, dfDensity_2, dfDensity_3)
  
  mode_1 <- posterior_density_1$x[which.max(posterior_density_1$y)]
  mode_2 <- posterior_density_2$x[which.max(posterior_density_2$y)]
  mode_3 <- posterior_density_3$x[which.max(posterior_density_3$y)]
  
  annotations <- data.frame(
    xpos = c(Inf),
    ypos =  c(Inf),
    annotateText = c(paste0(name_1, ":", round(mode_1, 0)), 
                     paste0(name_2, ":", round(mode_2, 0)), 
                     paste0(name_3, ":", round(mode_3, 0))),
    hjustvar = c(1.5),
    vjustvar = c(4,6, 8),
    width =  1.5,
    height = 1.5,
    Tolerance = c(name_1, name_2, name_3))
  
  
  ggplot(dfDensityMerged, aes(x = Posterior, y = Density, fill = Tolerance)) +
    theme_bw() +
    geom_density(stat = "identity", linetype = "dotted", alpha = 0.1, size = .6, color = "black") +
    geom_vline(xintercept = mode_1, col = "darkred", linetype = "twodash", size = .4) +
    geom_vline(xintercept = mode_2, col = "darkgreen", linetype = "twodash", size = .4) +
    geom_vline(xintercept = mode_3, col = "darkblue", linetype = "twodash", size = .4) +
    xlim(0, 100000) +
    theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
          legend.position = "none")  +
    labs(title = title, 
         x = "Posterior Effective Size",
         y = "Density") +
    geom_textbox(data=annotations,
                 aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar, label=annotateText, fill = Tolerance),
                 alpha = 0.6, 
                 width = unit(0.2, "npc"), 
                 color = "black")
  
}
