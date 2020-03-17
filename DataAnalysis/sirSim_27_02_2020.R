#!/usr/bin/env Rscript

# Script written by G. V. DiRenzo, T. S. Tunstall, & L. M. Browne
# Simulation script


# How to calculate proportion of mortality? How to deal with un-whole numbers and negative proportions
# Setting parameter values for abiotic resevior

# Load packages -----------------------------------------------------------
  library(deSolve)
  library(gdata)
  library(synchrony)

# Read in arguments -------------------------------------------------------

  # first argument is number of interacting species
  # Second argument is number of simulations
  args = commandArgs(trailingOnly=TRUE)
  cat("Arguments are:", args, "\n")
  
  # Set number of interactions amongst families
  interacting_families = as.numeric(args[1])
  
  # Number of simulations
  nsim <- as.numeric(args[2])
  
  # Whether in super spreader mode or not
  super_spreader <- as.logical(args[3])

  # Transmission rate in abiotic scenario
  abiotic_trans_rate <- as.numeric(args[4])
  
  # Whether abiotic reservoir or not
  abiotic <- ifelse(abiotic_trans_rate == 0, yes = FALSE, no = TRUE)
  
  # Simulation ID
  sim_id <- as.numeric(args[5])
  
  ## For testing
      # interacting_families <- 2
      # nsim <- 10
      # super_spreader <- FALSE
      # abiotic <- TRUE
      # abiotic_trans_rate <- 0.1
      # sim_id <- 1
  
  # Verbose - printing out things?
  verbose = FALSE
  
  # Plot output for figure 1?
  plot = FALSE
  
  ## Number of focal families
  nsp <- 10
  
  # Set number of days in simulation
  n_days <- 90
  times <- seq(0, n_days, by = 1) 
  
  ### parameters for transmission
  ###### Mortality
  mu_mort <- 1e-1
  K_mort <-  1e1
  theta_mort <- mu_mort/K_mort
  var_mort <- K_mort * theta_mort^2
  sha_mort <- mu_mort/theta_mort
  sca_mort <- mu_mort/K_mort
  
  ###### Inter-taxa transmission
  mu_inter <- 1e-10
  K_inter <-  10
  theta_inter <- mu_inter/K_inter
  var_inter <- K_inter * theta_inter^2
  sha_inter <- mu_inter/theta_inter
  sca_inter <- mu_inter/K_inter
  
  ##### Intra-taxa transmission
  mu_intra <- 1e-2
  K_intra <-  10
  theta_intra <- mu_intra/K_intra
  var_intra <- K_intra * theta_intra^2
  sha_intra <- mu_intra/theta_intra
  sca_intra <- mu_intra/K_intra
  
  ###### Abitoic reseroir
  # Time to infection by abiotic reservoir
  # Each value is the day of infection for families 2-10
  # Each family gets infected a number of times according to some distribution
  mu_inf <- abiotic_trans_rate
  K_inf <- 10
  theta_inf <- mu_inf/K_inf
  var_inf <- K_inf * theta_inf^2
  sha_inf <- mu_inf/theta_inf
  sca_inf <- mu_inf/K_inf
  
  # Matrix to store output
  res1 <- matrix(NA, nrow = nsim, ncol = 12)
  

# Set up equations --------------------------------------------------------

  
  if(abiotic == FALSE){
  # No abiotic reservoir
    Sv <- c()
    Iv <- c()
    for(r in 1:nsp){   # For r in 1 to n species
        v <- c()       # Create an empty variable v
        for(j in 1:nsp){  # for j in 1 to n species
            b <- paste0("beta", r, j) # Paste betarj
            s <- paste0("S", r) # paste Sr
            i <- paste0("I", j) # paste Ij
            v <- c(v, paste(b, s, i, sep="*")) # c(betarj*Sr*Ij, betarj*Sr*Ij, betarj*Sr*Ij, ...)
        }
        eq <- paste(v, collapse=" + ") # "betarj*Sr*Ij + betarj*Sr*Ij + betarj*Sr*Ij+ ..."
        paste0(paste0("dS",r), " <- ", "-(", eq, ")") # dSr <- -(betarj*Sr*Ij + betarj*Sr*Ij + ...)
        dSexp <- paste0(paste0("dS",r), " <- ", "-(", eq, ")") # dSexp <- dSr <- -(betarj*Sr*Ij + betarj*Sr*Ij + ...)
        dIexp <- paste0(paste0("dI",r), " <- ", eq, " - ", paste0("g",r), "*", paste0("I",r)) # dIexp <- dIr <- betarj*Sr*Ij + betarj*Sr*Ij + ... - g * Ir
        Sv <- c(Sv, dSexp) 
        Iv <- c(Iv, dIexp)
    }
# Add in abiotic reservoir dynamics ---------------------------------------
    } else if(abiotic == TRUE){
    Sv <- c()
    Iv <- c()
    for(r in 1:nsp){   # For r in 1 to n species
      v <- c()       # Create an empty variable v
      for(j in 1:nsp){  # for j in 1 to n species
        b <- paste0("beta", r, j) # Paste betarj
        s <- paste0("S", r) # paste Sr
        i <- paste0("I", j) # paste Ij
        d <- paste0("delta", r)
        v <- c(v, paste(b, s, i, sep="*")) # c(betarj*Sr*Ij, betarj*Sr*Ij, betarj*Sr*Ij, ...)
      }
      eq <- paste(v, collapse=" + ") # "betarj*Sr*Ij + betarj*Sr*Ij + betarj*Sr*Ij+ ..."
      ab <- paste0(d, "*", s)
      paste0(paste0("dS",r), " <- ", "-(", eq, "+", ab, ")") # dSr <- -(betarj*Sr*Ij + betarj*Sr*Ij + ...)
      dSexp <- paste0(paste0("dS",r), " <- ", "-(", eq, "+", ab, ")") # dSexp <- dSr <- -(betarj*Sr*Ij + betarj*Sr*Ij + ...)
      dIexp <- paste0(paste0("dI",r), " <- ", eq, "+", ab, " - ", paste0("g",r), "*", paste0("I",r)) # dIexp <- dIr <- betarj*Sr*Ij + betarj*Sr*Ij + ... - g * Ir
      Sv <- c(Sv, dSexp) 
      Iv <- c(Iv, dIexp)
    }
    
    ## host infection rate by abiotic reservoir
    delta <- rgamma(n=nsp, shape=sha_inf, scale=sca_inf)
  }
    
    ## sir func
    sir <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
            eval(parse(text = Sv[1])) # parse returns the parse but unevaluated form of an expression
                                      # eval- evaluates an expression in an enviornment
            eval(parse(text = Iv[1]))
            eval(parse(text = Sv[2]))
            eval(parse(text = Iv[2]))
            eval(parse(text = Sv[3]))
            eval(parse(text = Iv[3]))
            eval(parse(text = Sv[4]))
            eval(parse(text = Iv[4]))
            eval(parse(text = Sv[5]))
            eval(parse(text = Iv[5]))
            eval(parse(text = Sv[6]))
            eval(parse(text = Iv[6]))
            eval(parse(text = Sv[7]))
            eval(parse(text = Iv[7]))
            eval(parse(text = Sv[8]))
            eval(parse(text = Iv[8]))
            eval(parse(text = Sv[9]))
            eval(parse(text = Iv[9]))
            eval(parse(text = Sv[10]))
            eval(parse(text = Iv[10]))
            list(c(dS1, dI1, dS2, dI2, dS3, dI3, dS4, dI4, dS5, dI5, dS6, dI6, dS7, dI7,
                   dS8, dI8, dS9, dI9,  dS10, dI10))
        })
    }
    
    ## params for intra-family transmission gamma distribution
    ## parameter names
    pnames <- c()
    for(r in 1:nsp){   # For species r in 1:7
        for(j in 1:nsp){ # For species j in 1:7
            pnames <- c(pnames, paste0("beta", r, j))
        }
    }
    # Host mortality
    for(r in 1:nsp){   # For species r in 1:7
       pnames <- c(pnames, paste0("g", r))
    }
    
    # for abiotic reservoir
    if(abiotic == TRUE){
      for(r in 1:nsp){   # For species r in 1:7
        pnames <- c(pnames, paste0("delta", r))
      }
    }
    
    ## initial state and time vectors
    # Only 1 infected host
    state      <- c(S1 = 100, I1 = 1,
                    S2 = 100, I2 = 0,
                    S3 = 100, I3 = 0,
                    S4 = 100, I4 = 0,
                    S5 = 100, I5 = 0,
                    S6 = 100, I6 = 0,
                    S7 = 100, I7 = 0,
                    S8 = 100, I8 = 0,
                    S9 = 100, I9 = 0,
                    S10 = 100, I10 = 0
    )

  ## function to calculate time to max prevalence
  ## Average difference between pairwise comparisons
  day_range1 <- matrix(NA, nrow = nsp, ncol= nsp)
  
  Maxprevtime <- function(x, mode){
    
    if(!mode %in% c("exclude_zeros", "include_zeros")){
      stop("Mode argument must be either exclude_zeros or include_zeros")
    }
  
    if((mode == "include_zeros") & (length(x) == length(which(x == 0)))){
        max_prev <- which(x == max(x))[length(which(x == max(x)))] # # If the family never gets infected, give the last day in the series
        
    } else if((mode == "exclude_zeros") & (length(x) == length(which(x == 0)))){
        max_prev <- NA  # # If the family never gets infected, return NA
    } else {
      # If the family does get infected, give the first day in the series
      max_prev <- which(x == max(x))[1]
    }
        
    return(max_prev)
  }
  
  PairwiseDiff <- function(x){
    n <- length(x)
    dat <- matrix(NA, nrow = n, ncol= n)
    for(i in 1:n){
      for(j in 1:n){
        dat[i,j] <- x[i] - x[j]
      }
    }
    return(dat)
  }


# Loop through each simulation --------------------------------------------

    
    # Initialize number of simulations
    sim = 1
    
    while(sim <= nsim){
    
      cat("On simulation: ", sim, " ... \n")
      
    ### Single pathogen invasion with limited inter-family transmission
      #  * Each family can only infect one other family
        ## host mortality rate
        hmort <- rgamma(n=nsp, shape=sha_mort, scale=sca_mort)
        
        ## betas for transmission - initialize
        betas <- matrix(0,nrow=nsp, ncol=nsp)
        diag(betas) <- NA # Set intra to NA
        
        ##intra specific transmission
        tintra <- rgamma(n=nsp, shape=sha_intra, scale=sca_intra)
    
        ##inter specific transmission
        ninter <- ((nsp^2 - nsp)/2) # Number of values needed for lower diagonal
        tinter <- rgamma(n=ninter, shape=sha_inter, scale=sca_inter)
        
      # Initialize pairwise matrix for unique interactions  
        betas_pairwise <- betas
        betas_pairwise[lower.tri(betas_pairwise)] <- 1:ninter
        betas_pairwise <- as.matrix(Matrix::forceSymmetric(betas_pairwise,
                                                          uplo = "L")) # Make symmetric
        diag(betas_pairwise) <- NA
        
        ## Initialize dataframe for number of interactions per family
        int_df <- data.frame(family = 1:nsp,
                             interactions = 0)
        
        # Loop over columns to set interactions
        stop = FALSE # Initialize flag to break out of loop
        for(fam in 1:nsp){
        
          if(verbose) cat("\n Working on col:", fam, " ... \n")
          
          # Check to see if family we are on has reached max interactions
          if(int_df$interactions[int_df$family == fam] == interacting_families){
            cat("Reached max total number of interactions for fam:",
                fam, "... \n")
            next
          }
          
          # Get list of unique pairwise interaction labels for families not being interacted with
          possible_interactions <- na.omit(betas_pairwise[which(betas[,fam] == 0),
                                                  fam])
          
          # If no possible interactions, start over
          if(length(possible_interactions) == 0){
            cat("Breaking out of loop because no possible choices of family... \n")
            sim = sim - 1
            stop = TRUE # Trigger flag
            break
          } 
          # If only one possible choice, choose that one
          else if(length(possible_interactions) == 1){
            chosen_interactions <-  possible_interactions
          } else if(length(possible_interactions) < (interacting_families -
                    int_df$interactions[int_df$family
                                        == fam])){
            cat("Breaking out of loop because too few possible choices of family... \n")
            sim = sim - 1
            stop = TRUE # Trigger flag
            break
          } else {
          # Choose interacations to fill in
         chosen_interactions <-  sample(x = possible_interactions, 
                                        size = interacting_families -
                                          int_df$interactions[int_df$family
                                                              == fam])
          }
          
          if(verbose){
          cat("Possible pairwise interactions:", possible_interactions, 
              "... \n")
          cat("Chosen pairwise interactions:", chosen_interactions, "... \n")
          }
          
         # Set interactions
         for(interaction in chosen_interactions){
           betas[betas_pairwise %in% interaction] <- tinter[interaction]
           
           count_interaction <- which(betas_pairwise == interaction, arr.ind = T)
           count_interaction <- count_interaction[1,]
           
           # Increment number of interactions in interaction dataframe
           int_df[count_interaction, 2] <- int_df[count_interaction,2] + 1
           
         }
         
         # Remove possible interactions after the families reaches the limit
         for(fam_limit in int_df$family[int_df$interactions == interacting_families]){
           betas_pairwise[fam_limit, ] <- NA
           betas_pairwise[, fam_limit] <- NA
         }
         
         # Error check
         if(any(int_df$interactions > interacting_families)){
           stop("Number of interacting families too high")
         }
    
          if(verbose){
          print(betas)
          print(betas_pairwise)
          print(int_df)
          }
    
        } # End family / column loop  
        
        # If inner loop was stopped, then go to next simulation iteration
        if(stop == TRUE){
          sim = sim + 1
          next 
        }
    
        # Add intra-specific values
          diag(betas) <- tintra
        
        ## Check for error
        if(!(sum(betas[, 1] > 0) != interacting_families)){
          stop("Number of interacting families is not correct")
        }
          
        # Check symmetry  
        if(!isSymmetric(betas)){
          stop("Betas matrix not symmetric ... \n")
        }  
          
          

# Add in superspreader dynamics -------
 # family 1 is always the super spreader 
 # Inter-taxa transmission becomes == to intra-taxa transmission

          if(super_spreader == TRUE){
            
            # Get rid of intra specific transmission as well
            to_replace <- which(betas[, 1] > 0)[-1]
            
            # Check for errors
            if(length(to_replace) != interacting_families){
              stop("In superspreader section, number of transmission values to modify does not equal the number of interacting families... \n")
            }
            
            # Replace beta transmission values for super spreader using intraspecific transmission values
            betas[which(betas[, 1] > 0)[-1], 1] <- rgamma(n=length(to_replace), 
                                                          shape=sha_intra, 
                                                          scale=sca_intra)  
          }
          
        ## assume symetric
        bv <- as.vector(betas)
        if(abiotic == TRUE){
          pv <- c(bv, hmort, delta)
        } else {pv <- c(bv, hmort)}
        names(pv) <- pnames

       # Solve ODE
          out_t1 <- ode(y = state, times = times, func = sir, parms = pv)
        
          
        # Infected population and susceptible population sizes
        infected <- out_t1[,seq(3, by=2, length.out=nsp)]
        susceptible <- out_t1[,seq(2, by=2, length.out=nsp)]  
        
        # Calculate population size
        pop <- infected + susceptible
        
        # Calculate prevalence
        iout_t1 <- infected / pop
        

      # Calculate the proportion of individuals that died each time step
          # Number of infected individuals at t-1 minus number of individuals remaining infected at t / (total population size at t-1)
        prop_mort <- array(NA, dim = c(nrow(pop)-1, ncol(pop)))
        for(i in 2:nrow(pop)){
          prop_mort[i-1,] <- (pop[i-1,] - pop[i,]) / pop[i-1,]
       #   if(any(prop_mort[i-1,] < 0)) stop()
        }
        
        # Test plot
       # matplot(prop_mort, type = "l")
        
 
  # Test for synchrony across families - only those that got infected
        
        ## Prevalence
          # infected <- iout_t1[, colSums(iout_t1)  > 0]
          # fam_sync = kendall.w(infected, nrands = 0, type = 1)
          # res1[sim, 7] <- fam_sync$w.corrected
          # matplot(iout_t1, type = "l")
          # matplot(infected, type = "l")

        # Only infected species
          infected_only <- iout_t1[, colSums(iout_t1) > 0]
          
          sync_prev_infected <- community.sync(infected_only, nrands = 2, 
                                               type = 1, method = "pearson", 
                                                quiet = TRUE)
          
          res1[sim, 7] <- sync_prev_infected$obs
          res1[sim, 8] <- sync_prev_infected$pval
          
          ## All species
          sync_prev_all <- community.sync(iout_t1, nrands = 0, 
                                          type = 1, method = "pearson", quiet = TRUE)
          res1[sim, 9] <- sync_prev_all$obs
          
        ## Probability of mortality 
          # mort_sync = kendall.w(100 - pop, nrands = 0, type = 1)
          # res1[sim, 8] <- mort_sync$w.corrected
         # matplot(100 - pop, type = "l")
          
        # Only infected species 
          sync_mort_infected <- community.sync(prop_mort[, colSums(iout_t1) > 0], 
                                               nrands = 2,
                                               quiet = TRUE,
                                               type = 1, method = "pearson")
          res1[sim, 10] <- sync_mort_infected$obs
          res1[sim, 11] <- sync_mort_infected$pval
          
        # All species
          res1[sim, 12] <- community.sync(prop_mort, 
                                          nrands = 0,
                                          type = 1, method = "pearson")$obs
        
## Plot output for figure 1 
        if(plot == TRUE){
          
            library(reshape2)
            library(ggplot2)
          
           I_out1 <- as.data.frame(out_t1) %>%
              select(time, contains("I")) %>%
              gather("family", "abundance", -time)
          
            I_out1$family <- gsub("I", "Family ", I_out1$family)
            
            ## Subtract one day because sim starts on day 0
            I_out1$time <- I_out1$time - 1
            
            # pal <- c("palegreen", "thistle", "tan1", "yellow2", "steelblue4", "violetred2", "sienna", "pink")
            
            ggplot()+
              geom_line(data = I_out1, aes(x = time, y = abundance, col = family, lty = family),
                        lwd = 0.8,
                        show.legend = T)+
             # facet_wrap(~scenario, ncol = 1)+
             # scale_color_manual(values= pal)+
              theme_bw(8)+
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    legend.position = "top") +
              scale_x_continuous(breaks = c(0, 30, 60 , 90),
                               labels = c(0, 30, 60 , 90)) + 
              geom_vline(xintercept = 90) + 
              xlab("Time (days)")+
              ylab("Number of infected hosts")
            
            ggsave(paste0("./figures/figure 1/", interacting_families, " interacting families.pdf")
                   , height = 3, width = 4)
            
        } # End plot if
      
        
    
        # Increment sim if successful
        sim = sim + 1
     
    }


# Save output -----------------------------------------------------------


  oframe <- data.frame(interacting_families = interacting_families,
                       super_spreader = super_spreader,
                       abiotic = abiotic,
                       abiotic_trans_rate = abiotic_trans_rate,
                       # mean_days_prev = res1[,1], 
                       # mean_days_prev_nozero = res1[,2],
                       # max_days_prev = res1[,3], 
                       # max_days_prev_nozero = res1[,4],
                       # mean_days_death50 = res1[,5], 
                       # max_days_death50 = res1[,6],
                       sync_prev_infected = res1[, 7],
                       sync_prev_infected_pval = res1[, 8],
                       sync_prev_all = res1[, 9],
                       sync_mort_infected= res1[, 10],
                       sync_mort_infected_pval = res1[, 11],
                       sync_mort_all = res1[, 12])
  
  # Write csv file with data
  write.csv(oframe, paste0("./output/simResults_", sim_id, "_", Sys.Date(), ".csv"),
            row.names = FALSE)

