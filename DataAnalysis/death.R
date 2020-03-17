# Script written by G. V. DiRenzo & T. S. Tunstall
# Analyzing Mortality data

# Load libraries
library(plyr)
library(reshape2)
library(ggplot2)
library(pscl)
library(gridExtra)

## Define function to get individual cum death counts
getCounts <- function(x, timeVar, idVar, coVar = NULL, countVar){
    time <- getElement(x,timeVar) ## time or date variable
    countv <- getElement(x,countVar) ## variable the corresponds to counts
    if(!is.null(coVar)) {cv <- getElement(x, coVar)} ## extra variable
    idvar <- getElement(x, idVar)
    
    ## create vector offset by 1 to get difference in count
    offset <- c(0,countv[-length(countv)])
    newDead <- countv-offset
    
    ndf <- data.frame(new=newDead, id=idvar, time=time)
    if(exists("cv")){ ndf$cv=cv}

    fdf <- ndf[FALSE,]

##    browser()
    for(i in 1:nlevels(ndf$id)){
        cid <- levels(ndf$id)[i]
        if(exists("cv")){
            for(k in 1:nlevels(ndf$cv)){
                ccv <- levels(ndf$cv)[k]
                sndf <- with(ndf, ndf[id==cid & cv==ccv,])
                sndf$csum <- cumsum(sndf$new)
                fdf <- rbind(fdf, sndf)
            }
        }
        else{
            sndf <- with(ndf, ndf[id==cid,])
            sndf$csum <- cumsum(sndf$new)
            fdf <- rbind(fdf, sndf)
        }
    }
    return(fdf)
}
  
#################
## mortality data
#################

# Load the data
#setwd("~/Dropbox/BdPapersTate/EcologyLetters/Data/")
#deathData <- read.csv("deadNew.csv")
deathData <- read.csv("./data/deadNew.csv")

deathData <- within(deathData, date <- as.Date(paste(Year, match(deathData$Month,month.abb), Day, sep='-')))

## Fix text corrections using substitutions
deathData$Genus <- factor(gsub(" $","",deathData$Genus))
deathData$Species <- factor(gsub("^ ","",deathData$Species))
deathData$Transect <- gsub("trail","",tolower(deathData$Transect))
deathData$Transect <- factor(gsub(" $", "", deathData$Transect))
deathData$new.species <- factor(gsub("^ | $", "", deathData$new.species))

# Keep these levels for each family
fam_levs <- c("Bufonidae", "Centrolenidae", "Craugastoridae", "Dendrobatidae", "Hylidae", "Plethodontidae", "Strabomantidae")
deathData <- deathData[deathData$family %in% fam_levs,]
deathData <- droplevels(deathData)
table(deathData$family)

## Combine genus and species
deathData <- within(deathData, speciesFull <- factor(paste(Genus, Species)))
deathData <- within(deathData, newSpeciesFull <- factor(paste(new.genus, new.species)))

## add julian days
deathData$julDay <- julian(deathData$date, origin=as.Date("2004-03-05"))

## Only keep 2004 and 2005 data
death0405 <- subset(deathData, Year=="2004"|Year=="2005")

## Only keep data after first detection
death45ff <- subset(death0405, julDay>=226)

## Obtain the counts for number of amphibians
deaths <- getCounts(death45ff, timeVar="julDay", idVar="family",coVar=NULL,countVar="cum..dead")
# Rename the column
colnames(deaths)[2] <- "family" 
# Remove the first column
deaths <- deaths[,-1]
# Subtract 226 to make the values smaller
# deaths$JD <- deaths$time - min(deaths$time)

# Match the date range with other data - max is Julian day 315
deaths <- deaths[deaths$time < 316,]


# Standardize julian day
JD_mean_deaths <- mean(deaths$time)
JD_sd_deaths <- sd(deaths$time)
deaths$JD <- (deaths$time - JD_mean_deaths)/JD_sd_deaths




##### figures

ds <- ggplot(deaths, aes(time, csum))

ds + geom_line(aes(colour=family), size=1.5) + ylab("Cumulative mortality") + xlab("Julian day") + 
  scale_color_brewer(type= "qual", name="Family")

## make zero counts for days when no dead animals were found in that animals habitat
final <- data.frame(family = character(), 
                    count=integer(), 
                    julDay = integer(),
                    JD = integer()) 

# For each day when an animal was encoutered
for(day in unique(as.character(deaths$time))){
  # Subset the data
    curstate <- deaths[deaths$time==day,]
    # Look at the number of unique families
    fam <- unique(curstate$family)
    # Subset the current data
    curstate <- data.frame(family = curstate$family, 
                           count = curstate$csum, 
                           julDay = curstate$time,
                           JD = curstate$JD)
    # If the number of families = 7, then add it to the final dataset
    if(length(fam)==7){
      final <- rbind(final, curstate)
    } else{
      # If there are missing families
        # locate the missing families
        missing_family <- fam_levs[!fam_levs %in% curstate$family]
        # And create a data frame
        newstate <- data.frame(family = as.character(missing_family), 
                               count = rep(0, times = length(missing_family)), 
                               julDay = rep(curstate$julDay[1], times = length(missing_family)),
                               JD = rep(curstate$JD[1], times = length(missing_family)))
        all_dat <- rbind(newstate, curstate)
        final <- rbind(final, all_dat)
    } 
}

final$count <- as.integer(final$count)
final <- droplevels(subset(final, family!=" "))
head(final)

# Calculate presence absence
final$PA <- ifelse(final$count > 0, 1, 0)

# Family levels
final$family <- factor(final$family, levels = c("Bufonidae", "Centrolenidae","Craugastoridae",
                                                "Dendrobatidae","Hylidae" , "Plethodontidae",
                                                "Strabomantidae"))
### Presence/ absence
summary(m1 <- glm(PA ~ I(JD^2)*family-1 + JD*family, data = final, family = "binomial"))

# Create table with estimates for count model 
est_rel <- summary(m1)$coefficients[, "Estimate"]
se_rel <- summary(m1)$coefficients[, "Std. Error"]
t_rel <- summary(m1)$coefficients[, "z value"]
p_rel <- summary(m1)$coefficients[, "Pr(>|z|)"]

model_table <- data.frame(Estimate = est_rel, SE = se_rel, z = t_rel, P = p_rel)
model_table <- round(model_table, dig = 2)
#write.csv(model_table, "Mortality_estimates.csv")


newdat2 <- expand.grid(family = c("Bufonidae", "Centrolenidae","Craugastoridae",
                       "Dendrobatidae","Hylidae" , "Plethodontidae",
                       "Strabomantidae"),
                      JD = seq(min(final$JD), max(final$JD), by = 0.01))
#newdat2$julDay <- newdat2$JD + min(deaths$time)

newdat2$julDay <- (newdat2$JD * JD_sd_deaths + JD_mean_deaths) 

newdat2$pred_resp <- predict(m1, newdat2, type = "response")

#---- Plot
library(ggplot2)

pal <- c("palegreen", "thistle", "tan1", "yellow2", "steelblue4", "violetred2", "sienna")

##------ plot

bottom <-  ggplot() + 
  geom_jitter(data = final, aes(x=julDay, y = PA, col = family), height = 0.05)+
  geom_line(data=newdat2, aes(x=julDay, y= pred_resp, col = family), lwd = 1.5) +
  scale_color_manual(values= pal)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.title.x = element_text(size = 12, color = "black"),
        plot.margin=unit(c(7,0,0,9),"mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Julian day") + 
  ylab("Probability of mortality") 

# Testing for synchrony in probability of mortality -------------------------------------

library(tidyverse)
library(synchrony)

# Take only columns we need  
fam_mort = newdat2 %>%
  dplyr::select(JD, family, pred_resp) %>%
  tidyr::spread(key = family, value = pred_resp)

# Make sure lines make sense
matplot(fam_mort[, -1], type = "l")

# Test for synchrony across fam_mortilies  

fam_mort_sync = community.sync(fam_mort[, -1], nrands = 10000, type = 2, 
                                   alternative = "greater")
fam_mort_sync

mean(fam_mort_sync$rands)

community.sync(fam_mort[, -1], nrands = 1000, type = 1, alternative = "less")
community.sync(fam_mort[, -1], nrands = 1000, type = 1, alternative = "greater")
community.sync(fam_mort[, -1], nrands = 1000, type = 2, alternative = "less")
community.sync(fam_mort[, -1], nrands = 1000, type = 2, alternative = "greater")


plot(fam_mort_sync)

# End script
