# Script written by G. V. DiRenzo & L. M. Browne
# Analyzing Bd epizootic prevalence and infection intensity

# Load libraries
library(pscl)
library(tidyverse)
library(cowplot)
library(plotrix)

# Read in the data
setwd("/Users/Cici/Dropbox/BdPapersTate/EcologyLetters/Data/")
dat <- read.csv("./family_epizootic.csv")[,-1]

## For lukes computer
dat <- read.csv("./Data/family_epizootic.csv")[,-1]

# Look at structure of the data
str(dat)

# Calculate number of unique surveys
length(unique(dat$date))
# First survey
levels(dat$date)[1]
# Last survey
levels(dat$date)[length(levels(dat$date))]

# Number of surveys per month per transect
library(plyr)
ddply(.data = dat, .variable = c("month", "year", "loc"), .fun = summarize,
      total = length(unique(date)))
ddply(.data = dat, .variable = c("genusSpecies"), .fun = summarize,
      total = length(unique(date)))
# Total number of observations in the data set
sum(dat$Num)

# Total number of families
length(unique(levels(dat$family)))

# Total number of species
length(unique(levels(dat$genusSpecies)))

# Total sample sizes per family
sample_size <- ddply(.data = dat, .variable = "family", .fun = summarize,
      sum = sum(Num))
# write.csv(sample_size, "sample_size.csv")

# Round the zoospore values
dat$R_zoo <- as.integer(dat$zsp)

# Total number of Bd positives by habitat type
dat$Bd_PA <- ifelse(dat$R_zoo > 0, 1, 0)
hab <- table(dat$habtype, dat$Bd_PA)
#write.csv(hab, "hab_sample_size.csv")

# Total number of captures pre- and post-Bd
dat$date <- as.Date(dat$date)
dat$Pre_post <- ifelse(dat$date > "2004-10-16", 1, 0)
  # 1 = Post (after october)
  # 0 = Pre (before october)
Pre_post <- ddply(.data = dat, .variable = c("family", "Pre_post"), .fun = summarize, sum = sum(Num))
#write.csv(Pre_post, "Pre_post.csv")

# Calculate the date of first Bd detection
Bd_positives <- dat[dat$R_zoo > 0,]
ddply(.data = Bd_positives, .variable = "family", .fun = summarize,
      date = min(date),
      JD = min(julDay))

# Julian day of first Bd positive detection
dat[dat$date == "2004-10-16",]$julDay[1]

## Keep the data past the first positive
date_pos <- min(dat[which(dat$R_zoo > 0),]$julDay)

dat <- dat[dat$julDay > date_pos | dat$julDay == date_pos,]

#---- Family model

# Standardize julian day
JD_mean <- mean(dat$julDay)
JD_sd <- sd(dat$julDay)
dat$JD <- (dat$julDay - JD_mean)/JD_sd

 
 summary(m2 <- hurdle(R_zoo ~ I(JD^2)*family-1 + JD*family | I(JD^2)*family-1 + JD*family
                      , 
                      data = dat,
                      dist = "poisson",
                      zero.dist = "binomial"
 ))
 

 # Create predictions
 newdat <- expand.grid(family =  c("Bufonidae", "Centrolenidae","Craugastoridae",
                                   "Dendrobatidae","Hylidae" , "Plethodontidae",
                                   "Strabomantidae"),
                       JD = seq(min(dat$JD), max(dat$JD), by = 0.01))
 newdat$julDay <- (newdat$JD * JD_sd) + JD_mean
 
 newdat$pred_count <- predict(m2, newdat, type = "count")
 newdat$pred_zero <- predict(m2, newdat, type = "zero")

 # Create table with estimates for count model 
 est_rel <- summary(m2)$coefficients$count[, "Estimate"]
 se_rel <- summary(m2)$coefficients$count[, "Std. Error"]
 t_rel <- summary(m2)$coefficients$count[, "z value"]
 p_rel <- summary(m2)$coefficients$count[, "Pr(>|z|)"]
 
 model_table <- data.frame(Estimate = est_rel, SE = se_rel, t = t_rel, P = p_rel)
 model_table <- round(model_table, dig = 2)
 #write.csv(model_table, "Count_estimates.csv")

 # Create table with estimates for zero model 
 est_rel <- summary(m2)$coefficients$zero[, "Estimate"]
 se_rel <- summary(m2)$coefficients$zero[, "Std. Error"]
 t_rel <- summary(m2)$coefficients$zero[, "z value"]
 p_rel <- summary(m2)$coefficients$zero[, "Pr(>|z|)"]
 
 model_table <- data.frame(Estimate = est_rel, SE = se_rel, t = t_rel, P = p_rel)
 model_table <- round(model_table, dig = 2)
 #write.csv(model_table, "Zero_estimates.csv")
 
 #--- calculate Max prevalence
 nfam <- 7 # number of families
 peak <- matrix(NA, nrow = 7, ncol = 6) # Matrix to store values
 rownames(peak) <- levels(newdat$family)
 colnames(peak) <- c("Max Inf Int", "J Day", "Max prev", "J Day", "Max death", "J Day")
 
 for(i in levels(newdat$family)){
   # Subset the data
   sub <- newdat[newdat$family == i,]
   
   # Max infection intensity
   row <- which(sub$pred_count == max(sub$pred_count))
   
   peak[i,1] <- sub[row,]$pred_count # Infection 
   peak[i,2] <- sub[row,]$julDay # Julian day
   
   # Max prevalence
   row <- which(sub$pred_zero == max(sub$pred_zero))
   
   peak[i,3] <- sub[row,]$pred_zero # Prevalence
   peak[i,4] <- sub[row,]$julDay  # Julian day
   
   # Max mortality
   sub <- newdat2[newdat2$family == i,]
   row <- which(sub$pred_resp == max(sub$pred_resp))
   
   peak[i,5] <- sub[row,]$pred_resp # Prob of mortality
   peak[i,6] <- sub[row,]$julDay  # Julian day
 }
 
 peak <- round(peak, dig = 2)
 
 #write.csv(peak, "peak.csv")
 
 
 #----- To make the multi-panel plot: Run the death.R file first
 #---- Plot

 
 text_size <- 15
 label = "*"
 label_size = 13
 
 
 pal <- c("palegreen", "thistle", "tan1", "yellow2", "steelblue4", "violetred2", "sienna")
 lwd = 1.25
 ##------ plot
 top <- ggplot() + 
#   geom_jitter(data = dat, aes(x=julDay, y = log10(R_zoo+1), col = family), width = .5, height = 0)+
   scale_y_continuous(breaks = c(log10(1), log10(10), log10(100), log10(1000), log10(10000),
                                 log10(100000), log10(1000000)), 
                      labels = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                      limits = c(0, log10(500000)))+
   scale_color_manual(values= pal)+
   geom_line(data=newdat, aes(x=julDay, y=log10(pred_count), col = family, lty = family),
             lwd = lwd) +
   theme_bw()+
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(size = text_size, color = "black"),
         axis.title.y = element_text(size = text_size, color = "black"), 
         axis.title.x = element_text(size = text_size, color = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "top",
         legend.text=element_text(size=text_size-3.5),
         legend.title=element_text(size=text_size)) + 
   guides(col=guide_legend(title="Family"), lty = guide_legend(title = "Family"))+
   xlab("") + 
   ylab(expression(paste(italic(Bd), " zoospores + 1")))+
   annotate("text", x = peak[1,2], y = log10(peak[1,1]+1), label =label, size = label_size, col = pal[1])+
   annotate("text", x = peak[2,2], y = log10(peak[2,1]+1), label =label, size = label_size, col = pal[2])+
   annotate("text", x = peak[3,2], y = log10(peak[3,1]+1), label =label, size = label_size, col = pal[3])+
   annotate("text", x = peak[4,2], y = log10(peak[4,1]+1), label =label, size = label_size, col = pal[4])+
   annotate("text", x = peak[5,2], y = log10(peak[5,1]+1), label =label, size = label_size, col = pal[5])+
   annotate("text", x = peak[6,2], y = log10(peak[6,1]+1), label =label, size = label_size, col = pal[6])+
   annotate("text", x = peak[7,2], y = log10(peak[7,1]+1), label =label, size = label_size, col = pal[7])

 
 middle <-  ggplot() + 
 #  geom_jitter(data = dat, aes(x=julDay, y = Bd_PA, col = family), height = 0.05,)+
   geom_line(data= newdat, aes(x=julDay, y=(pred_zero), col = family, lty = family), lwd = lwd) +
   scale_color_manual(values= pal)+
   theme_bw()+
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(size = text_size, color = "black"),
         axis.title.y = element_text(size = text_size, color = "black"), 
         axis.title.x = element_text(size = text_size, color = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none") +
   xlab("") +
   ylim(c(0, 1.05)) +
   ylab(expression(paste(italic(Bd), " prevalence"))) +
   annotate("text", x = peak[1,4], y = peak[1,3], label =label, size = label_size, col = pal[1])+
   annotate("text", x = peak[2,4], y = peak[2,3], label =label, size = label_size, col = pal[2])+
   annotate("text", x = peak[3,4], y = peak[3,3], label =label, size = label_size, col = pal[3])+
   annotate("text", x = peak[4,4], y = peak[4,3], label =label, size = label_size, col = pal[4])+
   annotate("text", x = peak[5,4], y = peak[5,3], label =label, size = label_size, col = pal[5])+
   annotate("text", x = peak[6,4], y = peak[6,3], label =label, size = label_size, col = pal[6])+
   annotate("text", x = peak[7,4], y = peak[7,3], label =label, size = label_size, col = pal[7])
 
 bottom <-  ggplot() + 
#   geom_jitter(data = final, aes(x=julDay, y = PA, col = family), height = 0.05)+
   geom_line(data=newdat2, aes(x=julDay, y= pred_resp, col = family, lty = family), lwd = lwd) +
   scale_color_manual(values= pal)+
   theme_bw()+
   theme(axis.text.x = element_text(size = text_size, color = "black"),
         axis.text.y = element_text(size = text_size, color = "black"),
         axis.title.y = element_text(size = text_size, color = "black"), 
         axis.title.x = element_text(size = text_size, color = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none") +
   xlab("Julian day") + 
   ylab("Probability of mortality") +
   ylim(c(0, 1.05)) +
   annotate("text", x = peak[1,6], y = peak[1,5], label =label, size = label_size, col = pal[1])+
   annotate("text", x = peak[2,6], y = peak[2,5], label =label, size = label_size, col = pal[2])+
   annotate("text", x = peak[3,6], y = peak[3,5], label =label, size = label_size, col = pal[3])+
   annotate("text", x = peak[4,6], y = peak[4,5], label =label, size = label_size, col = pal[4])+
   annotate("text", x = peak[5,6], y = peak[5,5], label =label, size = label_size, col = pal[5])+
   annotate("text", x = peak[6,6], y = peak[6,5], label =label, size = label_size, col = pal[6])+
   annotate("text", x = peak[7,6], y = peak[7,5], label =label, size = label_size, col = pal[7])
 

 plot_grid(top,middle,bottom,nrow=3, labels = c("(a)", "(b)", "(c)"), label_size = 16,
           vjust = c(.95, 0, 0))
 
  ggsave(paste0("./figures/Family_epizootic ", Sys.Date(), ".pdf"), height = 9, width = 8)
 
  
  
  
  
  
 # Testing for synchrony in prevalence -------------------------------------
 
 library(tidyverse)
 library(synchrony)
 
# Take only columns we need  
 fam = newdat %>%
   dplyr::select(JD, family, pred_zero) %>%
   tidyr::spread(key = family, value = pred_zero)
 
 # Make sure lines make sense
 matplot(fam[, -1], type = "l")

 # Test for synchrony across families in disease prevalence

 fam_sync = community.sync(fam[, -1], nrands = 10000, type = 2, alternative = "greater")
 fam_sync
 
 mean(fam_sync$rands)
 
 community.sync(fam[, -1], nrands = 1000, type = 1, alternative = "less")
 community.sync(fam[, -1], nrands = 1000, type = 1, alternative = "greater")
 community.sync(fam[, -1], nrands = 1000, type = 2, alternative = "less")
 community.sync(fam[, -1], nrands = 10000, type = 2, alternative = "greater")
 
 
 
 # Testing for synchrony in infection intensity -------------------------------------
 
 # Take only columns we need  
 famC = newdat %>%
   dplyr::select(JD, family, pred_count) %>%
   tidyr::spread(key = family, value = pred_count)
 
 # Make sure lines make sense
 matplot(famC[, -1], type = "l")
 
 # Test for synchrony across families  
 
 famC_sync = community.sync(log10(famC[, -1] + 1), 
                            nrands = 10000, type = 2, 
                            alternative = "less")
 famC_sync
 
 mean(famC_sync$rands)
 
 community.sync(log10(famC[, -1] + 1), nrands = 1000, type = 1, alternative = "less")
 community.sync(log10(famC[, -1] + 1), nrands = 1000, type = 1, alternative = "greater")
 community.sync(log10(famC[, -1] + 1), nrands = 1000, type = 2, alternative = "less")
 community.sync(log10(famC[, -1] + 1), nrands = 1000, type = 2, alternative = "greater")
 
 # To make full plot, need to run death.R first

 
 # Custom function for plotting
 plot_sync <- function (x, main = "",
                        xlab = "Community-wide synchrony metric", 
                        ylab = "Frequency", line.col = "#377eb8",
                        lty = 1, lwd = 3, col = "grey90", 
                        ...) 
 {
   if (!is.null(x$rands) & class(x) == "synchrony") {
     if (!is.null(x$w.corrected)) 
       x$obs = x$w.corrected
     hist(x$rands, main = main, xlab = xlab, ylab = ylab, xlim = c(0, 1),
          col = col, breaks = 50, ...)
     abline(v = x$obs, col = line.col, lty = lty, lwd = lwd, 
            ...)
    # box()
   }
   else {
     stop("No permutation data to plot")
   }
 }
 
 
 ## Plot for supplementary information
 
 par(mfrow = c(1, 3))

 plot_sync(fam_sync, las = 0, main = "Prevalence", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
 mtext(side = 3, text = "(a)", adj = 0, padj = -1.05, cex = 1.5)
 plot_sync(famC_sync, las = 0, main = "Infection intensity", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
 mtext(side = 3, text = "(b)", adj = 0, padj = -1.05, cex = 1.5)
 plot_sync(fam_mort_sync, las = 0, main = "Probability of mortality", cex.axis = 1.5, cex.lab = 1.5)
 mtext(side = 3, text = "(c)", adj = 0, padj = -1.05, cex = 1.5)
 
 # Save plot
 dev.copy(pdf, paste0("./figures/Synchrony plots ", Sys.Date(), ".pdf"), height = 4, width = 15)
 dev.off()
 
 
 # End script
 
 
 ## Testing
 
 peak2 <- as.data.frame(peak)
 
 peak2$status <- c("c","e", "e", "e", "e","c","e")
 peak2$sample_size <- c(183, 215, 162,96, 58, 55,259)
 
 boxplot(peak2[, 2] ~ peak2$status, main = "Infection intensity")
 boxplot(peak2[, 4] ~ peak2$status, main = "Prevalence")
 boxplot(peak2[, 6] ~ peak2$status, main = "Mortality")
 
 boxplot(peak2[, 1] ~ peak2$status, main = "Infection intensity")
 boxplot(peak2[, 3] ~ peak2$status, main = "Prevalence")
 boxplot(peak2[, 5] ~ peak2$status, main = "Mortality")
 
 plot(peak2$sample_size, peak2[, 2], main = "Infection intensity", pch = 19)
 plot(peak2$sample_size, peak2[, 4], main = "Prevalence", pch = 19)
 plot(peak2$sample_size, peak2[, 6], main = "Mortality", pch = 19)
 
 plot(peak2$sample_size, peak2[, 1], main = "Infection intensity", pch = 19)
 plot(peak2$sample_size, peak2[, 3], main = "Prevalence", pch = 19)
 plot(peak2$sample_size, peak2[, 5], main = "Mortality", pch = 19)
 
 
 
 
 