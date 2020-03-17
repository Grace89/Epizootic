# Script written by G. V. DiRenzo & T. S. Tunstall
# Cleaning the data for the epizootic analysis

# Load the library
library(plyr)

## Data from March, July [19-25] 2004
forrest <- read.csv('/Users/Cici/Dropbox/BdPapersTate/EcologyLetters/old/Data/1_KL Mar & Jul 04 Panama samples spreadsheet 12-4-04ttcor_GVD.csv')
# Look at data structure
str(forrest)
# Format the date
forrest$date <- paste(forrest$Year, forrest$Month, forrest$Day, sep = "-")
forrest$date <- as.Date(forrest$date, "%Y-%m-%d")
# New data frame with the information we want
fc <- data.frame(genus=forrest$Genus, 
                 species=forrest$Species, 
                 zsp=forrest$Bd, 
                 loc = forrest$Transect, 
                 park=forrest$park, 
                 month = forrest$Month, 
                 day=forrest$Day, 
                 year=forrest$Year, 
                 status=NA, 
                 date=forrest$date, 
                 family = forrest$Family)
# Look at the structure
str(fc)
# How many families are there?
length(levels(fc$family))
##############

## Data from June 2004 to 2010
rawdata <- read.csv(file="/Users/Cici/Dropbox/BdPapersTate/EcologyLetters/old/Data/pcr_GVD.csv")
# Look at structure
str(rawdata)
# Take out the data that you need that matches the 2004 data
cdata <- with(rawdata, 
              data.frame(genus=New_Genus, 
                         species=New_Species,  
                         zsp=Zoospore_Load_max, 
                         loc=Location, 
                         park = Park, 
                         month=Mo, 
                         day=Day, 
                         year=Year, 
                         status=Disease_status, 
                         date=Date, 
                         family = New_Fam))
# Look at the structure
str(cdata)
# How many families are there?
length(levels(cdata$family))
# Was the sample collected in the park?
cdata$park <- as.factor(with(cdata, ifelse(park=="Yes","y","n")))
## take out original date
cdata <- cdata[,-10]

## convert to recognizable date in R
cdata$date <- with(cdata, (paste(year, month, day, sep='-')))
cdata$date <- as.Date(cdata$date, "%Y-%m-%d")
## combine data together
cdata <- rbind(cdata, fc)
## add julian day starting from 1st sample period
cdata$julDay <- julian(cdata$date, origin=min(cdata$date, na.rm=T))

##### Fix species and genus names
## change zeteki to varius
cdata$species <- as.factor(gsub('zeteki', 'varius', cdata$species))
## panamansis
cdata$species <- as.factor(gsub('inguinalis', 'panamansis', cdata$species))
cdata$species <- as.factor(gsub('inguiNAlis', 'panamansis', cdata$species))
## talamancae
cdata$species <- as.factor(gsub('talamanca', 'talamancae', cdata$species))
cdata$species <- as.factor(gsub('talamancaee', 'talamancae', cdata$species))
## haemotiticus
cdata$species <- as.factor(gsub('heamotiticus', 'haematiticus', cdata$species))
## paradalis
cdata$species <- as.factor(gsub('paradalis', 'pardalis', cdata$species))
# Bolitaglossa
cdata$genus <- as.factor(gsub('Bolitaglossa', 'Bolitoglossa', cdata$genus))
# marina
cdata$species <- as.factor(gsub('mariNA', 'marina', cdata$species))
cdata$species <- as.factor(gsub('marina', 'marianus', cdata$species))
# espadarana
cdata$genus <- as.factor(gsub('EspadaraNA', 'Espadarana', cdata$genus))
# rana
cdata$genus <- as.factor(gsub('RaNA', 'Rana', cdata$genus))
# dendrobates
cdata$genus <- as.factor(gsub('Dendrobades','Dendrobates', cdata$genus))
# Leptodactylus
cdata$genus <- as.factor(gsub('Leptotyphlops','Leptodactylus', cdata$genus))
# phaeota
cdata$species <- as.factor(gsub('phaota','phaeota', cdata$species))

# Combine genus and species into one column
cdata$genusSpecies <- paste(cdata$genus , cdata$species, sep = "_")

# Remove any rows without a family assignment
cdata <- cdata[is.na(cdata$family) == FALSE,]

## add in data on habitat type
locList <- cdata$loc
# Repeat NA the length of cdata$loc
locType <- rep(NA, length(locList))
# Copy the levels
locs <- levels(locList)
# Identify the trails
terr <- locs[c(10,11,17,18,19,23,24,26)]
# Identify the streams
stream <- locs[c(1,2,3,4,5,6,7,9,20,21,22)]

# Locations in park, but that are not stream or trail: 
# locs[c(8, 12, 13, 14, 15, 16, 25)]

# Assign a "t", "s", or "p"
for(i in 1:length(locList))
{
    test <- locList[i]
    if(is.na(test))
        {
        }
    else
        {
            ifelse(test%in%terr, locType[i] <- 't',ifelse(test%in%stream, locType[i] <- 's', locType[i] <- 'p'))
    }
    
}
# Add habitat type column
cdata$habtype <- as.factor(locType)
# Remove anything with a p or NA
cdata <- cdata[cdata$habtype != "p",]
cdata <- cdata[is.na(cdata$habtype) == FALSE,]
cdata <- droplevels(cdata)

# Add a column with 1's
cdata$Num <- 1
#- Subset the years for the analysis
cdata45 <- subset(cdata, (year==2004 | year==2005))
# Subset to data inside the park
data0405temp <-  subset(cdata45, park=="y")
# Drop levels that aren't needed
d45all <- droplevels(data0405temp)

# Calculate the sample size for each family
sample_size <- ddply(.data = d45all, .variable = "family", .fun = summarize,
      sum = sum(Num))
# How many families are there?
nrow(ddply(.data = d45all, .variable = "family", .fun = summarize,
             sum = sum(Num)))
# How many species are there?
nrow(ddply(.data = d45all, .variable = "genusSpecies", .fun = summarize,
      sum = sum(Num)))
# Export the csv file
#write.csv(ddply(.data = d45all, .variable = c("genus", "species", "family"), .fun = summarize,
#                sum = sum(Num)), file = "/Users/Cici/Dropbox/BdPapersTate/EcologyLetters/tables/species_sample_size.csv")
#write.csv(ddply(.data = d45all, .variable = c("family"), .fun = summarize,
#                sum = sum(Num)), file = "/Users/Cici/Dropbox/BdPapersTate/EcologyLetter/tables/family_sample_size.csv")
# What is the median number of swabs collected?
median(ddply(.data = d45all, .variable = "genusSpecies", .fun = summarize,
           sum = sum(Num))[,2])
# How many species had less than 25 swabs collected?
length(which(ddply(.data = d45all, .variable = "genusSpecies", .fun = summarize,
            sum = sum(Num))[,2] > 20))
# How many swabs are there?
nrow(d45all)
# Remove families with less than 20 samples
remove <- which(d45all$family %in% sample_size[which(sample_size[,2] < 20),]$family)

d45all <- d45all[-remove,]
# How many swabs are left?
nrow(d45all)
# Calcualte sample size for each habitat type
ddply(.data = d45all, .variable = "habtype", .fun = summarize,
      sum = sum(Num))

# Calculate sample size for each family
ddply(.data = d45all, .variable = "family", .fun = summarize,
      sum = sum(Num))

# Calculate total number of observations
sum(sample_size$sum)

#Plot captures through time

time_series <- ddply(.data = d45all, .variables = c("julDay", "family"), .fun = summarize, 
      total = sum(Num))

library(ggplot2)
ggplot(data = time_series, aes(x = julDay, y = total, col = family))+ geom_line()+ geom_jitter()+
  ylab("Total number of captures")+
  xlab("Julian day")+
  theme_bw()+
  scale_color_discrete(name = "Family")+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 12, color = "black"),
        legend.text =element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#setwd("/Users/Cici/Dropbox/BdPapersTate/EcologyLetters/Figures")
#ggsave("Raw data.pdf", width = 8, height = 6)

ggplot(data = d45all, aes(x = julDay))+ 
  geom_histogram(bins = 50, col = "black")+
  ylab("Total number of captures")+
  xlab("Julian day")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 12, color = "black"),
        legend.text =element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
#setwd("/Users/Cici/Dropbox/BdPapersTate/EcologyLetters/Figures")
#ggsave("Raw data histogram.pdf", width = 8, height = 6)

#write.csv(d45all, "family_epizootic.csv")
