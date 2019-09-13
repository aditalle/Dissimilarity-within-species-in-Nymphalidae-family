#####LIBRARIES##### 
library(tidyverse)
library(data.table)
library(ape)
library(DECIPHER)
library(ggplot2)
library(Biostrings)

#####DATA RETRIVAL#####

#It takes around 5 minutes to retrieve the data from BOLD. The columns were also filtered to keep recordID, bin_uri, species name, markercode, latitude, longitude, country, and nucleotide sequence
#dfNymphalidae <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Nymphalidae&geo=all&format=tsv")[, c("recordID", "bin_uri", "species_name", "markercode", "lat", "lon" , "country" ,  "nucleotides")]

#The data frame was written to a text file
#write_tsv(dfNymphalidae, 'Nymphalidae.txt')

#The data was read from the text file made
dfNymphalidae <- read_tsv('Nymphalidae.txt')

#Making the main data frame into a data table 
dtNymphalidae <- as.data.table(dfNymphalidae)

##### FILTERING AND DATA MANIPULATION #####

#Checking if there are NAs in bin_uri 
sum(is.na(dtNymphalidae$bin_uri))

#Filtering out entries without bin_uris
dtNymphalidae <- dtNymphalidae[bin_uri %like% "[:]"]

#filtering out entries without nucleotides
dtNymphalidae <- dtNymphalidae[nucleotides %like% "[ACTG]"]

#Adding new column for checking internal gap in nucleotide sequence
dtNymphalidae[, internal_gapN := str_count(nucleotides, c("[N-]"))]

#Filtering out internal gap more than 1% 
dtNymphalidae <- dtNymphalidae[, percentage_gapN := internal_gapN/nchar(nucleotides)][!percentage_gapN > 0.01]

#Removing columns 
dtNymphalidae[, c("internal_gapN", "percentage_gapN") := NULL]

#Filtering out nucleotide sequences that are less than 640 base pair and greater than 1000 base pair. 
dtNymphalidae <- dtNymphalidae[nchar(gsub("-", "", nucleotides)) %between% c(640, 1000)]

#Filtering out entries with no species name
dtNymphalidae <- dtNymphalidae[species_name %like% "[A-Z]"]

#Checking for entries without countries 
sum(is.na(dtNymphalidae$country))

#Filtering out entries without lat or lon coordinates 
dtNymphalidae <- dtNymphalidae[lat %like% "[0123456789]"]
dtNymphalidae <- dtNymphalidae[lon %like% "[0123456789]"]

#Picking COI-5P as marker of choice 
dtNymphalidae <- dtNymphalidae[markercode %like% "COI-5P"]

#Finding unique species in data table
species <- unique(dtNymphalidae$species_name)

#Creating list for use in filtering loop and for creating a data table later
Species_Name <- c()
Dissimilarity <- c()
lon <- c()
lat <- c()
countries <- c()

#This loop is for filtering. It loops through the unique species.
for (i in species){
  #Creates a data table from the main data table for a species each iteration 
  dfcurrent <- as.data.frame(dtNymphalidae[species_name %like% i])
  dfcurrent$nucleotides <- DNAStringSet(dfcurrent$nucleotides)
  dtcurrent <- as.data.table(dfcurrent)
  #Only species with more than 10 entries will be used for analysis. This decision was made in order to ensure that a proper distance matrix can be made.
  if (nrow(dtcurrent) > 10){
  #Running alignment. diags1 was used because it allows muscle to adopt k-mer extension to firstly find the "diagonals" which is high similarity in the first iteration.
  current_alignment <- DNAStringSet(muscle::muscle(dtcurrent$nucleotides, maxiters = 2, diags1 = TRUE))
  names(current_alignment) <- dtcurrent$recordID
  dna.bin_alignment <- as.DNAbin(current_alignment)
  #A distance matrix is contructed where the values produces are the dissimilarity score between two sequences. The analysis will be based on those dissimilarity score. The model used was TN93 because it assumes the rates of transition between A-G vs C-T as well as transversions. 
  dist_matrix <- dist.dna(dna.bin_alignment, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
    #This for loop will run in the range of the current data table and append to several lists
    for (x in c(1:nrow(dtcurrent))){
      Species_Name <- c(Species_Name, i)
      #For each row in the distance matrix, the maximum dissimilarity is picked. So each sequence will have its highest dissimilarity picked. 
      Dissimilarity <- c(Dissimilarity, max(dist_matrix[,x]))
      #Appending longitude and latitude to lists
      lon <- c(lon, dtcurrent$lon[x])
      lat <- c(lat, dtcurrent$lat[x]) 
      countries <- c(countries, dtcurrent$country[x])
    }
  }
}

#Making a data table from the lists made in the filtering step. This will be the main data table for analysis 
dtdissimilarity <- data.table(Species_Name, Dissimilarity, lon, lat, countries)

#Checking the distribution of the dissimilarity scores from the new data table
hist(dtdissimilarity$Dissimilarity)

#Using a boxplot to find out the dissimilarity scores that are considered outliers which are the points of interest since we want to find where there is the most dissimilarity
boxplot(dtdissimilarity$Dissimilarity)

#Creating a new table table of where the dissimilarity is greater than 0.075
dtOutliers <- dtdissimilarity[Dissimilarity > 0.075]

#Here, I am setting Dissimilarity as the key and the order is set to descending
setkey(dtdissimilarity, Dissimilarity)
setorder(dtdissimilarity, -Dissimilarity)

#Checking that the highest values are at the top
head(dtdissimilarity)

#Picking only the highest dissimilarity for each species
dtMaxDissSpecies <- dtdissimilarity[!duplicated((dtdissimilarity$Species_Name)),]

#Picking only the highest dissimilarity for each country
dtMaxDissCountry <- dtdissimilarity[!duplicated((dtdissimilarity$countries)),]


#Using ggplot, the points from dtOutliers can be plotted. 
#A colour gradient was used to show the dissimilarity scores. A wide gradient was used in order to better distinguish the scores. The alpha option was also used to make each point more transparent for clarity. The limit of the gradient is set by the min and max of the data. 
# Also, the greater the dissimilarity is, the bigger the point gets. 
ggplot() + borders("world", colour="white", fill="lightgray") + geom_point(data = dtOutliers, aes(x = lon, y = lat, color = Dissimilarity, size = Dissimilarity), alpha = 0.15) + scale_colour_gradientn(colours=c('blue' ,'cyan','green', 'yellow','orange','red'), limits=c(min(dtOutliers$Dissimilarity),max(dtOutliers$Dissimilarity)), oob = scales::squish) + labs(title ="Significant Dissimilarity Score in species of Nymphalidae family")

#Using ggplot, the points from dtMaxDissCountry was plotted. 
#The settings for this plot uses the same gradient colours but the size does not increase.
ggplot() + borders("world", colour="white", fill="lightgray") + geom_point(data = dtMaxDissCountry, aes(x = lon, y = lat, color = Dissimilarity) , size = 1, alpha = 1) + scale_colour_gradientn(colours=c('blue' ,'cyan','green', 'yellow','orange','red'),limits=c(min(dtMaxDissCountry$Dissimilarity),max(dtMaxDissCountry$Dissimilarity)), oob = scales::squish) + labs(title ="Maximum Dissimilarity Score in Nymphalidae family per country")

#Using ggplot, the points from dtMaxDissSpecies was plotted. 
#It uses the same colour gradient as the other plots. And the point gets bigger as Dissimilarity increases. 
ggplot() + borders("world", colour="white", fill="lightgray") + geom_point(data = dtMaxDissSpecies, aes(x = lon, y = lat, color = Dissimilarity, size = Dissimilarity), alpha = 0.25) + scale_colour_gradientn(colours=c('blue' ,'cyan','green', 'yellow','orange','red'),limits=c(min(dtMaxDissSpecies$Dissimilarity),max(dtMaxDissSpecies$Dissimilarity)), oob = scales::squish) + labs(title ="Maximum Dissimilarity Score per species of Nymphalidae family")

citation(package = "ape")
citation(package = "ggplot2")
citation(package = "Biostrings")
