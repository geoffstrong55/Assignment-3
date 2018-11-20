###Assignment 1C ----

#Tardigrada are water-dwelling micro-animals that live in a variety of different environments. Tardigrades have been found in areas such as mountain tops, volcanoes, the deep sea, tropical rainforests, and even the Antarctic. 

#Since these species live in a variety of different environments, I am interested in seeing if there are any differences in the samples' DNA structure from different environments. The first analysis I want to perform is to calculate AT proportion in species from different elevations to see if there is a significant structure in the DNA of the different samples. 

#Once I have completed that analysis, I want to compare the species richness between different geographical areas.

###Loading required packages ----
#I will load the packages that I will be using for my data analysis

library(tidyverse)
library(vegan)
library(ape)
library(stringi)
library(RSQLite)
library(iNEXT)
library(devtools)
library(Biostrings)
library(ggplot2)

###Loading data and data exploration ----
#The species that I will be observing is Tardigrada. I will load the data from the BOLD API directly

Tardigrada <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Tardigrada&format=tsv")

#I want to check the class of Tardigrada to see waht structure I am dealing with.

class(Tardigrada)

#Since Tardigrada is a data frame I want to change it to tibble. I want to do this because I want to make data wrangling later in my analysis, a little easier

Tardigrada <-as.tibble(Tardigrada)

#The first thing I want to analyze from the dataset is the number of unique species. I will use the bin_uri as a species identifier.

length(unique(Tardigrada$bin_uri))

#There is a total of 162 different species in the data set.


### Data filtering and subsetting ----
#As mentioned above, I want to analyze the nucleotide sequences of samples in different elevations to see if there is a difference in DNA structure. To do this I will be looking at the AT proportion of each samples DNA. I want to subset the data so I can look at the samples that have elevation values and their nucleotides sequenced. For this new data set, I will continue to use the BIN as the identifier for each sample. I will name this new dataset Tardigrada2. 

Tardigrada2 <- Tardigrada %>%
  select(bin_uri, elev, nucleotides) %>%
  filter(!is.na(elev)) %>%
  filter(!is.na(bin_uri)) %>%
  filter(!is.na(nucleotides))

###Analysis based on Elevation ----

#After this round of subesetting the data, there are now 186 observations after filtering for samples containing both elevation values and a nucleotide sequence. Next, I want to create a histogram to look at how the elevation values are distributed.  

hist(Tardigrada2$elev, col = 2, main = "Frequency of Species by Elevation",xlab = "Elevation", ylab = "Frequency")

#Based on the histogram, there is a large proportion of samples that are close to sea level. When looking at Tardigrada2, I noticed many samples are in elevations under 250. I want to know the exact number of species at different elevations. I will separate elevation values into three different groups. The three groups will be 0-250, 250-1000, and 1000-2100. Next, I will filter Tardigrada2 into 3 new tibbles named Tardi.0_250, Tardi.250_1000, and Tardi.1000_2100, respectively. I will use nrow to calculate the number of samples in each group.

Tardi.0_250 <- Tardigrada2 %>%
  filter(elev <= 250)

nrow(Tardi.0_250)

Tardi.250_1000 <- Tardigrada2 %>%
  filter(elev > 251, elev <= 1000)

nrow(Tardi.250_1000)

Tardi.1000_2100 <- Tardigrada2 %>%
  filter(elev > 1000, elev <= 2100)

nrow(Tardi.1000_2100)

#There are 108 samples in Tardi.0_250, 40 in Tardi.250_1000, and 38 in Tardi.1000_2100.

#Now that the species samples are separated into 3 groups based on elevation, I want to look at the AT proportion of their nucleotides. I want to find this value because I want to see if there is a difference in the DNA of the samples at different elevations. To look at the AT proportions, I will have to change nucleotides to biostrings character using the function DNAStringSet. 

Tardi.0_250$nucleotides <- DNAStringSet(Tardi.0_250$nucleotides)

Tardi.250_1000$nucleotides <- DNAStringSet(Tardi.250_1000$nucleotides)

Tardi.1000_2100$nucleotides <- DNAStringSet(Tardi.1000_2100$nucleotides)

#I will confirm the class of the nucleotides column.

class(Tardi.0_250$nucleotides)
class(Tardi.250_1000$nucleotides)
class(Tardi.1000_2100$nucleotides)

#Since the nucleotides columns are now "Biostrings", I will calculate the AT proportion for each sample in the 3 groups. I will use the letter frequency function to calculate the frequency of each nucleotide within each sample's DNA. I will create a new tibble for each of the three groups. After creating a new tibble, I will calculate the AT proportion of each BIN. I will also create a histrogram for each group to visualize the frequency.

Tardi.0_250.Freq <- as_tibble(letterFrequency(Tardi.0_250$nucleotides, letters = c("A", "C", "G", "T"))) %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C))) 

Tardi.250_1000.Freq <- as_tibble(letterFrequency(Tardi.250_1000$nucleotides, letters = c("A", "C", "G", "T"))) %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C)))

Tardi.1000_2100.Freq <- as.tibble(letterFrequency(Tardi.1000_2100$nucleotides, letters = c("A", "C", "G", "T"))) %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C)))

#We can clean up the global environment by removing some objects that are no longer needed. 
rm(Tardi.0_250, Tardi.1000_2100, Tardi.250_1000)

#To see a comparison of histograms generated above based on the elevation groups determined, we can use par() to stack them together. This allows for an easier visualization of the frequency of AT proportions based on the subsetted elevations. You can generate each histogram separately or you can also highlight the entire portion of code below to stack all 3 histograms vertically to compare. 

par(mfrow=c(3,1))
hist(Tardi.0_250.Freq$ATproportion, col = 3, main = "AT Proportion for group Tardi.0_250",xlab = "AT Proportion", ylab = "Frequency", xlim = c(0.58, 0.72), ylim = c(0, 50))
hist(Tardi.250_1000.Freq$ATproportion, col = 4, main = "AT Proportion for group Tardi.250_1000",xlab = "AT Proportion", ylab = "Frequency", xlim = c(0.58, 0.72), ylim = c(0, 50))
hist(Tardi.1000_2100.Freq$ATproportion, col = 5, main = "AT Proportion for group Tardi.1000_2100",xlab = "AT Proportion", ylab = "Frequency", xlim = c(0.58, 0.72), ylim = c(0, 50))

#Now that the AT proportion has been calculated for each group, I want to compare each group to one another to see if there is a staistical significant difference between the groups. This will be done by running a t-test comparing each group.

t.test(Tardi.0_250.Freq$ATproportion, Tardi.250_1000.Freq$ATproportion)

t.test(Tardi.0_250.Freq$ATproportion, Tardi.1000_2100.Freq$ATproportion)

t.test(Tardi.250_1000.Freq$ATproportion, Tardi.1000_2100.Freq$ATproportion)

#When comparing Tardi.0_250.Freq to Tardi.250_1000.Freq the p-value is 0.1271. The p-value when comparing Tardi.0_250.Freq to Tardi.1000_2100.Freq is 0.236. Lastly, when comparing Tardi.250_1000.Freq to Tardi.1000_2100.Freq the p-value is 0.7677. Based on the p-value obtained from each group we can say that DNA structure based on AT proportion is similar when comparing the three groups to one another. Future research looking at how certain genes are expressed and the types of proteins produced by species in different elevations would be interestng to see if there are any differences that allow each species to thrive in their environment.

#This was a very thorough analysis based on separated elevations. Instead of performing 3 different t-tests, perhaps we can look at using the ANOVA statistical test to analyze the variance between elevation and AT proportion. 

#First, to analyze the AT proportions of the Tardigrada samples, we need to convert Tardigrada2 into a DNA stringset. 
Tardigrada2$nucleotides <- DNAStringSet(Tardigrada2$nucleotides)

#Next, we can calculate the AT proportions of each sample and use the mutate function to add a column to the new data frame Tardigrada2.Freq. This data frame shows the frequency of each nucleotide plus the AT proportions for the samples. 
Tardigrada2.Freq <- as_tibble(letterFrequency(Tardigrada2$nucleotides, letters = c("A", "C", "G", "T"))) %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C))) 

#To perform the ANOVA analysis between elevation and AT proportions, we need to combine the two data frames so we have elevation data and AT proportions in a single data frame. This can be done using the cbind function. 
Tardigrada.elev <- cbind(Tardigrada2, Tardigrada2.Freq)

#Finally, we can perform the analysis of variance using the aov() from the stats package and assign it to a varible "Tardi_aov". From there, you can print a summary of the aov using summary.aov() to show the F statistic and p-value for statistical significance. 
Tardi_aov <- aov(Tardigrada.elev$elev~Tardigrada.elev$ATproportion)
summary.aov(Tardi_aov)

#For the next analysis, we won't be needing some of the objects in the global environment anymore. 
rm(Tardi.0_250.Freq, Tardi.250_1000.Freq, Tardi.1000_2100.Freq, Tardi_aov, Tardigrada.elev, Tardigrada2.Freq)

###Analysis of Species Richness----
#The next thing I want to analyze is the species abundance of this data set. To do this I will create a new tibble that contains the country that each sample is taken from, the BIN, and the count of samples for each BIN, which will be represented by n. I also wanto know the country that has the largest number of samples. 

Tardigrada %>%
  group_by(country) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))

#The countries with the highest number of samples are Antarctica, Spain, Italy, the United Kingdom, and Chile. This countries will be used for a later analysis. I will then create a new tibble with the country names, bin_uri, and the value n representing the number of samples collected.

Tardigrada3 <- Tardigrada %>%
  group_by(country, bin_uri) %>%
  count(bin_uri)

#Now that I have the species counted by BIN, I need to put the data into a community object so it can be used for analysis to calculated the Fisher's Alpha. I will make a new tibble called Tardigrada4. I will make a community object using the function spread().

Tardigrada4 <- spread(Tardigrada3, bin_uri, n)

#To be able to run this new data frame in a function in vegan, I will have to replace the NAs in the data frame to a numeric value of 0. 

Tardigrada4[is.na(Tardigrada4)] <- 0

#Next, I will create a new data set Tardigrada5, removing the row names, and switching the names of the columns to country. 

Tardigrada5 <- Tardigrada4 %>%
  remove_rownames %>%
  column_to_rownames(var="country")

#To analyze the species abundance of Tardigrada, I am going to calculate the Fisher's alpha, to look at the relative adundance in each country and plot the Fisher's alpha of each country.

fisher.alpha(Tardigrada5)
plot(fisher.alpha(Tardigrada5))

#I am going to use a Fisher's fit model. This will calculate the expected abudnace of each species in each country. 

plot(fisherfit(Tardigrada5))

#The plot produced shows thats there is a large number of speices expected to be only sampled 1 time, then a small number of speices will be sample with larger frequency. A higher number closer the the value one shows larger species diversity compared to a plot that would have a larger number of species with higher sample frequency. 

#The next analysis I want ot perform, is to look at the species richness of the phylum Tardigardia. I will do this by using the iNEXT package. I will look only at the 5 countries with the highest number of samples. This number of 5 was chosen for simplicity of analysis and the visualization of the data. As stated before the 5 countries with the highest number of samples are Antarctica, Spain, Italy, the United Kingdom, and Chile. 

#I will have the variable country be the columns and I will have to remove the NAs from the data set.

Tardigrada6 <- spread(Tardigrada3, country, n)
Tardigrada6[is.na(Tardigrada6)] <- 0

#As stated previously, I am only interested in the 5 countries with the most number of samples. To run the iNEXT function also values must be numeric. I will change the first column to row names by using the column_to_rownames function

Tardigrada7 <- Tardigrada6 %>%
  select(Antarctica, Italy, Spain, Chile, `United Kingdom`) %>%
  column_to_rownames(var="bin_uri")

#Next, I need to confirm that the object Tardigrada7 is a tibble. The reason I need to confirm this is to due iNEXT not taking tibbles as an argument.

class(Tardigrada7)

#Since Tardigrada7 is a tibble, I need to change it to a matrix.

Tardigrada8 <- as.matrix(Tardigrada7)

#lastly, I will use the letter m to create a vector for the number of species needed for the iNEXT function. I will create a new object called Tardi.Richness with the species richness values for Tardigrada8. I will then use ggiNEXT to create a plot showing the expected species richness values for each of the 5 countries selected. 

m <- c(1, 5, 20, 50, 100, 200, 400)

Tardi.Richness <- iNEXT(Tardigrada8, q = 0, "abundance", size = m)    

ggiNEXT(Tardi.Richness, type=1)

#The plot created shows that Antarctica has the highest species richness and estimated species richness. Italy is second for both of those value. Chile, United Kingdom, and Spain have the lowest species richness values, respectively. This analysis shows that Antarctica has the highest number of sample species richenss and the highest diversity of the phylum Tardigrada.
