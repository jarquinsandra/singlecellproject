
## load Sebastian Gibb's MALDIquant and MALDIquantForeign libraries
## Gibb S, Strimmer K. MALDIquant: a versatile R package for the analysis of mass spectrometry data. 
## Bioinformatics. 2012 Sep 1;28(17):2270-1. doi:10.1093/bioinformatics/bts447. 
## Epub 2012 Jul 12. PubMed PMID: 22796955.
#Script modified and complemented to analyze differences in features between two setup conditions for an ambient ionization source

library("MALDIquant")
library("MALDIquantForeign")
library("vegetarian")
library("pheatmap")
#library("corrplot")
#library("PerformanceAnalytics")
#library("FactoMineR")
#library("factoextra")
library(tidyverse)
#library("ggpubr")
#library(RVAideMemoire)
library(RColorBrewer)
library(viridis)
library(Rtsne)
#library(randomForest)
#require(caTools)
#library(caret)


## get names of all .mzML files in current working directory
#windows
path <- "~/Sandra_MS_data/technicalnote/positive"
#mac
path <-"/Users/sandramartinez/Downloads/technical_note_data/positive"

#import spectrum
spectra<-importMzMl(path = path)
#normalization to TIC
spectra <- calibrateIntensity(spectra, method="TIC")
## Peak picking, with signal to noise filtering
peaks<-detectPeaks(spectra, SNR=3)
## Peak alignment using the warping function of Maldiquant
warpingFunctions<-determineWarpingFunctions(peaks,allowNoMatches=TRUE)
peaks<-warpMassPeaks(peaks,warpingFunctions)
#Peak binning and creation of feature matrix
peaks<-binPeaks(peaks)
featureMatrix<-intensityMatrix(peaks)
## Replace NAN values by 0
featureMatrix[is.na(featureMatrix)] <- 0
#remove all features with more than 20% values by class
features<-as_tibble(featureMatrix,.name_repair = ~ make.names(., unique = TRUE))
#Remove columns that have more than 20% 0 values
zerocutoff<-nrow(features)*0.2
features[colSums(features > 0) <= zerocutoff]<- NULL

 #This is filled with the number of spectrum by class for either negative or positive mode, the numbers here correspond to the ones from the files used in the paper 
#Negative spectra
ASPC<-replicate(1268, "ASPC-1")
HeLa<-replicate(1264, "HeLa")
Panc10<-replicate(1267, "Panc 10.05")
alllabels<-c(ASPC,HeLa,Panc10)
 #positive spectrum
ASPC<-replicate(1265, "ASPC-1")
HeLa<-replicate(1262, "HeLa")
Panc10<-replicate(1263, "Panc 10.05")
 
 #add labels to the spectra, they are read in alphabetic order
alllabels<-c(ASPC,HeLa,Panc10)
features$sample<-alllabels
#add id number as unique identifier per cell, this is needed for some loops
features<-tibble::rowid_to_column(features, "id")
#moving important columns to the first two places
features <- features %>%
select(id,sample, everything())
#Round features name to 4 digits
names<-colnames(features)
names<-sub('X', '', names)
names<-names[3:ncol(features)]
names<-as.numeric(names)
names<-round(names,digits=4)
names<-as.character(names)
#add X to characters
names<-paste0("X", names)
labels<-c("id","sample",names)
colnames(features)<-labels
#Filter cells positive 
cells<-features[!(features$X732.5554==0|features$X760.5894)==0,]
#Filter cells Negative 
cells<-features[!(features$X804.5752==0|features$X830.5872)==0,]


# Curating the database for analysis with both t-SNE and PCA
Labels<-cells$sample
#convert sample type to factor 
cells$sample<-as.factor(cells$sample)
## for plotting
colors = rainbow(length(unique(cells$sample)))
names(colors) = unique(cells$sample)
##modified from https://datavizpyr.com/how-to-make-tsne-plot-in-r/
## Executing the algorithm on curated data
tsne <- Rtsne(cells[,3:ncol(features)], dims = 2, perplexity=50, verbose=TRUE, max_iter = 5000)
#get labels
spectranumcol<-cells$id
samplecol<-cells$sample
tSNE_df <- tsne$Y %>% 
  as.data.frame() %>%
  dplyr::rename(tsne1="V1",
                tsne2="V2") %>%
  mutate(sample=samplecol, spectranum=spectranumcol)


tSNE_df %>%
  ggplot(aes(x = tsne1, 
             y = tsne2,
             color = sample))+
  geom_point(size=0.5)+
  theme_classic() +
  theme(legend.position = "top",text = element_text(size = 20),panel.border = element_rect(colour = "black", fill=NA, size=1))+ guides(colour = guide_legend(override.aes = list(size=5)))
#change file name or it will overwrite the existing file
ggsave("technicalnotepositive.png")

##to create the heatmap sampling data was done
#get names from the mz values
set.seed(123)
rfdata<-cells[,2:ncol(cells)]
rfdata$sample<-as_factor(rfdata$sample)
ind <- sample(2, nrow(rfdata), replace = TRUE, prob = c(0.95, 0.05))
train <- rfdata[ind==1,]
test <- rfdata[ind==2,]
#get the dimensions on the subsample set
dim(test)

#####The heatmaps can be customised, for positive mode we selected 150 features based on their intensities, for negative mode the total number of features is 104 so it was not necessary to select features.
###For positive data heatmap 
columnmax<-apply(test[2:ncol(test)], 2, max)
#sort and select the 150 most intense values
colmax<-sort(columnmax,decreasing=TRUE)
colmaxselection<-as.matrix(colmax[0:150])
colmaxselection=as.vector(rownames(colmaxselection))
#limiting matrix to selection
max_normout_featureMatrix<-test[,colmaxselection]
#transpose matrix
t_max_normout_featureMatrix<-t(max_normout_featureMatrix)
# Create heatmap, you can change the name of the pdf file
pheatmap(t_max_normout_featureMatrix,color=inferno(21),scale="row",center='TRUE',cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_methode="complete",cellwidth = 2, cellheight = 2, fontsize=1, filename="cells_heatmap_positive.pdf")


#For negative data heatmap 
t_max_normout_featureMatrix<-t(test[-1])

colnames(t_max_normout_featureMatrix)<-test$sample

#negative
pheatmap(t_max_normout_featureMatrix,color=inferno(9),scale="row",center='TRUE',cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_methode="complete",cellwidth = 2, cellheight = 2, fontsize=1, filename="cells_heatmap_negative.pdf")

