
## load Sebastian Gibb's MALDIquant and MALDIquantForeign libraries
## Gibb S, Strimmer K. MALDIquant: a versatile R package for the analysis of mass spectrometry data. 
## Bioinformatics. 2012 Sep 1;28(17):2270-1. doi:10.1093/bioinformatics/bts447. 
## Epub 2012 Jul 12. PubMed PMID: 22796955.
#Script modified and complemented to analyze differences in features between two setup conditions for an ambient ionization source

library("MALDIquant")
library("MALDIquantForeign")
library("vegetarian")
library("pheatmap")
library("corrplot")
library("PerformanceAnalytics")
library("FactoMineR")
library("factoextra")
library(tidyverse)
library("ggpubr")
library(RVAideMemoire)
# reduce matrix size, using a summarising function (default, mean)
## get names of all .mzML files in current working directory
#windows
path <- "~/Sandra_MS_data/technicalnote/positive"
#mac
path <-"/Users/sandramartinez/Downloads/drive-download-20220715T123036Z-001"
#import spectrum
spectra<-importMzMl(path = path)
#normalization to TIC
spectra <- calibrateIntensity(spectra, method="TIC")
## Peak picking, with signal to noise filtering
peaks<-detectPeaks(spectra, SNR=3)
## Peak alignment using the warping function of Maldiquant
warpingFunctions<-determineWarpingFunctions(peaks,allowNoMatches=TRUE)
peaks<-warpMassPeaks(peaks,warpingFunctions)
# binPeaks tolerance is expressed in delta mass / mass, we want it in ppm:
#tol.ppm <- 5
#tol.maldiquant <- tol.ppm / 1e5

# The new common masses will overwrite the original ones. The intensities remain unmodified
#peaks <- binPeaks(spectra, method = 'strict', tolerance = tol.ppm)
#Peak binning and creation of feature matrix
peaks<-binPeaks(peaks)
featureMatrix<-intensityMatrix(peaks)

## Replace NAN values by 0
featureMatrix[is.na(featureMatrix)] <- 0
#find a way to add labels
## Normalize spectra
#remove all features with more than 20% values by class
#remove all features with more than 20% values by class
features<-as_tibble(featureMatrix,.name_repair = ~ make.names(., unique = TRUE))
#Remove columns that have more than 20% 0 values
zerocutoff<-nrow(features)*0.2

features[colSums(features > 0) <= zerocutoff]<- NULL
#apply( features , 2 , function(x) sum ( x == 0 ) )

 #1positive 
# spectra[[1998]]@metaData[["file"]]
 ASPC<-replicate(1265, "ASPC-1")
 HeLa<-replicate(1262, "HeLa")
 Panc10<-replicate(1263, "Panc 10.05")
 alllabels<-c(ASPC,HeLa,Panc10)
 
#add labels
features$sample<-alllabels
features<-tibble::rowid_to_column(features, "id")
#moving sample to first column
features <- features %>%
  select(id,sample, everything())
#filter rows by mass to select cells positive

names<-colnames(features)
names<-sub('X', '', names)
names<-names[3:ncol(features)]
names<-as.numeric(names)
names<-round(names,digits=4)
names<-as.character(names)
#add X to characters
names<-paste0("X", names)
labels<-c("sample","id",names)
colnames(features)<-labels

cells<-features[!(features$X732.5554==0|features$X760.5894)==0,]




library(Rtsne)

## Curating the database for analysis with both t-SNE and PCA
Labels<-cells$sample
cells$sample<-as.factor(cells$sample)
## for plotting
colors = rainbow(length(unique(cells$sample)))
names(colors) = unique(cells$sample)

## Executing the algorithm on curated data
tsne <- Rtsne(cells[,3:ncol(features)], dims = 2, perplexity=50, verbose=TRUE, max_iter = 5000)
#get labels
spectranumcol<-cells$sample
samplecol<-cells$id
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
ggsave("technicalnote120_3000.png")


library(plotly)
library(ggfortify)

df <- cells[3:ncol(cells)]
pca_res <- prcomp(df, scale. = TRUE)

p <- autoplot(pca_res, data = cells, colour = 'id')

ggplotly(p)


columnmax<-apply(featureMatrix, 2, max)

#sort and select the 75 most intense values
colmax<-sort(columnmax,decreasing=TRUE)
colmaxselection<-as.matrix(colmax[1:50])
colmaxselection=as.vector(rownames(colmaxselection))

#limiting matrix to selection
max_normout_featureMatrix<-featureMatrix[,colmaxselection]
#data<-- max_normout_featureMatrix

## Create heatmap

t_max_normout_featureMatrix<-t(max_normout_featureMatrix)
colnames(t_max_normout_featureMatrix)<-alllabels

pheatmap(t_max_normout_featureMatrix,scale="row",cluster_cols = TRUE, cluster_rows = TRUE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_methode="complete",cellwidth = 10, cellheight = 10, fontsize=5, filename="metabolic-heatmapSNR1_all_50norm.pdf")
#get names from the mz values
library(randomForest)
require(caTools)
library(caret)
set.seed(123)
rfdata<-cells[,2:ncol(cells)]
rfdata$sample<-as_factor(rfdata$sample)
ind <- sample(2, nrow(rfdata), replace = TRUE, prob = c(0.95, 0.05))
train <- rfdata[ind==1,]
test <- rfdata[ind==2,]


dim(train)
dim(test)

columnmax<-apply(test[2:ncol(test)], 2, max)

#sort and select the 75 most intense values
colmax<-sort(columnmax,decreasing=TRUE)
colmaxselection<-as.matrix(colmax[0:150])
colmaxselection=as.vector(rownames(colmaxselection))

#limiting matrix to selection
max_normout_featureMatrix<-test[,colmaxselection]
#data<-- max_normout_featureMatrix

## Create heatmap

t_max_normout_featureMatrix<-t(max_normout_featureMatrix)
colnames(t_max_normout_featureMatrix)<-test$sample

pheatmap(t_max_normout_featureMatrix,color=rainbow(50),scale="row",cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_methode="complete",cellwidth = 1, cellheight = 1, fontsize=1, filename="metabolic-heatmapSNR1_all_50norm.pdf")
#pheatmap(t_max_normout_featureMatrix,scale="row",cluster_cols = TRUE, cluster_rows = TRUE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_methode="complete",cellwidth = 10, cellheight = 10, fontsize=5, filename="metabolic-heatmapSNR1_all_50norm.pdf")

library(ggplot2)
library(patchwork)
# Wide to long transformation
data_for_ggplot <- as.data.frame(max_normout_featureMatrix) %>% 
  mutate(row = rownames(.)) %>% 
  tidyr::pivot_longer(-row, names_to = "col") %>%
  mutate(row = as.numeric(row), col = readr::parse_number(col))

# with geom_tile()
p1 <- ggplot(data_for_ggplot, aes(x = col, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "white", mid = "darkorange", high = "black",
    limits = c(0, 3), midpoint = 1.5, oob = scales::squish
  ) +
  labs(title = "geom_tile") +
  theme_void() +
  theme(legend.position = "none")
# with geom_raster()
p2 <- ggplot(data_for_ggplot, aes(x = col, y = row, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "white", mid = "darkorange", high = "black",
    limits = c(0, 3), midpoint = 1.5, oob = scales::squish
  ) +
  labs(title = "geom_raster") +
  theme_void() +
  theme(legend.position = "none")
p1 + p2

