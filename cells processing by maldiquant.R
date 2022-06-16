
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

## get names of all .mzML files in current working directory
path <- "~/Documents/MSData/mice_wtkras_20220405/raw"
spectra<-importMzMl(path = path)
#normalization to TIC

spectra <- calibrateIntensity(spectra, method="TIC")

## Peak picking, you can change SNR however for orbitrap data the noise is very low

peaks<-detectPeaks(spectra, SNR=3)
## Peak alignment using the warping function of Maldiquant

warpingFunctions<-determineWarpingFunctions(peaks)
peaks<-warpMassPeaks(peaks,warpingFunctions)

#Peak binning and creation of feature matrix
peaks<-binPeaks(peaks)
featureMatrix<-intensityMatrix(peaks)

## Replace NAN values by 0
featureMatrix[is.na(featureMatrix)] <- 0
####uncomment this to get the features matrix as a csv
#write.csv(featureMatrix,"featureMatrix.csv")

## Normalize spectra

normout_featureMatrix<-normalize.rows(featureMatrix)
samplenames<-gsub(".mzML","",mzML_list)
rownames(normout_featureMatrix)<-samplenames
colnames(normout_featureMatrix)<-colnames(featureMatrix)

###########################################################################
##For heatmaps construction
## Sort by columns with maximum intensity and eliminate others

#find maximal value for each column (mz-bin)
columnmax<-apply(normout_featureMatrix, 2, max)

#sort and select the x most intense values
colmax<-sort(columnmax,decreasing=TRUE)
colmaxselection<-as.matrix(colmax[1:100])
colmaxselection=as.vector(rownames(colmaxselection))

#limiting matrix to selection
max_normout_featureMatrix<-normout_featureMatrix[,colmaxselection]

## Create heatmap

t_max_normout_featureMatrix<-t(max_normout_featureMatrix)
pheatmap(t_max_normout_featureMatrix,scale="row",cluster_cols = TRUE, cluster_rows = TRUE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_methode="complete",cellwidth = 10, cellheight = 10, fontsize=5, filename="metabolic-heatmapSNR0_all_100.pdf")

###########################################################################

####for PCA construction

df<-as.data.frame(normout_featureMatrix)
flow<-c("1.5","1.5","1.5","1.0","1.0","1.0","0.5","0.5","0.5","1.5","1.5","1.5","1.0","1.0","0.5","0.5","0.5")
lid<-c("sealed","sealed","sealed","sealed","sealed","sealed","sealed","sealed","opened","opened","opened","opened","opened","opened","opened","opened")
voltage<-c("2.44","3.00","3.44","2.44","3.00","3.44","2.44","3.00","3.44","2.44","3.00","3.44","2.44","3.44","2.44","3.00","3.44")
samples <-c("a","b","c","d","e","f","g","h","i","j","k","l","m","o","p","q","r")

#Here you calculate the pca components
res.pca <- prcomp(df,  scale = TRUE)
#Plot PCA components change habillage to specify the parameter you want to evaluate and the tiff file
tiff(filename = "N_flow.tiff",
     width = 15, height = 10, units = "cm",res=1500)
fviz_pca_ind(res.pca, label="none", habillage=flow, geom = c("point"),
             addEllipses=TRUE, ellipse.type = "convex", ellipse.level=.8)
  #guides(fill=guide_legend(title="Flow (L/min)"))
dev.off()

##########PCAS with X most intense features


df<-as.data.frame(max_normout_featureMatrix)
flow<-c("1.5","1.5","1.5","1.0","1.0","1.0","0.5","0.5","0.5","1.5","1.5","1.5","1.0","1.0","0.5","0.5","0.5")
lid<-c("sealed","sealed","sealed","sealed","sealed","sealed","sealed","sealed","opened","opened","opened","opened","opened","opened","opened","opened")
voltage<-c("2.44","3.00","3.44","2.44","3.00","3.44","2.44","3.00","3.44","2.44","3.00","3.44","2.44","3.44","2.44","3.00","3.44")
samples <-c("a","b","c","d","e","f","g","h","i","j","k","l","m","o","p","q","r")

#Here you calculate the pca components
res.pca <- prcomp(df,  scale = TRUE)
#Plot PCA components change habillage to specify the parameter you want to evaluate and the tiff file
tiff(filename = "100_flow.tiff",
     width = 15, height = 10, units = "cm",res=1500)
fviz_pca_ind(res.pca, label="none", habillage=flow, geom = c("point"),
             addEllipses=TRUE, ellipse.type = "convex", ellipse.level=.8)
#guides(fill=guide_legend(title="Flow (L/min)"))
dev.off()

###Test for specific masses to evaluate how a compound changes in a group of samples in normalized data

normout_featureMatrix<-as.tibble(normout_featureMatrix)
stats <- data.frame( )
counter<-1
pdf("features_normalized.pdf")
while (counter<ncol(normout_featureMatrix)) {
  
  plotdt<-normout_featureMatrix%>%dplyr::pull(counter)
  
  boxplot(plotdt~lid,
          main= names(normout_featureMatrix[counter]),
          xlab="Evaluated parameter",
          ylab="Normalized intensity",
          col="orange",
          border="brown"
  )
  start= plotdt[1:(length(plotdt)/2)]
  end= plotdt[((length(plotdt)/2)+1):length(plotdt)]
  #res <- wilcox.test(start,end,paired=TRUE)
  res<-kruskal.test(list(start,end))
  print(res)
  
  stats[counter,1] <- res[1]
  stats[counter,2] <- res[2]
  stats[counter,3] <- res[3]
  stats[counter,4]<-names(normout_featureMatrix[counter])
  counter=counter+1
}
dev.off()
#selects significatively different features with a p value smaller than 0.001, this parameter can be modified
sigfeatures<-subset(stats, p.value<0.001)
sigfeaturesmasses<-sigfeatures[4]


counter<-1
pdf("significative_features_normalized.pdf")
while (counter<nrow(sigfeatures)) {
  
  mass <- sigfeaturesmasses[counter,1]
  plotdt<-normout_featureMatrix %>% select(mass)
  my_data <- data.frame(feature=plotdt,
                        lid=lid)
  colnames(my_data)[1]<-"feature"
  group_by(my_data, lid) %>%
    summarise(
      count = n(),
      median = median(feature, na.rm = TRUE),
    )
  my_data<-as.tibble(my_data)
  boxplot(my_data$feature~my_data$lid,
          main= sigfeaturesmasses[counter,1],
          xlab="Evaluated parameter",
          ylab="Normalized intensity",
          col="orange",
          border="brown"
  )
  
  counter=counter+1
}
dev.off()
##################################
#Same in non normalized data
###Test for specific masses to evaluate how a compound changes in a group of samples in normalized data

featureMatrix<-as.tibble(featureMatrix)
stats <- data.frame( )
counter<-1
pdf("features.pdf")
while (counter<ncol(featureMatrix)) {
  
  plotdt<-featureMatrix%>%dplyr::pull(counter)
  
  boxplot(plotdt~lid,
          main= names(normout_featureMatrix[counter]),
          xlab="Evaluated parameter",
          ylab="Absolut intensity",
          col="orange",
          border="brown"
  )
  start= plotdt[1:(length(plotdt)/2)]
  end= plotdt[((length(plotdt)/2)+1):length(plotdt)]
  #res <- wilcox.test(start,end,paired=TRUE)
  res<-kruskal.test(list(start,end))
  print(res)
  
  stats[counter,1] <- res[1]
  stats[counter,2] <- res[2]
  stats[counter,3] <- res[3]
  stats[counter,4]<-names(normout_featureMatrix[counter])
  counter=counter+1
}
dev.off()
#selects significatively different features with a p value smaller than 0.05, this parameter can be modified
sigfeatures<-subset(stats, p.value<0.05)



counter<-1
pdf("significative_features.pdf")
while (counter<nrow(sigfeatures)) {
  
  mass <- sigfeaturesmasses[counter,1]
  plotdt<-normout_featureMatrix %>% select(mass)
  my_data <- data.frame(feature=plotdt,
                        lid=lid)
  colnames(my_data)[1]<-"feature"
  group_by(my_data, lid) %>%
    summarise(
      count = n(),
      median = median(feature, na.rm = TRUE),
    )
  my_data<-as.tibble(my_data)
  boxplot(my_data$feature~my_data$lid,
          main= sigfeaturesmasses[counter,1],
          xlab="Evaluated parameter",
          ylab="Absolute intensity",
          col="orange",
          border="brown"
  )
  
  counter=counter+1
}
dev.off()
