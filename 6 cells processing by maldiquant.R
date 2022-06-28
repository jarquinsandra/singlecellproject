
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
path <- "~/Sandra_MS_data/maldiquant_cells/neg_raw"
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
#find a way to add labels
## Normalize spectra
#remove all features with more than 20% values by class
#remove all features with more than 20% values by class
features<-as_tibble(featureMatrix,.name_repair = ~ make.names(., unique = TRUE))
#Remove columns that have more than 20% 0 values
zerocutoff<-nrow(features)*0.2

features[colSums(features > 0) <= zerocutoff]<- NULL
apply( features , 2 , function(x) sum ( x == 0 ) )

#select all features in original tibble
features<-features%>% select(all_of(ft))
#1negative 
#spectra[[1998]]@metaData[["file"]]
ASPC<-replicate(1268, "ASPC")
BxPC<-replicate(1264, "BxPC")
CFPAC<-replicate(1337, "CFPAC")
HeLa<-replicate(1264, "HeLa")
Panc1<-replicate(1484, "Panc1")
Panc10<-replicate(1267, "Panc10")

alllabels<-c(ASPC,BxPC,CFPAC,HeLa,Panc1,Panc10)
#add labels
features$sample<-alllabels
features<-tibble::rowid_to_column(features, "id")
#moving sample to first column
features <- features %>%
  select(id,sample, everything())

#filter rows by mass to select cells Negative
cells<-features[!(features$X804.574332386052==0|features$X830.606791452377)==0,]



library(randomForest)
require(caTools)
library(caret)
set.seed(222)
rfdata<-features[,2:115]
rfdata$sample<-as_factor(rfdata$sample)
ind <- sample(2, nrow(rfdata), replace = TRUE, prob = c(0.7, 0.3))
train <- rfdata[ind==1,]
test <- rfdata[ind==2,]


dim(train)
dim(test)


rf <- randomForest(
  sample ~ .,
  data=test
)
#Prediction & Confusion Matrix – train data
p1 <- predict(rf, train)
confusionMatrix(p1, train$sample)
#Prediction & Confusion Matrix – test data
p2 <- predict(rf, test)
confusionMatrix(p2, test$sample)


hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")
#Variable Importance
tiff(filename = "All_adypocites_gini.tiff",
     width = 30, height = 20, units = "cm",res=1500)
varImpPlot(rf,
           sort = T,
           n.var = 30,
           main = "Top 30 - Variable Importance")
dev.off()
importance(rf)

capture.output(rf, file = "rf.txt", append = TRUE)

feat_imp_df <- importance(rf) %>% 
  data.frame() %>% 
  mutate(feature = row.names(.)) 


library(Rtsne)

## Curating the database for analysis with both t-SNE and PCA
Labels<-features$sample
features$sample<-as.factor(features$sample)
## for plotting
colors = rainbow(length(unique(features$sample)))
names(colors) = unique(features$sample)

## Executing the algorithm on curated data
tsne <- Rtsne(features[,3:115], dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
#get labels
spectranumcol<-features$id
samplecol<-features$sample
tSNE_df <- tsne$Y %>% 
  as.data.frame() %>%
  dplyr::rename(tsne1="V1",
                tsne2="V2") %>%
  mutate(sample=samplecol, spectranum=spectranumcol)


tSNE_df %>%
  ggplot(aes(x = tsne1, 
             y = tsne2,
             color = sample))+
  geom_point()+
  theme(legend.position="bottom")
ggsave("tSNE_all_maldiquant.png")

#Boruta
library(Boruta)
boruta <- Boruta(sample ~ ., data = rfdata, doTrace = 1, maxRuns = 100)

print(boruta)


plot(boruta, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta$ImpHistory),function(i)
  boruta$ImpHistory[is.finite(boruta$ImpHistory[,i]),i])
names(lz) <- colnames(boruta$ImpHistory)
Labels <- sort(sapply(lz,median))

axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta$ImpHistory), cex.axis = 0.7)
#Select features
selected_features<-getSelectedAttributes(boruta, withTentative = F)
boruta.df <- attStats(boruta)
class(boruta.df)



features<-features %>% 
  mutate_if(is.numeric, round,digits=4)

#transpose the matrix including column names as a column called featur
trans_test2<- data.table::transpose(features,keep.names = "feature")
trans_test2<-trans_test2[-2,]
trans_test2<-as_tibble(trans_test2)

####Get selected features from the table
keep<- c("sample", "id", selected_features)

df <- features[keep]

#get only confirmed features
boruta.df<-filter(boruta.df, decision=='Confirmed')
mz<-rownames(boruta.df)
boruta.df$mz<-mz
boruta.df$mz<-gsub("X","",boruta.df$mz)
boruta.df$mz<-as.numeric(boruta.df$mz)
boruta.df<-as_tibble(boruta.df)
write.csv(boruta.df, file='featuresimportance.csv') 

####clases file from python
names<-colnames(clusters)
names<-sub('X', '', names)
names<-names[3:161]
names<-as.numeric(names)
names<-round(names,digits=4)
names<-as.character(names)
#names<-names[names > 400]
#add X to characters
names<-paste0("X", names)
labels<-c("sample","id",names,"hdb_prob", "hdb_cluster")
colnames(clusters)<-labels



