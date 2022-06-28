
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
#find a way to add labels
## Normalize spectra
#remove all features with more than 20% values by class
features<-as_tibble(featureMatrix,.name_repair = ~ make.names(., unique = TRUE))
#Remove columns that have more than 20% 0 values
zerocutoff<-nrow(features)*0.2

features[colSums(features > 0) <= zerocutoff]<- NULL
apply( features , 2 , function(x) sum ( x == 0 ) )

####test for automatic label extraction
i=1
labelst<-tibble()
for (i in 1:length(spectra)) {
  label<-spectra[[i]]@metaData[["file"]]
  labelst<-rbind(labelst,label)
  i=+1
}

labels <- mutate_if(labelst, 
                is.character, 
                str_replace_all, 
                pattern = "/Users/sandramartinez/Documents/MSData/mice_wtkras_20220405/raw/", 
                replacement = "")
labels <- mutate_if(labels, 
                is.character, 
                str_replace_all, 
                pattern = ".mzML", 
                replacement = "")

#3679 KC + 3933 WT 
#KC<-replicate(3679, "KC")
#WT<-replicate(3933, "WT")
#alllabels<-c(KC,WT)
#add labels
features$sample<-labels
features<-tibble::rowid_to_column(features, "id")
#moving sample to first column
features <- features %>%
  select(id,sample, everything())

#filter rows by mass to select cells
cells<-features[!(features$X782.565591685428==0|features$X760.581056660418)==0,]
#df2<-df1[!(df1$Name=="George" | df1$Name=="Andrea"),]

library(randomForest)
require(caTools)
library(caret)
set.seed(222)
rfdata<-cells[,2:765]
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
Labels<-cells$sample
cells$sample<-as.factor(cells$sample)
## for plotting
colors = rainbow(length(unique(cells$sample)))
names(colors) = unique(cells$sample)

## Executing the algorithm on curated data
tsne <- Rtsne(cells[,3:765], dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
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
  geom_point()+
  theme(legend.position="bottom")
ggsave("tSNE_all_maldiquant_cellsfilter.png")

#Boruta
library(Boruta)
boruta <- Boruta(sample ~ ., data = rfdata, doTrace = 1, maxRuns = 100)

print(boruta)
tiff(filename = "borutatest.tiff",
     width = 55, height = 25, units = "cm",res=1500)

plot(boruta, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta$ImpHistory),function(i)
  boruta$ImpHistory[is.finite(boruta$ImpHistory[,i]),i])
names(lz) <- colnames(boruta$ImpHistory)
Labels <- sort(sapply(lz,median))
dev.off()


axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta$ImpHistory), cex.axis = 0.7)
#Select cells
selected_cells<-getSelectedAttributes(boruta, withTentative = F)
boruta.df <- attStats(boruta)
class(boruta.df)



cells<-cells %>% 
  mutate_if(is.numeric, round,digits=4)

#transpose the matrix including column names as a column called featur
trans_test2<- data.table::transpose(cells,keep.names = "feature")
trans_test2<-trans_test2[-2,]
trans_test2<-as_tibble(trans_test2)

####Get selected cells from the table
keep<- c("sample", "id", selected_cells)

df <- cells[keep]

#get only confirmed cells
boruta.df<-filter(boruta.df, decision=='Confirmed')
mz<-rownames(boruta.df)
mz<-sub('X', '', mz)
mz<-as.numeric(mz)
mz<-round(mz,digits=4)
mz<-as.character(mz)
mz<-paste0("X", mz)
boruta.df$mz<-mz
boruta.df$mz<-gsub("X","",boruta.df$mz)
boruta.df$mz<-as.numeric(boruta.df$mz)
boruta.df<-as_tibble(boruta.df)

#Remove weird extra numbers after the decimal point

names<-colnames(df)
names<-sub('X', '', names)
names<-names[3:166]
names<-as.numeric(names)
names<-round(names,digits=4)
names<-as.character(names)
#names<-names[names > 400]
#add X to characters
names<-paste0("X", names)
labels<-c("sample","id",names)
colnames(df)<-labels
write.csv(boruta.df, file='cellsimportance.csv') 
write.csv(df, file='cellsborutaselectedfeatures_filtered.csv') 


####clases file from python
# names<-colnames(clusters)
# names<-sub('X', '', names)
# names<-names[3:161]
# names<-as.numeric(names)
# names<-round(names,digits=4)
# names<-as.character(names)
#names<-names[names > 400]
#add X to characters
# names<-paste0("X", names)
# labels<-c("sample","id",names,"hdb_prob", "hdb_cluster")
# colnames(clusters)<-labels



