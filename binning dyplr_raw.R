library(tidyverse)
library("corrplot")
library("PerformanceAnalytics")
library("FactoMineR")
library("factoextra")
library(data.table)
library(ROCR)

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)


library(plyr)
library(readr)
#Mac
mydir <- "raw"
myfiles <- list.files(path=mydir, pattern="*.csv", full.names=TRUE)
dat_csv <- ldply(myfiles, read_csv)
cells<-dat_csv
#Windows

temp = list.files(pattern="csv$")

cells <- lapply(temp, read_csv) %>% 
  bind_rows() 
#Add id as unique identifier for all every spectrum
cells$id<-gsub(" ", "_", paste(cells$sample,cells$spectranum))

#cells$id <- paste(cells$sample, "_", cells$spectranum)
#remove all white spaces
#cells$id <- gsub('\\s+', '', cells$id)
cells<-filter(cells, between(mz, 782.5616, 782.5695)) 
########get the spectrum number from the spectrum with a signal for adipo
cellsidx<-cells$id

#####Get complete spectrum from the spectrum with a adipo signal
cellspectrum<-tibble()
for (i in 1:length(cellsidx)) {
  idx<-cellsidx[i]
  spectraadipoidx<-filter(cells, spectranum == idx)
  cellspectrum = rbind(cellspectrum,spectraadipoidx)
}
####round intensity to 4 digits

cellspectrum$mz<-round(adipospectrum$mz,digits = 4)
max_int<- cellspectrom %>% 
 dplyr::group_by(id) %>% 
dplyr::summarise(int_max = max(intensity)) 

normspectramatrix<-cellspectrum %>% inner_join(max_int, by="spectranum")

normspectramatrix$int_norm <- normspectramatrix$intensity/normspectramatrix$int_max




res <- 0.007
mzmin<- min(cells$mz) - res
mzmax <- max(cells$mz) + res
len <-(mzmax - mzmin) / res
mz <- seq( from=mzmin, to=mzmax, length=len )
binned_data<-matrix()

cells$mz<-round(cells$mz, digits = 4)
bins2<-cells %>%
  #group_by(sample)%>%
  mutate(bin = cut(mz, breaks = seq(mzmin, mzmax, by = 0.001))) 
spectrumpersample<-plyr::count(cells$id)
n<-nrow(spectrumpersample)
cutoff <- 0.2 *n

frequency<-bins2%>%dplyr::count(bin)

binnedmatrix<-bins2 %>% inner_join(frequency, by="bin")
#Remove bins not present in 20% of the samples
binnedmatrix<-filter(binnedmatrix, n>cutoff)
binnedmatrix<-filter(binnedmatrix, mz>400)

binnedmatrix<-separate(binnedmatrix,bin,c("start","end"),sep = ",")

binnedmatrix$start<-as.numeric(str_remove(binnedmatrix$start, "\\("))

binnedmatrix$end<-as.numeric(str_remove(binnedmatrix$end, "\\]"))

binnedmatrix$averagemass<-(averagemass = (binnedmatrix$start + binnedmatrix$end)/2)

features<-select(binnedmatrix,averagemass,intensity,sample,spectranum)

featuresmz<-dplyr::count(features, averagemass)
#res.pca = PCA(features[,1:2], scale.unit=TRUE, ncp=5, graph=T)
#Matrix for features
samplesn<-as.numeric(nrow(featuresmz))
featuresn<-as.numeric(nrow(spectrumpersample))

#featurematrix<-matrix(nrow=samplesn,ncol=featuresn)
#featurematrix
#counter<-1
#for (sample in 1:samplesn) {
 # if (features$spectranum==sample){
  #  print(sample)
   # counter+1
  #}
#}


features$V1 <- paste("X", features$averagemass, sep = "")
features<-features[-1]
features<-dplyr::rename(features, "averagemass"=V1)

#Transpose the matrix goup by sample, if the bin has two values it averages the intensity
test2<-features %>% group_by(sample) %>%
  pivot_wider(id_cols = c(spectranum,sample), names_from = averagemass, values_from = intensity,values_fn = list(intensity= mean),values_fill = 0)

#standarize values if necesary
#test2<-test2 %>%   mutate_at(vars(-("sample"),-("spectranum")),scale2)


res.pca = PCA(test2[,3:1173], scale.unit = TRUE)
#plot.PCA(res.pca, axes=c(1, 2), choix="ind")
test2$sample<-factor(test2$sample)
#fviz_pca(res.pca, geom="point",habillage= test2$sample)
tiff(filename = "All_PANC_pos.tiff",
     width = 30, height = 20, units = "cm",res=1500)
fviz_pca_ind(res.pca, label="none", habillage=test2$sample,
             addEllipses=TRUE, ellipse.level=0.95)
dev.off()
#print(p)

cutoff <- nrow(features)



PCAprcomp <- prcomp(test2[,3:1173], scale.=TRUE, center = TRUE)
p <- pca(test2[,3:1173], scale.=TRUE, center = TRUE,removeVar = 0.1)
sd <- PCAprcomp$sdev
loadings <- PCAprcomp$rotation
rownames(loadings) <- colnames(test2[,3:1173])
scores <- PCAprcomp$x

var <- sd^2
varPercent <- var/sum(var) * 100

#Plot components and the % of variance explained by each one
barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')



varPercent[1:3]
sum(varPercent[1:3])

loadings
sqrt(1/ncol(test2[,3:50])) # cutoff for 'important' loadings
###10 PCA
PCA1<-sort(loadings[,1])
PCA1<-data.frame(lapply(PCA1, type.convert), stringsAsFactors=FALSE)
PCA2<-sort(loadings[,2])
PCA2<-data.frame(lapply(PCA2, type.convert), stringsAsFactors=FALSE)
#0.09759001
#biplot(scores[, 1:2], loadings[, 1:2], cex=0.09759001)
#biplot(scores[, 1:2], loadings[, 1:2], cex=0.7, pc.biplot=TRUE)


########Other method
PCAprcomp <- prcomp(test2[,3:1175], scale.=TRUE)

sd <- PCAprcomp$sdev
loadings <- PCAprcomp$rotation
rownames(loadings) <- colnames(test2[,3:1175])
scores <- PCAprcomp$x

var <- sd^2
varPercent <- var/sum(var) * 100

#Plot components and the % of variance explained by each one
barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')

abline(h=1/ncol(test2[,3:1175])*100, col='red')


varPercent[1:10]
sum(varPercent[1:10])

loadings
sqrt(1/ncol(test2[,3:1175])) # cutoff for 'important' loadings
###10 PCA
PCA1<-sort(loadings[,1])
PCA1<-data.frame(lapply(PCA1, type.convert), stringsAsFactors=FALSE)
PCA2<-sort(loadings[,2])
PCA2<-data.frame(lapply(PCA2, type.convert), stringsAsFactors=FALSE)
#0.09759001
biplot(scores[, 1:2], loadings[, 1:2], cex=0.02919786)
biplot(scores[, 1:2], loadings[, 1:2], cex=0.7, pc.biplot=TRUE)

# Estimate preprocessing parameters

preproc.param <- train %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train)
test.transformed <- preproc.param %>% predict(test)

library(MASS)
model <- lda(sample~., data = train.transformed)
model
# Make predictions
predictions <- model %>% predict(test.transformed)

p2 <- predict(model, test.transformed)$class
# Model accuracy
mean(predictions$class==test.transformed$sample)
lda.data <- cbind(train.transformed, predict(model)$x)

tab1 <- table(Predicted = p2, Actual = test.transformed$sample)

sum(diag(tab1))/sum(tab1)

wdbc_raw.lda.predict <- predict(model, newdata = test.transformed)

lda_plot <- cbind(train.transformed, predict(model)$x)

ggplot(lda_plot, aes(LD1, LD2)) +
  geom_point(aes(color = Species))

library(randomForest)
require(caTools)
library(caret)
set.seed(222)
rfdata<-test2[,2:225]
rfdata$sample<-as_factor(rfdata$sample)
ind <- sample(2, nrow(rfdata), replace = TRUE, prob = c(0.7, 0.3))
train <- rfdata[ind==1,]
test <- rfdata[ind==2,]


dim(train)
dim(test)


rf <- randomForest(
  sample ~ .,
  data=rfdata
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
tiff(filename = "All_cells_gini.tiff",
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
Labels<-test$sample
test$label<-as.factor(test$sample)
## for plotting
colors = rainbow(length(unique(test$sample)))
names(colors) = unique(test$sample)

## Executing the algorithm on curated data
tsne <- Rtsne(test2[,3:225], dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
#exeTimeTsne<- system.time(Rtsne(train[,-1,-2], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500))

#tiff(filename = "all_ady_tSNe.tiff",
 #    width = 30, height = 20, units = "cm",res=1500)
## Plotting
#plot(tsne$Y, t='n', main="tsne")
#text(tsne$Y, labels=test2$sample, col=colors[test2$sample])
#text(tsne$Y,  col=colors[test2$sample])

#dev.off()
#https://datavizpyr.com/how-to-make-tsne-plot-in-r/
#The tSNE result object contains two tSNE components that we are interested in. We can extract the components and save it in a dataframe.
#get labels
spectranumcol<-test2$spectranum
samplecol<-test2$sample
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
ggsave("tSNE_plot_example100.png")







####Select features by random forest and do PCA and tSNE again for positive

test2<-test2 %>% 
  mutate_if(is.numeric, round,digits=4)

#transpose the matrix including column names as a column called featur
trans_test2<- data.table::transpose(test2,keep.names = "feature")
trans_test2<-trans_test2[-2,]
trans_test2<-as_tibble(trans_test2)
binnedft<-trans_test2 %>% inner_join(selected_features, by="feature")
binnedft<-binnedft %>% select(feature,MeanDecreaseGini, everything())
spectranumtag<-test2$spectranum
sampletag<-test2$sample
binnedft<-arrange(binnedft, desc(MeanDecreaseGini))
importantft<-binnedft [1:200,]
trans_importantft<- data.table::transpose(importantft)
colnames(trans_importantft) <- trans_importantft[1, ]
trans_importantft<-trans_importantft[-c(1,2),]

#trans_importantft<-as.numeric(trans_importantft)

#trans_importantft %>% as.integer

trans_importantft<-trans_importantft %>% 
  mutate_all(as.numeric)

trans_importantft$sample<-sampletag
trans_importantft$spectranum<-spectranumtag
trans_importantft<-trans_importantft %>% select(spectranum,sample, everything())
trans_importantft<-as_tibble(trans_importantft)

# transposed_ft<-t(test2)
# transposed_ft <- transposed_ft %>% 
#   data.frame() %>% 
#   mutate(feature = row.names(.))
# transposed_ft<-as_tibble(transposed_ft)
# transposed_ft<-transposed_ft[-2,]
# colnames(transposed_ft) <- transposed_ft[1, ]
# transposed_ft<-transposed_ft[-1,]
# transposed_ft<-transposed_ft %>% select(feature, everything())
# binnedft<-features_invivo %>% join(feat_imp_df, by="feature")
# binnedft<-binnedft %>% select(MeanDecreaseGini, everything())
# binnedft<-arrange(binnedft, desc(MeanDecreaseGini))
# binnedft$feature <-  sub("X", "", binnedft$feature)
# write_csv(binnedft,file = "features_invivo.csv")
# 
# 

## Curating the database for analysis with both t-SNE and PCA
Labels<-trans_importantft$sample
trans_importantft$label<-as.factor(trans_importantft$sample)
## for plotting
colors = rainbow(length(unique(trans_importantft$sample)))
names(colors) = unique(trans_importantft$sample)

## Executing the algorithm on curated data
tsne <- Rtsne(trans_importantft[,1:200], dims = 2, perplexity=10, verbose=TRUE, max_iter = 1000)

#https://datavizpyr.com/how-to-make-tsne-plot-in-r/
#The tSNE result object contains two tSNE components that we are interested in. We can extract the components and save it in a dataframe.
#get labels
spectranumcol<-trans_importantft$spectranum
samplecol<-trans_importantft$sample
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
ggsave("tSNE_plot_200fts.png")


res.pca = PCA(trans_importantft[,1:200], scale.unit = TRUE)
#plot.PCA(res.pca, axes=c(1, 2), choix="ind")
trans_importantft$sample<-factor(trans_importantft$sample)
#fviz_pca(res.pca, geom="point",habillage= test2$sample)
tiff(filename = "All_cells_200ft.tiff",
     width = 30, height = 20, units = "cm",res=1500)
fviz_pca_ind(res.pca, label="none", habillage=trans_importantft$sample,
             addEllipses=TRUE, ellipse.level=0.95)
dev.off()
#print(p)

cutoff <- nrow(features)



PCAprcomp <- prcomp(test2[,3:322], scale.=TRUE, center = TRUE)

sd <- PCAprcomp$sdev
loadings <- PCAprcomp$rotation
rownames(loadings) <- colnames(test2[,3:322])
scores <- PCAprcomp$x

var <- sd^2
varPercent <- var/sum(var) * 100

#Plot components and the % of variance explained by each one
barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')



varPercent[1:5]
sum(varPercent[1:5])

loadings
sqrt(1/ncol(test2[,3:322])) # cutoff for 'important' loadings
###10 PCA
PCA1<-sort(loadings[,1])
PCA1<-data.frame(lapply(PCA1, type.convert), stringsAsFactors=FALSE)
PCA2<-sort(loadings[,2])
PCA2<-data.frame(lapply(PCA2, type.convert), stringsAsFactors=FALSE)


####LDA

# Estimate preprocessing parameters
preproc.param <- train %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train)
test.transformed <- preproc.param %>% predict(test)

# library(MASS)
# # Fit the model
# model <- lda(sample~., data = train.transformed)
# # Make predictions
# predictions <- model %>% predict(test.transformed)
# # Model accuracy
# mean(predictions$class==test.transformed$sample)
# # Predicted classes
# head(predictions$class, 6)
# # Predicted probabilities of class memebership.
# head(predictions$posterior, 6) 
# # Linear discriminants
# head(predictions$x, 3) 
# 
# lda.data <- cbind(train.transformed, predict(model)$x)
# ggplot(lda.data, aes(LD1, LD2)) +
#   geom_point(aes(color = sample))
# library(mda)
# # Fit the model
# modelmda <- mda(sample~., data = train.transformed)
# 
# # Make predictions
# predicted.classesmda <- modelmda %>% predict(test.transformed)
# # Model accuracy
# mean(predicted.classesmda == test.transformed$sample)
# 
# ####QDA
# # Fit the model
# modelqda <- qda(sample~., data = train.transformed)
# 
# # Make predictions
# predictionsqda <- modelqda %>% predict(test.transformed)
# # Model accuracy
# mean(predictionsqda$class == test.transformed$sample)
# datPred <- data.frame(sample=predict(model)$class,predict(model)$x)
# 
# ggplot(datPred, aes(LD1,color=sample))

# ensure the results are repeatable
set.seed(7)
# load the library
library(mlbench)
library(caret)
# load the data
#data(PimaIndiansDiabetes)
# calculate correlation matrix
correlationMatrix <- cor(test2[,3:1173])
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
# print indexes of highly correlated attributes
print(highlyCorrelated)

set.seed(7)
# load the library
library(mlbench)
library(caret)
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(sample~., data=test2, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)



#####K means

library(ClusterR)
library(cluster)


test2_k <- test2[, c(-1,-2)]

set.seed(240) # Setting seed
kmeans.re <- kmeans(test2_k, centers = 6, nstart = 20)
kmeans.re

# Cluster identification for 
# each observation
kmeans.re$cluster

# Confusion Matrix
cm <- table(test2$sample, kmeans.re$cluster)
cm

# Model Evaluation and visualization
plot(test2_k[c("AsPC-1-2Neg", "Sepal.Width")])
plot(test2_k, 
     col = kmeans.re$cluster, 
     main = "K-means with 6 clusters")

## Plotiing cluster centers
kmeans.re$centers
kmeans.re$centers[, c("Sepal.Length", "Sepal.Width")]

# cex is font size, pch is symbol
points(kmeans.re$centers[, c("Sepal.Length", "Sepal.Width")], 
       col = 1:3, pch = 8, cex = 3) 

## Visualizing clusters
y_kmeans <- kmeans.re$cluster
clusplot(iris_1[, c("Sepal.Length", "Sepal.Width")],
         y_kmeans,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         labels = 2,
         plotchar = FALSE,
         span = TRUE,
         main = paste("Cluster iris"),
         xlab = 'Sepal.Length',
         ylab = 'Sepal.Width')
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



test2<-test2 %>% 
  mutate_if(is.numeric, round,digits=4)

#transpose the matrix including column names as a column called featur
trans_test2<- data.table::transpose(test2,keep.names = "feature")
trans_test2<-trans_test2[-2,]
trans_test2<-as_tibble(trans_test2)

i=1
while(i<length(selected_features)){
  print(selected_features[i])
  i=i+1
}

####Get selected features from the table
keep<- c("sample", "spectranum", selected_features)

df <- test2[keep]

library(Rtsne)

## Curating the database for analysis with both t-SNE and PCA
Labels<-df$sample
df$label<-as.factor(df$sample)
## for plotting
colors = rainbow(length(unique(df$sample)))
names(colors) = unique(df$sample)


tsne <- Rtsne(df[,3:127], dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

spectranumcol<-df$spectranum
samplecol<-df$sample
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
ggsave("tSNE_plot_boruta_neg_perplexity100.png")

res.pca = PCA(df[,3:126], scale.unit = TRUE)
#plot.PCA(res.pca, axes=c(1, 2), choix="ind")
df$sample<-factor(df$sample)
#fviz_pca(res.pca, geom="point",habillage= test2$sample)
tiff(filename = "boruta_neg_PCA.tiff",
     width = 30, height = 20, units = "cm",res=1500)
fviz_pca_ind(res.pca, label="none", habillage=helat$sample,
             addEllipses=TRUE, ellipse.level=0.95)
dev.off()
#print(p)

correlationMatrix <- cor(df[,3:126])
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
# print indexes of highly correlated attributes
print(highlyCorrelated)


count(test2, 'sample')


helatest <- filter(df,sample=="cfpanc-1"| sample=="HeLa-1Neg")
res.pca = PCA(helatest[,3:126], scale.unit = TRUE)
helatest$sample<-factor(helatest$sample)
#fviz_pca(res.pca, geom="point",habillage= test2$sample)
tiff(filename = "PCA_Neg_HeLavsCFPAC.tiff",
     width = 30, height = 20, units = "cm",res=1500)
fviz_pca_ind(res.pca, label="none", habillage=helatest$sample,
             addEllipses=TRUE, ellipse.level=0.95)
dev.off()
#print(p)


helatest <- filter(df,sample=="panc-10-05Neg"| sample=="AsPC-1-2Neg"| sample=="BxPC-3-1Neg")
res.pca = PCA(helatest[,3:126], scale.unit = TRUE)
helatest$sample<-factor(helatest$sample)
#fviz_pca(res.pca, geom="point",habillage= test2$sample)
tiff(filename = "PCA_Neg_pancvsAsPCvsBxPC3.tiff",
     width = 30, height = 20, units = "cm",res=1500)
fviz_pca_ind(res.pca, label="none", habillage=helatest$sample,
             addEllipses=TRUE, ellipse.level=0.95)
dev.off()
#print(p)



if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

  p <- pca(test2[,3:225], removeVar = 0.1)

  screeplot(p, axisLabSize = 18, titleLabSize = 22)
  pairsplot(p)
  plotloadings(p, labSize = 3)

  pairsplot(p,
            components = getComponents(p, c(4,33,11,1)),
            triangle = FALSE,
            hline = 0, vline = 0,
            pointSize = 0.8,
            gridlines.major = FALSE, gridlines.minor = FALSE,
            colby = 'ER',
            title = 'Pairs plot', titleLabSize = 22,
            axisLabSize = 14, plotaxes = TRUE,
            margingaps = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
  
  cor.mat <- round(cor(df[,3:126]),2)
  library("corrplot")
  corrplot(cor.mat, type="upper", order="hclust", 
           tl.col="black", tl.srt=45, labSize = 1)
  
  
#####RF after boruta
  
  set.seed(222)
  rfdata<-df[,-2]
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
  
  
  preproc.param <- train %>% 
    preProcess(method = c("center", "scale"))
  # Transform the data using the estimated parameters
  train.transformed <- preproc.param %>% predict(train)
  test.transformed <- preproc.param %>% predict(test)
  
  library(MASS)
  model <- lda(sample~., data = train.transformed)
  model
  # Make predictions
  predictions <- model %>% predict(test.transformed)
  
  p2 <- predict(model, test.transformed)$class
  # Model accuracy
  mean(predictions$class==test.transformed$sample)
  lda.data <- cbind(train.transformed, predict(model)$x)
  
  tab1 <- table(Predicted = p2, Actual = test.transformed$sample)
  
  sum(diag(tab1))/sum(tab1)
  
  wdbc_raw.lda.predict <- predict(model, newdata = test.transformed)
  
  lda_plot <- cbind(train.transformed, predict(model)$x)
  
  ggplot(lda_plot, aes(LD1, LD2)) +
    geom_point(aes(color = Species))
  
  
