
library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)
library(gplots)
library(cluster)
library(factoextra)
library(tidyverse)
library(hrbrthemes)
library(viridis)

## read the raw data from a mzML file
raw_data <- readMSData(files = "/Users/sandramartinezjarquin/Downloads/220706mice_csv/KC.mzML", mode = "onDisk")
#extract metadata and adquisition information from the samples
info<-fData(raw_data)
#Remove NA values
info2<-info[, colSums(is.na(info)) != nrow(info)]
#Unlist to convert to dataframe
info2 <- data.frame(matrix(unlist(info2), nrow=length(info2), byrow=TRUE))
#transpose tab;e
info2<-t(info2)
#column 9 in RT and 21 is spectrum
cols <- c(7:11, 21)
info2<-info2[,cols]
#colnames(info) <- c("OriginalPeaksCount","TotalIonCurrent","basePeakMZ","basePeakMZ","basePeakIntensity","RetTime")
#info<-as_tibble(info)

#######TIC graph
ggplot(data=info, aes(y=totIonCurrent, x=retentionTime, group=1)) +
  geom_line(color="darkgreen")+labs(title="Total Ion Current", x="Aquisition time [s]", y = "Absolute Intensity")
#Obtain useful information from metadata
#Number of peaks per spectrum
PeaksNumber<-as.numeric(info2[,1])
#retention times for each spectrum
RT<-as.numeric(info2[,3])
#spectra num id will be used later
spectranum<-as.numeric(info2[,6])
#TIC per spectrum
tic<-as.numeric(info2[,2])
#Base peak per spectrum
bp<-as.numeric(info2[,4])
#Base peak intensity per spectrum
bpint<-as.numeric(info2[,5])
#A data frame with useful metadata 
spectradf<-data.frame(RT,spectranum,PeaksNumber,tic,bp,bpint)
#Convert to tibble
spectradf<-as_tibble(spectradf)
#read samples aquisition times, copy them from an excel or cvs table and run the following line. Format used: start time (min)//end time (min)//start time(s)	//end	time (s)//aquisition time//condition (sample name, will be used to identify in following processing)
#This table should be ordered alphabetically by condition (sample name)
data<-read.table(pipe("pbpaste"),sep="\t",header=T)
#get sample names
samples = data[6]
#Set counter
counter<-1
#make space for the data
processed_data<-data.frame()
#get spectra correspondent to each sample based on the start and end times input in the table
#while the counter is less than the number of samples
while (counter<=nrow(data)) {
  #gets start aquisition time
  rt1<-data[counter,3]
  #gets end aquisition time
  rt2<-data[counter,4]
  #gets sample name
  name <-data[counter,6]
  #Extracts the spectra correspondent to a sample based on start and end times, and adds the sample name as a column, all spectrum from samples are now in the table
  processed_data = rbind(processed_data, data.frame(filter(spectradf, between(RT, rt1,rt2)),name))
  #add to counter to continue the loop with a new sample
  counter=counter+1
}
#converts the table to a tibble
processed_data<-as_tibble(processed_data)
#count spectrum per sample
spectrumpersample<-processed_data%>%count(name)
#Get all spectra number per sample
#We need three counters for this loop to obtain the spectrum number that correspond to each sample. spectra num is a unique id.
i<-1
counter<-1
counter2<-0
while (i<=nrow(spectrumpersample)) {
  sample<-spectrumpersample[i,1]
  sample<-as.character(sample)
  start <-counter
  spectrumnumber<-spectrumpersample[i,2]
  spectrumnumber<-as.numeric(spectrumnumber)
  end<-counter2+spectrumnumber
  assign(sample,processed_data$spectranum[start:end])
  counter2<-end
  counter<-end+1
  i<-i+1
}
counter<-1
spectramatrix<-data.frame()
#############If this loop is not working
#The problem is that they are not ordered in the same way, the metadata information needs to be ordered from A to Z by name!!
#This loop takes a long time to run for larger acquisition times, so you can leave it running at night for samples with more than 30 min acquisitions, or get a snack :), you might use a walk
for (i in 1:nrow(samples)) {
  spectrumnum<-as.numeric(spectrumpersample[i,2])
  sample<-get(samples[i,1])
  counter<-1
  while (counter<=spectrumnum) {
    n<-sample[counter]
    spectramatrix = rbind(spectramatrix, data.frame(mz(raw_data[[n]]),intensity(raw_data[[n]]),spectrumpersample[i,1],n))
    counter<-counter+1
  }
}
###################start here on monday

spectramatrix<- spectramatrix[complete.cases(spectramatrix), ]

spectramatrix<-as_tibble(spectramatrix)

colnames(spectramatrix)<-c("mz","intensity","sample","spectranum")
distinct(spectramatrix, spectranum)
spectramatrix%>% group_by(sample)%>%distinct(spectranum)

spectramatrix<-spectramatrix %>% filter(intensity > 0)
peaksperspectra<-spectramatrix %>%                    # take the data.frame "data"
  filter(!is.na(mz)) %>%    # Using "data", filter out all rows with NAs in aa 
  group_by(spectranum) %>%          # Then, with the filtered data, group it by "bb"
  summarise(peaks = n_distinct(mz))   # Now summarise with unique elements per group

################Raw data
library(plyr)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

####log2 scale data
scaledlog2<-spectramatrix
scaledlog2$intensity<-log2(scaledlog2$intensity)
####Filter spectrum in the adipo mass range
#in vivo new markers
#adipo_log2<-filter(scaledlog2, between(mz, 848.761512, 848.778488))
#adipo_log2<-filter(scaledlog2, between(mz, 876.791232,876.808768))
#adipo_log2<-filter(scaledlog2, between(mz, 786.592134, 786.607866)) 
#1 ppm
adipo_log2<-filter(scaledlog2, between(mz, 786.564617, 786.566183)) 

#adipo_log2<- filter(scaledlog2, between(mz, 804.5719, 804.588046)|between(mz, 830.5816, 830.5983)) 
#adipo_log2<- filter(scaledlog2, between(mz, 804.5719, 804.588046)) 


#adipo_log2 %>%                    # take the data.frame "data"
 # group_by(spectranum) %>%          # Then, with the filtered data, group it by "bb"
  #summarise(Unique_Elements = n_distinct(spectranum))   # Now summarise with unique elements per group


adipo_log2$nomalizationmethod<-"log2"
tiff(filename = "scaledlog2.tiff",
     width = 15, height = 10, units = "cm",res=1500)
ggplot(adipo_log2, aes(sample, intensity)) +
  geom_boxplot(color="darkblue", fill="skyblue", alpha=0.2)+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("scaled to log2") +
  xlab("Sample id")+ylab("log2 scaled intensity ")
dev.off()


########get the spectrum number from the spectrum with a signal for adipo
adipoidx<-adipo_log2$spectranum
#####Get complete spectrum from the spectrum with a adipo signal
adipospectrum<-tibble()
for (i in 1:length(adipoidx)) {
  idx<-adipoidx[i]
  spectraadipoidx<-filter(scaledlog2, spectranum == idx)
  adipospectrum = rbind(adipospectrum,spectraadipoidx)
}
####round intensity to 4 digits

adipospectrum$mz<-round(adipospectrum$mz,digits = 4)

#max_int<- normspectramatrix %>% 
 # dplyr::group_by(spectranum) %>% 
  #dplyr::summarise(int_max = max(intensity)) 
  
#normspectramatrix<-normspectramatrix %>% inner_join(max_int, by="spectranum")

#normspectramatrix$int_norm <- normspectramatrix$intensity/normspectramatrix$int_max

#write file
write_csv(adipospectrum, file="WT.csv")


