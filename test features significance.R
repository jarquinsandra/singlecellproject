library(tidyverse)
library(rstatix)
library(ggpubr)
library(viridis)
#select each feature by group

testdata<-clusters[-1]
testdata<-testdata[-1]
testdata<-testdata[-160]
stats <- data.frame( )
counter<-1
#pdf("features_normalizedscaled.pdf")
while (counter<ncol(testdata)-1) {
  
  plotdt<-testdata%>%dplyr::pull(counter)
  labels<-testdata$hdb_cluster
  testvalues<-tibble(plotdt,labels)
  testdata %>% sample_n_by(hdb_cluster, size = 1)
  feature<-colnames(testdata[counter])
    pwc <- testvalues %>%
    pairwise_t_test(plotdt ~ labels, p.adjust.method = "bonferroni")
  print(pwc)
  testvalues$labels<-as.factor(testvalues$labels)
  # p1<-testvalues %>%
  #   ggplot( aes(y=plotdt, x=labels)) +
  #   geom_boxplot() +
  #   scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #   geom_jitter(color="black", size=0.4, alpha=0.9) +
  #   theme_ipsum() +
  #   theme(
  #     legend.position="none",
  #     plot.title = element_text(size=11)
  #   ) +
  #   ggtitle(feature) +
  #   xlab("Group") +
  #   ylab("Intensity")
  testvalues2<- filter(testvalues, plotdt!=0)
  testvalues2$plotdt<-scale(testvalues2$plotdt)
  p2<-testvalues2 %>%
    ggplot( aes(y=plotdt, x=labels)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    #theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(feature) +
    xlab("Group") +
    ylab("Intensity")
    counter=counter+1
    #print(p1)
    print(p2)
}
#dev.off()

#res <- wilcox.test(start,end,paired=TRUE)
res<-kruskal.test(list(start,end))
print(res)

stats[counter,1] <- res[1]
stats[counter,2] <- res[2]
stats[counter,3] <- res[3]
stats[counter,4]<-names(normout_featureMatrix[counter])

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
