t.test(KC$X118.0857, WT$X118.0857, alternative = "two.sided", var.equal = FALSE)
group1<-filter(clusters,hdb_cluster==-1)

KC<-filter(clusters,hdb_cluster==0)
WT<-filter(clusters,hdb_cluster==1)
mix<-filter(clusters,hdb_cluster==-1)
KC<-KC %>% filter(across(everything(), ~ !str_detect(., "WT")))
WT<-WT %>% filter(across(everything(), ~ !str_detect(., "KC")))
mixKC<-mix %>% filter(across(everything(), ~ !str_detect(., "WT")))
mixWT<-mix %>% filter(across(everything(), ~ !str_detect(., "KC")))
mixlabelsKC<-replicate(212, "KC_3")
mixlabelsWT<-replicate(276, "WT_3")
mixKC$sample<-mixlabelsKC
mixWT$sample<-mixlabelsWT


WTKC<-rbind(KC,WT)
allcells<-rbind(KC,WT,mixKC,mixWT)

comparison1<-rbind(group1,group2)
comparison1<-as_tibble(comparison1)
stats <- data.frame( )
counter<-3
pdf("features_group1and2.pdf")
while (counter<ncol(comparison1)) {
  
  plotdt<-comparison1%>%dplyr::pull(counter)
  
  boxplot(plotdt~comparison1$hdb_cluster,
          main= names(comparison1[counter]),
          xlab="Group",
          ylab="Normalized intensity",
          col="orange",
          border="brown"
  )
  start=
  start= plotdt[1:553]
  end= plotdt[554:length(plotdt)]
  #res <- wilcox.test(start,end,paired=TRUE)
  res<-t.test(start,end, alternative = "two.sided", var.equal = FALSE)
  #res<-kruskal.test(list(start,end))
  print(res)
  
  stats[counter,1] <- res[1]
  stats[counter,2] <- res[2]
  stats[counter,3] <- res[3]
  stats[counter,4]<-names(comparison1[counter])
  counter=counter+1
}
dev.off()
#selects significatively different features with a p value smaller than 0.001, this parameter can be modified
sigfeatures<-subset(stats, p.value<0.05)
sigfeatures<-filter(sigfeatures, p.value!=0)

sigfeaturesmasses<-sigfeatures[4]


comparison2<-rbind(group1,group0)
comparison2<-as_tibble(comparison2)
stats <- data.frame( )
counter<-3
pdf("features_group1and0.pdf")
while (counter<ncol(comparison2)) {
  
  plotdt<-comparison2%>%dplyr::pull(counter)
  
  boxplot(plotdt~comparison2$hdb_cluster,
          main= names(comparison2[counter]),
          xlab="Group",
          ylab="Normalized intensity",
          col="orange",
          border="brown"
  )
  start=
    start= plotdt[1:553]
  end= plotdt[554:length(plotdt)]
  #res <- wilcox.test(start,end,paired=TRUE)
  res<-t.test(start,end, alternative = "two.sided", var.equal = FALSE)
  #res<-kruskal.test(list(start,end))
  print(res)
  
  stats[counter,1] <- res[1]
  stats[counter,2] <- res[2]
  stats[counter,3] <- res[3]
  stats[counter,4]<-names(comparison1[counter])
  counter=counter+1
}
dev.off()
#selects significatively different features with a p value smaller than 0.001, this parameter can be modified
sigfeatures2<-subset(stats, p.value<0.05)
sigfeatures2<-filter(sigfeatures2, p.value!=0)

sigfeaturesmasses2<-sigfeatures2[4]

write_csv(sigfeatures2, file = "sigfeaturesmassescomparison2.csv")



comparison3<-rbind(group0,group2)
comparison2<-as_tibble(comparison3)
stats <- data.frame( )
counter<-3
pdf("features_group2and0.pdf")
while (counter<ncol(comparison3)) {
  
  plotdt<-comparison3%>%dplyr::pull(counter)
  
  boxplot(plotdt~comparison3$hdb_cluster,
          main= names(comparison3[counter]),
          xlab="Group",
          ylab="Normalized intensity",
          col="orange",
          border="brown"
  )
  start=
    start= plotdt[1:3207]
  end= plotdt[3208:length(plotdt)]
  #res <- wilcox.test(start,end,paired=TRUE)
  res<-t.test(start,end, alternative = "two.sided", var.equal = FALSE)
  #res<-kruskal.test(list(start,end))
  print(res)
  
  stats[counter,1] <- res[1]
  stats[counter,2] <- res[2]
  stats[counter,3] <- res[3]
  stats[counter,4]<-names(comparison1[counter])
  counter=counter+1
}
dev.off()
#selects significatively different features with a p value smaller than 0.001, this parameter can be modified
sigfeatures3<-filter(stats, p.value<0.05)
sigfeatures3<-filter(sigfeatures3, p.value!=0)
sigfeaturesmasses3<-sigfeatures3[4]

write_csv(sigfeatures3, file = "sigfeaturesmassescomparison3.csv")


sigfeaturesmasses<-head(sigfeaturesmasses, -1)

sigfeaturesmasses2<-head(sigfeaturesmasses2, -1)
sigfeaturesmasses3<-head(sigfeaturesmasses3, -1)

venn.diagram(
  x = list(sigfeaturesmasses, sigfeaturesmasses2, sigfeaturesmasses3),
  category.names = c("comparison 1" , "comparison 2 " , "comparison 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)





