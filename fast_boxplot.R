# library
library(ggplot2)

# create a data frame

# grouped boxplot
tiff(filename = "X162.1123_2.tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(allcells, aes(y=X162.1123, x=sample)) + 
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = c(0, 1))+
  ylab('intensity')+ theme_classic()
dev.off()

