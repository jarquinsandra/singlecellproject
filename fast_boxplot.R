# library
library(ggplot2)

# create a data frame

# grouped boxplot
tiff(filename = "X302.3044.tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(WTKC, aes(y=X302.3044, x=sample)) + 
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = c(0, 0.3))+
  ylab('intensity')+ theme_classic()
dev.off()

