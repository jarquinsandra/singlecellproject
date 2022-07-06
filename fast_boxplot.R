# library
library(ggplot2)

# create a data frame

# grouped boxplot
tiff(filename = "X786.6027.tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(PC, aes(y=X786.6027, x=sample)) + 
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = c(0, 1))+
  ylab('intensity')+ theme_classic()
dev.off()

