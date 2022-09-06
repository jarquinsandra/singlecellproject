# library
library(ggplot2)

# create a data frame

# grouped boxplot
#X810.5985 
tiff(filename = "Fig 5k PC361.tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(allcells, aes(y=X810.5985, x=sample)) + 
  geom_boxplot(outlier.shape = NA, color='#CCCCCC', fill = c('#9e0242','#5e4fa2'))+
  coord_cartesian(ylim = c(0, 1.5))+
  ylab('Norm.int.')+ theme_minimal()+ theme(text = element_text(size = 20))
dev.off()

tiff(filename = "Fig 5i PC362.tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(allcells, aes(y=X786.5983, x=sample)) + 
  geom_boxplot(outlier.shape = NA, color='#CCCCCC', fill = c('#9e0242','#5e4fa2'))+
  coord_cartesian(ylim = c(0, 2))+
  ylab('Norm.int.')+ theme_minimal()+ theme(text = element_text(size = 20))
dev.off()

tiff(filename = "Fig 5e valine.tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(allcells, aes(y=X118.0857, x=sample)) + 
  geom_boxplot(outlier.shape = NA, color='#CCCCCC', fill = c('#9e0242','#5e4fa2'))+
  coord_cartesian(ylim = c(0, 2))+
  ylab('Norm.int.')+ theme_minimal()+ theme(text = element_text(size = 20))
dev.off()

tiff(filename = "Fig 5c creatine .tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(allcells, aes(y=X132.0762, x=sample)) + 
  geom_boxplot(outlier.shape = NA, color='#CCCCCC', fill = c('#9e0242','#5e4fa2'))+
  coord_cartesian(ylim = c(0, 1))+
  ylab('Norm.int.')+ theme_minimal()+ theme(text = element_text(size = 20))
dev.off()


tiff(filename = "Fig 5g carnitine.tiff",
     width = 10, height = 10, units = "cm",res=500)
ggplot(allcells, aes(y=X162.1123, x=sample)) + 
  geom_boxplot(outlier.shape = NA, color='#CCCCCC', fill = c('#9e0242','#5e4fa2'))+
  coord_cartesian(ylim = c(0, 0.4))+
  ylab('Norm.int.')+ theme_minimal()+ theme(text = element_text(size = 20))
dev.off()





