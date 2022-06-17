library(statsExpressions)
library(ggplot2)
library(ggstatsplot)

dat<-clusters[-1]
dat<-dat[-1]
dat<-dat[-160]

# Edit from here
x <- which(names(dat) == "hdb_cluster") # name of grouping variable
y <- which(names(dat) == "X415.2100" # names of variables to test
           | names(dat) == "X416.2260" |
             names(dat) == "X419.3135" |
             names(dat) == "X427.2494"|
             names(dat) == "X431.2453"|
             names(dat) == "X441.2654"|
             names(dat) == "X442.2716"|
             names(dat) == "X448.2304"|
             names(dat) == "X452.2230"|
             names(dat) == "X460.2525"|
             names(dat) == "X474.2671"|
             names(dat) == "X476.3562"|
             names(dat) == "X481.3111"|
             names(dat) == "X490.2623"|
             names(dat) == "X498.8019"|
             names(dat) == "X518.3656"|
             names(dat) == "X540.4229"|
             names(dat) == "X585.3351"|
             names(dat) == "X599.3884"|
             names(dat) == "X663.4515"|
             names(dat) == "X680.4765"|
             names(dat) == "X681.4799"|
             names(dat) == "X685.4313"|
             names(dat) == "X701.4063"|
             names(dat) == "X713.4182"|
             names(dat) == "X734.5664"|
             names(dat) == "X758.5656"|
             names(dat) == "X782.5656"|
             names(dat) == "X786.5983"|
             names(dat) == "X810.5985")

method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
#my_comparisons <- list(c(0, 1), c(1, 2), c(2, 0)) # comparisons for post-hoc tests
# Edit until here


# Edit at your own risk
for (i in y) {
  for (j in x) {
    p <- ggboxplot(dat,
                   x = colnames(dat[j]), y = colnames(dat[i]),
                   color = colnames(dat[j]),
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(dat[, i], na.rm = TRUE)
      )
      #+ stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}



# Comparison between species

# edit from here
x <- "hdb_cluster"
cols <- 3:159 # the 4 continuous dependent variables
type <- "parametric" # given the large number of observations, we use the parametric version
paired <- FALSE # FALSE for independent samples, TRUE for paired samples
# edit until here

# edit at your own risk
plotlist <-
  purrr::pmap(
    .l = list(
      data = list(as_tibble(dat)),
      x = x,
      y = as.list(colnames(dat)[cols]),
      plot.type = "box", # for boxplot
      type = type, # parametric or nonparametric
      pairwise.comparisons = TRUE, # to run post-hoc tests if more than 2 groups
      pairwise.display = "significant", # show only significant differences
      bf.message = FALSE, # remove message about Bayes Factor
      centrality.plotting = FALSE # remove central measure
    ),
    .f = ifelse(paired, # automatically use ggwithinstats if paired samples, ggbetweenstats otherwise
                ggstatsplot::ggwithinstats,
                ggstatsplot::ggbetweenstats
    )
  )

# print all plots together with statistical results
pdf('comparativeplots.pdf',width = 14)

for (i in 1:length(plotlist)) {
  print(plotlist[[i]] +
          labs(caption = NULL)) # remove caption
}
dev.off()

ggbetweenstats(testdata,"X441.7854") %>%
  extract_stats()