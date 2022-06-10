dat<-features6cells_pos[-119]
dat<-dat[-1]

# Edit from here
x <- which(names(dat) == "sample") # name of grouping variable
y <- which(names(dat) == "X413.267" # names of variables to test
           | names(dat) == "X414.271" |
             names(dat) == "X415.212" |
             names(dat) == "X419.316"|
             names(dat) == "X432.24"|
             names(dat) == "X434.782"|
             names(dat) == "X441.154"|
             names(dat) == "X453.788"|
             names(dat) == "X471.308"|
             names(dat) == "X481.82"|
             names(dat) == "X485.362")

     
       "X485.814"
"X499.189" "X590.428" "X600.285" "X634.456" "X663.457" "X678.479" "X680.485"
] "X681.485" "X685.439" "X703.58"  "X722.509" "X732.557" "X733.561" "X734.575"
[29] "X746.576" "X746.607" "X758.574" "X772.628" "X782.568" "X784.589" "X785.594"
[36] "X786.601" "X787.61"  "X788.621" "X806.576" "X811.608" "X813.69"  "X834.603"
[43] "X414.297" "X416.217" "X420.32"  "X429.242" "X433.242" "X436.31"  "X441.788"
[50] "X447.349" "X448.236" "X453.168" "X458.322" "X458.348" "X462.148" "X480.338"
[57] "X480.839" "X483.819" "X486.815" "X487.818" "X492.325" "X502.35"  "X502.376"
[64] "X515.843" "X517.84"  "X518.843" "X521.774" "X526.363" "X540.429" "X546.403"
[71] "X635.459" "X664.459" "X665.462" "X682.487" "X686.443" "X706.542" "X723.515"
[78] "X735.576" "X744.593" "X756.558" "X759.577" "X760.59"  "X761.589" "X766.533"
[85] "X780.555" "X789.626" "X794.609" "X808.587" "X810.561" "X810.605" "X811.564"
[92] "X898.611" "X437.196" "X483.984" "X485.815" "X502.351" "X502.853" "X518.37" 
[99] "X854.585" "X901.558" "X430.244" "X524.362" "X679.488" "X720.591" "X767.536"
[106] "X458.827" "X762.598" "X436.813" "X701.414" "X707.543" "X414.272" "X718.58" 
[113] "X583.379" "X680.486" "X430.245" "X682.488" "X433.243"
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("AsPC1_pos", "HeLa_pos"), c("BxPC-3-pos", "PANC-1_pos"), c("HeLa_pos", "cfpan-1-pos")) # comparisons for post-hoc tests
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
x <- "sample"
cols <- 3:117 # the 4 continuous dependent variables
type <- "nonparametric" # given the large number of observations, we use the parametric version
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
pdf('Examplenonp.pdf',width = 14)

for (i in 1:length(plotlist)) {
  print(plotlist[[i]] +
          labs(caption = NULL)) # remove caption
}
dev.off()