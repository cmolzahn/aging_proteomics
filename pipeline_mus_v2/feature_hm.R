library(tidyverse)


dat <- read.delim("c2_pub_features_z.txt", header = T, sep = "\t")
pdat <- read.delim("c2_pub_features_p.txt", header = T, sep = "\t")
## confirm that your features are in the same order in both files, otherwise this will be inaccurate 
dat.gat <- dat %>% gather(key = "cls", value = "zscore", -var)
p.gat <- pdat %>% gather(key = "cls", value = "pvalue", -var)
dat.gat$pval <- p.gat$pvalue

## this orders the x and y axis features based on choices I made manually. You will likely want to change the order.
p1 <- ggplot(dat.gat,aes(factor(cls, level = c("hpa", "cpa", "hpy", "cpy")), factor(var, level = c("dis","rbp","phs","mfc","ssc","csl","chg", "pip","len","ssa","plc","tsl","pfm","gvy","ssb")), fill = zscore))+
  geom_tile()+
  geom_text(aes(label = format(pval, scientific = TRUE)))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_gradient2(low ="slateblue", mid = "white", high = "tomato", midpoint = 0, limits = c(-1, 1), na.value = "dimgray")


pdf("cohort2_pellet_hm.pdf")
p1
dev.off()
