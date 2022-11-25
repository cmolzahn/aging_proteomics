library(FactoMineR)
library(factoextra)
library(corrplot)
library(tidyr)
library(cowplot)
library(stringr)

dirs <- dir("pipeline_mus/src")
dirs <- dirs[-c(10,14)]
dflist <- list()
for (i in dirs) {
  df <- read.delim(file.path("pipeline_mus/src", i), header = FALSE, sep = "\t")
  dflist[[i]] <- df
}

nam <- as.data.frame(dirs)
colnames(nam) <- "nam"

# Protein feature data
dat <- read.table("mouse_dat.txt", header=T)

out <- c()
##############this loop filters dat for TRUE values in each of the data frames and outputs the mean value for each feature for each data frame
for (i in dflist) {
  dat$var <-  dat$nam %in% i[,1]
  dat.fil <- dat%>%
    filter(dat$var ==TRUE)
  mn <- colMeans(dat.fil[2:ncol(dat.fil)], na.rm = TRUE)
  out <- rbind(out, mn)
  print(out)
}
out<- cbind(out, nam)

#############feature data frame

dat.out <- out
my.sub = c()
my.sub$nam <- dat.out$nam
my.sub$dis <- dat.out$dis
my.sub$len <- dat.out$len
# my.sub$run <- dat.out$run
#  my.sub$max <- dat.out$max
my.sub$chg <- dat.out$chg
# my.sub$net <- dat.out$net 
# my.sub$abd <- dat.out$abd
my.sub$gvy   <- dat.out$gvy
my.sub$pip <- dat.out$pip
# my.sub$tgo <- dat.out$tgo
# my.sub$iac <- dat.out$iac
my.sub$mfc <- dat.out$mfc
my.sub$csl <- dat.out$csl
# my.sub$sol <- dat.out$sol
# my.sub$sto <- dat.out$sto
# my.sub$stc <- dat.out$stc
my.sub$ssa <- dat.out$ssa
my.sub$ssb <- dat.out$ssb
my.sub$ssc <- dat.out$ssc
# my.sub$sft <- dat.out$sft
# my.sub$scn <- dat.out$scn
# my.sub$sbb <- dat.out$sbb
# my.sub$pol   <- dat.out$pol
# my.sub$pap <- dat.out$pap
my.sub$plc <- dat.out$plc
# my.sub$rna <- dat.out$rna
my.sub$rbp <- dat.out$rbp
# my.sub$tmx <- dat.out$tmx
# my.sub$tom <- dat.out$tom
# my.sub$toh <- dat.out$toh
# my.sub$yrd <- dat.out$yrd
my.sub$phs   <- dat.out$phs
# my.sub$mag <- dat.out$mag
# my.sub$oxi <- dat.out$oxi
# my.sub$prw <- dat.out$prW
# my.sub$coc <- dat.out$coc
my.sub$ppf <- dat.out$ppf
# my.sub$pfm <- dat.out$pfm

my.sub <- as.data.frame(my.sub)
my.wrk <- my.sub[,-1]
rownames(my.wrk) <- my.sub$nam
nams <- unique(my.sub$nam)
my.pca <- PCA(my.wrk, graph=F)
my.var <- get_pca_var(my.pca)
regexp <- "[[:digit:]]+"

p2 <-fviz_pca_var(my.pca, col.var="cos2", gradient.col=c("#C2C2C2", "#BB2200"),repel=F)
p4 <- fviz_pca_ind(my.pca, col.ind = as.numeric(str_extract(dat.out$nam, regexp)), gradient.col=c("slateblue","#C2C2C2","tomato"))+
  scale_color_gradient2(low = "slateblue", mid = "#C2C2C2", high = "tomato", midpoint = 40 )
p1 <-corrplot(my.var$cos2,is.corr=F)


plot_grid(p2,p4)