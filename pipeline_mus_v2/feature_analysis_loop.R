library(tidyverse)



df1 <- read.delim("data1.txt", header =  TRUE, sep = "\t")
df2 <- read.delim("data2.txt", header = TRUE, sep = "\t")
dat <- read.table("mouse_dat_pfam_corrected.txt", header = T)%>%
  dplyr::select(nam,dis, len, run, max, chg, net, plc, abd_log, gvy, pip, iac,mfc,csl,tsl,ssa,ssb,ssc,pol,tmx,toh,sst, phs, pfm, mag,rbp, pap_cnt,plc_cnt,rbp_cnt)



df1$cls <- "cond1"
df2$cls <- "cond2"


df <- rbind(df1, df2)

dat.out <- c()
for (i in unique(df$cls)) {
  df.fil <- df %>% filter(cls == i)
  dat.tmp <-dat[which(dat$nam %in% df.fil$Entry),]
  dat.tmp$cls <- paste(i)
  dat.out <- rbind(dat.out, dat.tmp)
}
dat$cls <- "pme"
dat.tmp <- rbind(dat.out, dat)
df.gat <- dat.tmp %>% gather(key = var, value= value,-cls, -nam)
p.out <-c()
n.out <- c()
z.out <- c()
plot_list <- list()
for (i in unique(df.gat$var)) {
  dat.fil <- filter(df.gat, var == i)
  if (i == "pap_cnt"|i=="prw_cnt"|i== "plc_cnt"| i=="rbp_cnt") {
    pdf(paste(i,".pdf"))
    barplot(c(length(dat.fil[which(dat.fil$cls=="cond1" & dat.fil$value==1),]$nam) / length(dat.fil[which(dat.fil$cls=="cond1"),]$nam),
              length(dat.fil[which(dat.fil$cls=="cond2" & dat.fil$value==1),]$nam) / length(dat.fil[which(dat.fil$cls=="cond2"),]$nam),
              length(dat.fil[which(dat.fil$cls=="pme" & dat.fil$value==1),]$nam)/
                length(dat.fil[which(dat.fil$cls=="pme"),]$nam)),
            ylab=paste(i),names=c("cond1","cond2", "pme"),
            col=c("slateblue","tomato","gainsboro"))
    dev.off()
    
    p <- c(i,
           cond1 = 
             fisher.test(matrix(c(length(dat.fil[which(dat.fil$cls=="cond1" & dat.fil$value==1),]$nam),(length(dat.fil[which(dat.fil$cls=="pme" & dat.fil$value==1),]$nam) - length(dat.fil[which(dat.fil$cls=="cond1" & dat.fil$value==1),]$nam)),length(dat.fil[which(dat.fil$cls=="cond1" & dat.fil$value!=1),]$nam),(length(dat.fil[which(dat.fil$cls=="pme" & dat.fil$value!=1),]$nam) - length(dat.fil[which(dat.fil$cls=="cond1" & dat.fil$value!=1),]$nam))),2,2),alternative='two.sided')$p.value,
           cond2 =
             fisher.test(matrix(c(length(dat.fil[which(dat.fil$cls=="cond2" & dat.fil$value==1),]$nam),(length(dat.fil[which(dat.fil$cls=="pme" & dat.fil$value==1),]$nam) - length(dat.fil[which(dat.fil$cls=="cond2" & dat.fil$value==1),]$nam)),length(dat.fil[which(dat.fil$cls=="cond2" & dat.fil$value!=1),]$nam),(length(dat.fil[which(dat.fil$cls=="pme" & dat.fil$value!=1),]$nam) - length(dat.fil[which(dat.fil$cls=="cond2" & dat.fil$value!=1),]$nam))),2,2),alternative='two.sided')$p.value
    )
    
    
  }else{
    p <- c(i,
           
           cond1 = wilcox.test(dat.fil[which(dat.fil$cls == "cond1"),]$value, dat.fil[which(dat.fil$cls == "pme"),]$value)$p.val,
           cond2 = wilcox.test(dat.fil[which(dat.fil$cls == "cond2"),]$value, dat.fil[which(dat.fil$cls == "pme"),]$value)$p.val
    )
    
    n <- c(var = i,
           
           cond1=length(which(!is.na(dat.fil[which(dat.fil$cls == "cond1"),]$value))),
           cond2 = length(which(!is.na(dat.fil[which(dat.fil$cls == "cond2"),]$value)))
    )
    z <- c(var = i, cond1= (mean(na.omit(dat.fil[which(dat.fil$cls == "cond1"),]$value)) - mean(na.omit(dat.fil[which(dat.fil$cls == "pme"),]$value)))/ sd(na.omit(dat.fil[which(dat.fil$cls == "cond1"),]$value)), 
           cond2= (mean(na.omit(dat.fil[which(dat.fil$cls == "cond2"),]$value)) - mean(na.omit(dat.fil[which(dat.fil$cls == "pme"),]$value)))/ sd(na.omit(dat.fil[which(dat.fil$cls == "cond2"),]$value))
    )
    
    z.out <- as.data.frame(rbind(z.out, z))
    n.out <-as.data.frame(rbind(n.out, n))
    p1 <- ggplot(dat.fil, aes(x=factor(cls, level=c( "cond1", "cond2","pme")),y=value)) + 
      geom_violin(aes(fill=cls),trim=FALSE,adjust=1.5) + 
      labs(x="Set", y=paste(i)) + 
      geom_boxplot(width=0.2,fill="light gray",outlier.shape=NA) + 
      theme_classic(base_size=22) + 
      theme(axis.text.x=element_text(angle=90,hjust=1)) + 
      ylim(c(min(dat.fil$value, na.rm = TRUE),max(dat.fil$value,na.rm = TRUE))) + 
      scale_fill_manual(values=c("tomato","slateblue","gainsboro")) + 
      theme(legend.position="none")
    
  }
  p.out <-as.data.frame(rbind(p.out, p))
  plot_list[[i]] = p1
}



p.table <-data.frame(var = p.out$V1, 
                     
                     cond1 = p.adjust(p.out$cond1, method = "fdr"),
                     cond2 = p.adjust(p.out$cond2, method = "fdr"))

write.table(p.table, "features_p-value.txt", row.names = FALSE, sep = "\t")
write.table(n.out, "features_n.txt", row.names = FALSE, sep = "\t")
pdf("plot_data_pel.pdf")
plot_list
dev.off()
