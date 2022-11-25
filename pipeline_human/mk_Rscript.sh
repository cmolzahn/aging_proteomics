#! /bin/bash

observable="dis len run max chg net abd_log gvy pip iac mfc csl tsl ssa ssb ssc pol tmx toh sst phs mag rbp"
count="pap_cnt plc_cnt rbp_cnt"
aas="R H K D E S T N Q C U G P A V I L M F Y W"
hmn_exp="drummond_AD hales_AD seyfried_FTD cherry_CTE"
compare_to="human"
col_list="\"firebrick\",\"firebrick\",\"firebrick\",\"firebrick\",\"gainsboro\""

plot_parms () {
  if [ $1 == "dis" ]; then echo "Disorder 0 100"; fi
  if [ $1 == "len" ]; then echo "Seq_Length 0 2000"; fi
  if [ $1 == "run" ]; then echo "n_IDRs 0 10"; fi
  if [ $1 == "max" ]; then echo "Max_IDR 0 600"; fi
  if [ $1 == "chg" ]; then echo "%Charged_AA 0 60"; fi
  if [ $1 == "net" ]; then echo "Net_Charge -50 50"; fi
  if [ $1 == "abd_log" ]; then echo "Log_Abundance -4 5"; fi
  if [ $1 == "gvy" ]; then echo "GRAVY -3 2"; fi
  if [ $1 == "pip" ]; then echo "PScore 0 8"; fi
  if [ $1 == "iac" ]; then echo "Biogrid_Interactions 0 35"; fi
  if [ $1 == "mfc" ]; then echo "n_MoRFs 0 20"; fi
  if [ $1 == "csl" ]; then echo "Camsol -10 10"; fi
  if [ $1 == "sto" ]; then echo "TANGO_Sections 0 50"; fi
  if [ $1 == "tsl" ]; then echo "%hi_TANGO 0.0 0.75"; fi
  if [ $1 == "ssa" ]; then echo "%Alpha 0 90"; fi
  if [ $1 == "ssb" ]; then echo "%Beta 0 90"; fi
  if [ $1 == "ssc" ]; then echo "%Coil 0 100"; fi
  if [ $1 == "scn" ]; then echo "Chain_softness 0.9 1.5"; fi
  if [ $1 == "sbb" ]; then echo "BB_softness 0.5 1.2"; fi
  if [ $1 == "pol" ]; then echo "Polarizabilty 60 80"; fi
  if [ $1 == "yrd" ]; then echo "YtoR_separation 0.0 0.75"; fi
  if [ $1 == "tmx" ]; then echo "Tm 30 70"; fi
  if [ $1 == "toh" ]; then echo "Turnover 0 100"; fi
  if [ $1 == "sst" ]; then echo "SuperSaturation -3 5"; fi
  if [ $1 == "phs" ]; then echo "Annotated_PSites 0 30"; fi
  if [ $1 == "mag" ]; then echo "MaGS -3 5"; fi
  if [ $1 == "oxi" ]; then echo "OxidationSites 0 15"; fi 
  if [ $1 == "rbp" ]; then echo "RBPPred 0 1"; fi
}

############################################################
## SCRIPT START                                    ^_~_^
## It's dangerous to go alone! Take this with you (='.'=)
#############################################################


echo "## LOAD LIBRARIES ##
library(plotrix)
library(grid)
library(gridExtra)
library(scales)
library(ggplot2)

## UNIVERSAL VARIABLES ##
mylevhmn=c("
for exp in $hmn_exp; do
  echo -n "\"$exp\","
done 
echo "\"human\")"

echo "## READ DATA ##
hmn_all = read.table(\"dat/human_dat.txt\", header=T)
hmn_aas = read.table(\"dat/human_pme.aa\", header=T)
hmn_pme = read.table(\"dat/human_pme.txt\")"

for exp in $hmn_exp; do
  echo "$exp = read.table(\"src/$exp.txt\")"
done

############
# CALC
############
echo "
hmn_all\$abd_log <- log(hmn_all\$abd,10)
hmn_all\$rbp_cnt <- lapply(hmn_all\$rbp, round)
hmn_all\$pap_cnt <- ifelse(hmn_all\$pap>=0.05,1,0)
hmn_all\$plc_cnt <- ifelse(hmn_all\$plc>0.0,1,0)
hmn_all\$abd_log_cen <- hmn_all\$abd_log-mean(na.omit(hmn_all\$abd_log))
hmn_all\$abd_log_nrm <- hmn_all\$abd_log_cen/sd(na.omit(hmn_all\$abd_log_cen))
hmn_all\$csl_cen <- hmn_all\$csl-mean(na.omit(hmn_all\$csl))
hmn_all\$csl_nrm <- hmn_all\$csl_cen/sd(na.omit(hmn_all\$csl_cen))
hmn_all\$tsl     <- hmn_all\$stc/hmn_all\$len
hmn_all\$sst     <- hmn_all\$abd_log - hmn_all\$csl_nrm

my.hmn      <- hmn_all[which(hmn_all\$nam %in% hmn_pme\$V1),]
my.hmn\$tag <- \"human\"
my.aam      <- hmn_aas[which(hmn_aas\$nam %in% hmn_pme \$V1),]
my.aam\$tag <- \"human\"
ref.hmn     <- my.hmn
ref.aam     <- my.aam

#######################
# DATA BLOCK - my god forgive my soul
#######################
"

for exp in $hmn_exp; do 
 echo "#
 my.tmp      <- ref.hmn[which(ref.hmn\$nam %in% $exp\$V1),]
 my.tmp\$tag <- \"$exp\"
 my.hmn      <- rbind(my.hmn,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam\$nam %in% $exp\$V1),]
 my.tmp\$tag <- \"$exp\"
 my.aam      <- rbind(my.aam,my.tmp)"
done

echo "#
###########
# PLOT   
###########
pdf(\"plot_data.pdf\")

my.hmn\$tag <-factor(my.hmn\$tag, levels=mylevhmn)
my.aam\$tag <-factor(my.aam\$tag, levels=mylevhmn)"

for obs in $observable; do 
tmp=$(plot_parms $obs)
cnt=0; ARGS=()
for el in $tmp; do ARGS[$cnt]=$el; ((cnt=cnt+1)); done
echo "## $obs ## 
ggplot(my.hmn, aes(x=tag,y=$obs)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x=\"Set\", y=\"${ARGS[0]}\") + geom_boxplot(width=0.2,fill=\"light gray\",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(${ARGS[1]},${ARGS[2]}) + scale_fill_manual(values=c($col_list)) + theme(legend.position=\"none\")
#"
done

for obs in $count; do
  echo "## $obs ##"
  echo -n "barplot(c("
  var=''
  for exp in $hmn_exp; do
    var="${var}length(my.hmn[which(my.hmn\$tag==\"$exp\" & my.hmn\$$obs==1),]\$nam) / length(my.hmn[which(my.hmn\$tag==\"$exp\"),]\$nam),"
  done
  var="${var}length(my.hmn[which(my.hmn\$tag==\"$compare_to\" & my.hmn\$$obs==1),]\$nam) / length(my.hmn[which(my.hmn\$tag==\"$compare_to\"),]\$nam)),ylab=\"$obs Fraction\",names=c("
  for exp in $hmn_exp; do var="${var}\"$exp\","; done
  var="${var}\"human\"),col=c($col_list))"
  echo "$var"
done

echo "dev.off()"

echo "#
###############
# AA Plot
###############
pdf(\"plot_aa.pdf\")"

for obs in $aas; do 
echo "## $obs ##
ggplot(my.aam, aes(x=tag,y=$obs)) + geom_violin(trim=FALSE, fill=c(\"firebrick2\"),adjust=1.5) + labs(x=\"Set\", y=\"$obs\") + geom_boxplot(width=0.2,fill=\"light gray\",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0,20)
#"
done

echo "dev.off()"

echo "#
###########
# P-Values
###########
#"
for obs in $observable; do
  echo -n "p_$obs <- c("
  var=''
  for exp in $hmn_exp; do
    var="${var}wilcox.test(my.hmn[which(my.hmn\$tag=='${exp}'),]\$$obs, my.hmn[which(my.hmn\$tag=='$compare_to'),]\$$obs)\$p.val,"
  done
  var="${var})"; var="${var%,)})"
  echo "$var"
done

for obs in $count; do
  echo -n "p_$obs <- c("
  var=''
  for exp in $hmn_exp; do
    hit_sample="length(my.hmn[which(my.hmn\$tag==\"$exp\" & my.hmn\$$obs==1),]\$nam)"
    hit_left="(length(my.hmn[which(my.hmn\$tag==\"human\" & my.hmn\$$obs==1),]\$nam) - length(my.hmn[which(my.hmn\$tag==\"$exp\" & my.hmn\$$obs==1),]\$nam))"
    miss_sample="length(my.hmn[which(my.hmn\$tag==\"$exp\" & my.hmn\$$obs!=1),]\$nam)"
    miss_left="(length(my.hmn[which(my.hmn\$tag==\"human\" & my.hmn\$$obs!=1),]\$nam) - length(my.hmn[which(my.hmn\$tag==\"$exp\" & my.hmn\$$obs!=1),]\$nam))"
   my_test="fisher.test(matrix(c($hit_sample,$hit_left,$miss_sample,$miss_left),2,2),alternative='two.sided')\$p.value"
   var="${var}$my_test,"
  done
  var="${var})"; var="${var%,)})"
  echo "$var"
done

echo -n "p_hmn <- rbind("
var=''
for obs in $observable; do 
  var="${var}p_${obs},"
done
for obs in $count; do
  var="${var}p_${obs},"
done
var="${var})"
var="${var%,)})"
echo "$var"


echo "#"
cnt=1
for exp in $hmn_exp; do 
  echo "$exp <- p.adjust(p_hmn[,$cnt],method=\"hochberg\")"
  cnt=$((cnt+1))
done
echo -n "corr_pval<-cbind("
var=''
for exp in $hmn_exp; do 
 var="${var}$exp,"
done 
var="${var})"
var="${var%,)})"
echo "$var"
echo " 
ptab<-signif(corr_pval,digits=2)
tab_top    <- head(ptab,20)
tab_bottom <- ptab[-(1:20),]

pdf(\"plot_pval.pdf\")
grid.table(tab_top)
grid.newpage()
grid.table(tab_bottom)
dev.off()
#"

