#! /bin/bash

observable="dis len run max chg net abd_log gvy pip iac mfc csl tsl ssa ssb ssc pol tmx toh sst phs mag oxi rbp pfm ppf coc"
count="pap_cnt plc_cnt rbp_cnt prW_cnt"
aas="R H K D E S T N Q C U G P A V I L M F Y W"
mus_exp="WT_11 HTT_12 AD_15 ALS_25 FTD_30 HTT_40 FTD_41 AD_52 AD_87 WT_100 AD_104 WT_104"
col_list="\"lightblue\", \"lightblue\", \"lightblue\", \"dimgray\", \"dimgray\", \"dimgray\", \"dimgray\", \"dimgray\", \"firebrick\", \"firebrick\", \"firebrick\", \"firebrick\", \"gainsboro\""

plot_parms () {
  if [ $1 == "dis" ]; then echo "Disorder 0 100"; fi
  if [ $1 == "len" ]; then echo "Seq_Length 0 2500"; fi
  if [ $1 == "run" ]; then echo "n_IDRs 0 20"; fi
  if [ $1 == "max" ]; then echo "Max_IDR 0 1000"; fi
  if [ $1 == "chg" ]; then echo "%Charged_AA 0 60"; fi
  if [ $1 == "net" ]; then echo "Net_Charge -50 50"; fi
  if [ $1 == "abd_log" ]; then echo "Log_Abundance -4 5"; fi
  if [ $1 == "gvy" ]; then echo "GRAVY -3 1"; fi
  if [ $1 == "pip" ]; then echo "PScore 0 10"; fi
  if [ $1 == "iac" ]; then echo "Binary_Interactions 0 25"; fi
  if [ $1 == "mfc" ]; then echo "n_MoRFs 0 20"; fi
  if [ $1 == "csl" ]; then echo "Camsol -10 10"; fi
  if [ $1 == "sto" ]; then echo "TANGO_Sections 0 50"; fi
  if [ $1 == "tsl" ]; then echo "%hi_TANGO 0.0 0.5"; fi
  if [ $1 == "ssa" ]; then echo "%Alpha 0 90"; fi
  if [ $1 == "ssb" ]; then echo "%Beta 0 90"; fi
  if [ $1 == "ssc" ]; then echo "%Coil 0 100"; fi
  if [ $1 == "scn" ]; then echo "Chain_softness 0.9 1.5"; fi
  if [ $1 == "sbb" ]; then echo "BB_softness 0.5 1.2"; fi
  if [ $1 == "pol" ]; then echo "Polarizabilty 60 80"; fi
  if [ $1 == "yrd" ]; then echo "YtoR_separation 0.0 0.75"; fi
  if [ $1 == "tmx" ]; then echo "Tm 30 70"; fi
  if [ $1 == "toh" ]; then echo "Turnover 0 50"; fi
  if [ $1 == "sst" ]; then echo "SuperSaturation -3 5"; fi
  if [ $1 == "phs" ]; then echo "Annotated_PSites 0 25"; fi
  if [ $1 == "mag" ]; then echo "MaGS -3 5"; fi
  if [ $1 == "oxi" ]; then echo "OxidationSites 0 15"; fi 
  if [ $1 == "rbp" ]; then echo "RBPPred 0 1"; fi
  if [ $1 == "pfm" ]; then echo "N_PFAM 0 20"; fi
  if [ $1 == "ppf" ]; then echo "%PFAM 0 100"; fi
  if [ $1 == "coc" ]; then echo "CorumComplexes 0 10"; fi
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
mylevmus=c("
for exp in $mus_exp; do
  echo -n "\"$exp\","
done 
echo "\"mouse\")"

echo "## READ DATA ##
mus_all = read.table(\"dat/mouse_dat.txt\", header=T)
mus_aas = read.table(\"dat/mouse_pme.aa\", header=T)
mus_pme = read.table(\"dat/mouse_pme.txt\")"

for exp in $mus_exp; do
  echo "$exp = read.table(\"src/$exp.txt\")"
done

############
# CALC
############
echo "
mus_all\$abd_log <- log(mus_all\$abd,10)
mus_all\$rbp_cnt <- lapply(mus_all\$rbp, round)
mus_all\$pap_cnt <- ifelse(mus_all\$pap>=0.05,1,0)
mus_all\$prW_cnt <- ifelse(mus_all\$prW>0.0,1,0)
mus_all\$plc_cnt <- ifelse(mus_all\$plc>0.0,1,0)
mus_all\$abd_log_cen <- mus_all\$abd_log-mean(na.omit(mus_all\$abd_log))
mus_all\$abd_log_nrm <- mus_all\$abd_log_cen/sd(na.omit(mus_all\$abd_log_cen))
mus_all\$csl_cen <- mus_all\$csl-mean(na.omit(mus_all\$csl))
mus_all\$csl_nrm <- mus_all\$csl_cen/sd(na.omit(mus_all\$csl_cen))
mus_all\$tsl     <- mus_all\$stc/mus_all\$len
mus_all\$sst     <- mus_all\$abd_log - mus_all\$csl_nrm

my.mus      <- mus_all[which(mus_all\$nam %in% mus_pme\$V1),]
my.mus\$tag <- \"mouse\"
my.aam      <- mus_aas[which(mus_aas\$nam %in% mus_pme \$V1),]
my.aam\$tag <- \"mouse\"
ref.mus     <- my.mus
ref.aam     <- my.aam

#######################
# DATA BLOCK - my god forgive my soul
#######################
"

for exp in $mus_exp; do 
 echo "#
 my.tmp      <- ref.mus[which(ref.mus\$nam %in% $exp\$V1),]
 my.tmp\$tag <- \"$exp\"
 my.mus      <- rbind(my.mus,my.tmp)
 my.tmp      <- ref.aam[which(ref.aam\$nam %in% $exp\$V1),]
 my.tmp\$tag <- \"$exp\"
 my.aam      <- rbind(my.aam,my.tmp)"
done

echo "#
###########
# PLOT   
###########
pdf(\"age_data.pdf\")

my.mus\$tag <-factor(my.mus\$tag, levels=mylevmus)
my.aam\$tag <-factor(my.aam\$tag, levels=mylevmus)"

for obs in $observable; do 
tmp=$(plot_parms $obs)
cnt=0; ARGS=()
for el in $tmp; do ARGS[$cnt]=$el; ((cnt=cnt+1)); done
echo "## $obs ## 
ggplot(my.mus, aes(x=tag,y=$obs)) + geom_violin(aes(fill=tag),trim=FALSE,adjust=1.5) + labs(x=\"Set\", y=\"${ARGS[0]}\") + geom_boxplot(width=0.2,fill=\"light gray\",outlier.shape=NA) + theme_classic(base_size=22) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(${ARGS[1]},${ARGS[2]}) + scale_fill_manual(values=c($col_list)) + theme(legend.position=\"none\")
#"
done
for obs in $count; do
  echo "## $obs ##"
  echo -n "barplot(c("
  var=''
  for exp in $mus_exp; do
    var="${var}length(my.mus[which(my.mus\$tag==\"$exp\" & my.mus\$$obs==1),]\$nam) / length(my.mus[which(my.mus\$tag==\"$exp\"),]\$nam),"
  done
  var="${var}length(my.mus[which(my.mus\$tag==\"mouse\" & my.mus\$$obs==1),]\$nam) / length(my.mus[which(my.mus\$tag==\"mouse\"),]\$nam)),ylab=\"$obs Fraction\",names=c("
  for exp in $mus_exp; do var="${var}\"$exp\","; done
  var="${var}\"mouse\"),col=c($col_list))"
  echo "$var"
done


echo "dev.off()"

echo "#
###############
# AA Plot
###############
pdf(\"age_aa.pdf\")"

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
  for exp in $mus_exp; do
    var="${var}wilcox.test(my.mus[which(my.mus\$tag=='${exp}'),]\$$obs, my.mus[which(my.mus\$tag=='mouse'),]\$$obs)\$p.val,"
  done
  var="${var})"; var="${var%,)})"
  echo "$var"
done

for obs in $count; do
  echo -n "p_$obs <- c("
  var=''
  for exp in $mus_exp; do
    hit_sample="length(my.mus[which(my.mus\$tag==\"$exp\" & my.mus\$$obs==1),]\$nam)"
    hit_left="(length(my.mus[which(my.mus\$tag==\"mouse\" & my.mus\$$obs==1),]\$nam) - length(my.mus[which(my.mus\$tag==\"$exp\" & my.mus\$$obs==1),]\$nam))"
    miss_sample="length(my.mus[which(my.mus\$tag==\"$exp\" & my.mus\$$obs!=1),]\$nam)"
    miss_left="(length(my.mus[which(my.mus\$tag==\"mouse\" & my.mus\$$obs!=1),]\$nam) - length(my.mus[which(my.mus\$tag==\"$exp\" & my.mus\$$obs!=1),]\$nam))"
   my_test="fisher.test(matrix(c($hit_sample,$hit_left,$miss_sample,$miss_left),2,2),alternative='two.sided')\$p.value"
   var="${var}$my_test,"
  done
  var="${var})"; var="${var%,)})"
  echo "$var"
done

echo -n "p_mus <- rbind("
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
for exp in $mus_exp; do 
  echo "$exp <- p.adjust(p_mus[,$cnt],method=\"hochberg\")"
  cnt=$((cnt+1))
done
echo -n "corr_pval<-cbind("
var=''
for exp in $mus_exp; do 
 var="${var}$exp,"
done 
var="${var})"
var="${var%,)})"
echo "$var"
echo " 
ptab<-signif(corr_pval,digits=2)
tab_top    <- head(ptab,20)
tab_bottom <- ptab[-(1:20),]

pdf(\"age_pval.pdf\")
grid.table(tab_top)
grid.newpage()
grid.table(tab_bottom)
dev.off()
#"

