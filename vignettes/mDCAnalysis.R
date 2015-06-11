## ----load,echo=FALSE-----------------------------------------------------
rm(list=ls())
suppressPackageStartupMessages({
  
library(GGally)
library(grid)
library(ggplot2)
library(reshape2)
library(org.Hs.eg.db)
library(plyr)
library(glasso)
library(data.table)
library(GO.db)
library(hom.Hs.inp.db)
library(MAST)
library(Matrix)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(org.Mm.eg.db)
library(GSEABase)
library(corpcor)
library(Rtsne)
library(MASTDataPackage)

})
data(MASTDataPackage)
data_dir <- "data"
if(packageVersion("MAST")>="0.930"){
  message("Version Okay")
}else{
  stop("Wrong SingleCellAssay Version")
}
FCTHRESHOLD<-log2(1.5)
plotheme<-theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))
knitr::opts_chunk$set(list(echo=FALSE,eval=TRUE,message=FALSE,error=FALSE,warning=FALSE,fig.width=10,fig.height=8,dev=c("png")))

## ----filtering_and_modules, cache=FALSE,echo=FALSE-----------------------
filtered<-sca_alex
percent_expressed<-0.1
thres<-"adapt"
tt <- thresholdSCRNACountMatrix(2^exprs(filtered) - 1, nbins = 20, min_per_bin = 30)

if (thres == "adapt") {
    exprs(filtered) <- tt$counts_threshold
} else if (thres == "fixed") {
    mat <- exprs(filtered)
    mat[mat < 1/log(2)] <- 0
    exprs(filtered) <- mat
}

expressed_genes <- colMeans(exprs(filtered) > 0) > percent_expressed
filtered <- filtered[, expressed_genes]
data <- exprs(filtered)
cd <- cData(filtered)
fd <- fData(filtered)

dt <- data.table(gene_id = fd$symbolid, t(data))
dt_long <- data.table(data.table:::melt.data.table(dt))

load("../inst/extdata/clusters_shalek.rda")
CORE_ANTIVIRAL<-as.character(subset(clusters,CLUSTER=="Id")$GENE)
PEAKED_INFLAM<-as.character(subset(clusters,CLUSTER=="IIIc")$GENE)
SUSTAINED_INFLAM<-as.character(subset(clusters,CLUSTER=="IIId")$GENE)

ids=featureData(filtered)$primerid
ids.idx <- 1:length(ids)

#get all the GO IDS for all the genes
sym2go<-AnnotationDbi:::select(org.Mm.eg.db,keys=gsub("^(.)","\\U\\1",tolower(ids),perl=TRUE),columns="GOALL",keytype="SYMBOL")
sym2go <- na.omit(sym2go)
sym2go <- data.table(sym2go)
ids <- data.table(cbind(ids,ids.idx))
setkey(ids,ids)
sym2go=sym2go[,ids:=SYMBOL]
setkey(sym2go,ids)
sym2go <- merge(ids,sym2go)
#module_member is the gene index for the given module
BP <- sym2go[ONTOLOGYALL%in%"BP"]
BP <- BP[,module_member:=ids.idx,list(GOALL)]

## ----zlm,dependson="filtering_and_modules", cache=FALSE------------------
cData(filtered)$Stim <- factor(cData(filtered)$Stim,levels=c("Unstimulated","LPS","PAM","PIC"))
cData(filtered)$Time <- factor(cData(filtered)$Time)
filtered_nobaseline <- subset(filtered,!Stim%in%c("Unstimulated")) #drop Time 0
cData(filtered_nobaseline)$Stim <- factor(cData(filtered_nobaseline)$Stim)
cData(filtered_nobaseline)$Time <- factor(cData(filtered_nobaseline)$Time)
options(mc.cores=8)
fit.bystim <- zlm.SingleCellAssay(~cngeneson+Stim/Time,sca=filtered_nobaseline,method="ridge",ebayes=TRUE,hook=deviance_residuals_hook,lambda=0.1)
fit<-fit.bystim

## ----no_ngeneson, cache=FALSE, dependson="filtering_and_modules"---------
options(mc.cores=7)
#Fit a model without ngeneson for comparison
fit.bystim.nongeneson <- zlm.SingleCellAssay(~Stim/Time,sca=filtered_nobaseline,method="ridge",ebayes=TRUE,hook=deviance_residuals_hook,lambda=0.1)

## ----any_time_effect,dependson="zlm", cache=FALSE------------------------
#Test for any time effect
M <- matrix(0,nrow=ncol(coef(fit,"D")))
rownames(M) <- colnames(coef(fit,"D"))
M[colnames(coef(fit,"D"))%like%"LPS",] <- 1
M.time <- M
M.time['(Intercept)',]<- 0
anyTime <- lrTest(fit, hypothesis=M.time)
anyTime.sorted <- na.omit(anyTime[order(anyTime[,'hurdle','Pr(>Chisq)']),'hurdle',])
lfc.anytime<-getLogFC(fit)[contrast%like%"Time"]
res_gene_hurdle<-data.table(dcast(melt(anyTime.sorted),primerid~metric))
res_gene_hurdle=res_gene_hurdle[,adj:=p.adjust(`Pr(>Chisq)`,"fdr")]

#test with no ngeneson
nong.test<-lrTest(fit.bystim.nongeneson, hypothesis=M.time[-2,,drop=FALSE])
nong.test <- na.omit(nong.test[order(nong.test[,'hurdle','Pr(>Chisq)']),'hurdle',])
nong.test<-data.table(dcast(melt(nong.test),primerid~metric))
nong.test=nong.test[,adj:=p.adjust(`Pr(>Chisq)`,"fdr")]


tg<-merge(lfc.anytime[contrast%like%"LPS"&contrast%like%"Time6h"],res_gene_hurdle,by="primerid")[order((logFC),adj)][1:100,primerid]
COMPASS::pheatmap(exprs(filtered_nobaseline[cData(filtered_nobaseline)$Time%in%c("1h","6h")&cData(filtered_nobaseline)$Stim%in%"LPS",as.character(tg)]))

## ----any_stim_effect,dependson="zlm",eval=FALSE,echo=FALSE---------------
#  ### Test for any stimulation effect
#  #Test for any Stim effect
#  M <- matrix(0,nrow=ncol(coef(fit,"D")))
#  rownames(M) <- colnames(coef(fit,"D"))
#  M[colnames(coef(fit,"D"))%like%"Stim",] <- 1
#  M.stim <- M
#  anyStim <- lrTest(fit, hypothesis=M.stim)
#  anyStim.sorted <- na.omit(anyStim[order(anyStim[,'hurdle','Pr(>Chisq)']),'hurdle',])

## ----mouse_go------------------------------------------------------------
library(data.table)
library(limma)
library(GO.db)
gene_association= fread("../inst/extdata/gene_association.mgi",skip=6)
geneset_id      = split(toupper(gene_association$V3), gene_association$V5)
geneset_terms   = Term(GOTERM)[names(geneset_id)]
geneset_index   = ids2indices(geneset_id, fData(filtered_nobaseline)$symbolid)

