#'Processed single cell RNA seqdata for our single cell RNA seq paper from Alex study.
#'
#'This data set contains the processed data for our fake 
#'assay, with 100 replicated measurements for each of four subjects.
#'@docType data
#'@format This is a SingleCellAssay object
#'\describe{
#' \itemize{
#' \item{cData}{Cell level feature annotation.}
#' \item{wellKey} unique identifier
#' \item{Stim} Stimulation
#' \item{Time} Time 
#' \item{BioRep} Boolean for biological replicates
#' \item{TechRep} Boolean for technical replicates
#' \item{OtherCondition} Boolean for other conditions 
#' \item{ngeneson} percentage of genes that are turned on
#' \item{cngeneson} centered ngeneson
#'	}
#'\item{fData}{Gene level feature annotation.}
#' \itemize{
#' \item{primerid} gene symbols
#' \item{GENES} gene symbols
#'	}
#'\item{exprs}{Expression matrix}
#'}
#'@source The data was from Matin
#'@note We used UCSC hg19 as the reference genome and used RSEM (bowtie) to align and calculate the count and the TPM matrix.
#'@author Jingyuan Deng, Masanao Yajima, Andrew McDavid, Greg Finak
#'@name sca_alex
#'@title Create Alex data 
NULL

suppressPackageStartupMessages({
library(org.Mm.eg.db)
library(GO.db)
library(plyr)
library(SingleCellAssay)
library(abind)
library(arm)
library(RColorBrewer)
library(data.table)
library(biomaRt)

})
.testPackageVersion <- function(package,minversion){
    version <- as.numeric(as.character(ldply(strsplit(as.character(packageVersion(package)),"\\."))[2]))
    if(version < minversion) {
        stop(paste0("Package ",package,"must be at least version ",minversion,"\n You have version ",version))
    }else{
        message("Found ",package," ",version)
    }
}

#.testPackageVersion(package="SingleCellAssay",minversion=914)

#The data path should be set appropriately
data_path <- "/shared/silo_researcher/Gottardo_R/SingleCellRNA_Data/Alex"
# Load Alex data
# Primary data
# AllSuppLibs <- read.table(file.path(data_path,"SCSII_SuppTable3_Final.txt"),sep="\t",header=TRUE)
AllSuppLibs <- fread(file.path(data_path,"SCSII_SuppTable3_Final.txt"))
# AllSuppLibs <- data.table(AllSuppLibs)
# construct SingleCellAssay and annotate

genes <- AllSuppLibs[,GENES]
tpm_matrix <- data.frame(AllSuppLibs)[,-1]

libname <- colnames(tpm_matrix)
rownames(tpm_matrix) <- genes

ids= genes
ids.idx  <- 1:length(ids)

# entrez id
entrez <- as.data.frame(org.Mm.egSYMBOL)
setDT(entrez)	
setkey(entrez, symbol)
entrez.order <- entrez[ids,,mult='first']
entrez.order[is.na(gene_id),gene_id:=symbol]
stopifnot(all(entrez.order$symbol==ids))
entrez.order[,primerid:=make.unique(gene_id)]

xx       <- as.list(org.Mm.egALIAS2EG)
names(xx)<- toupper(names(xx))
xx.id    <- xx[ids]
names( xx.id ) <-ids
xx.id.1  <- lapply(xx.id,`[`,1)
xx.id.fix<- lapply( xx.id.1, function(x) ifelse( is.null(x), NA, x ) )
entrez<-unlist(xx.id.fix)

fdat <- data.frame(primerid=genes, entrez = entrez, symbolid=genes, stringsAsFactors=FALSE)
rownames(fdat) <- genes
#annotate
CD <- data.table(wellKey=libname)
CD[,libname:=gsub("Unstimulated","Unstimulated_0h",libname)]
CD <- data.table(CD,Stim=ldply(strsplit(as.character(CD$libname),"_"),function(x)x[1]),Time=ldply(strsplit(as.character(CD$libname),"_"),function(x)x[2]))
CD[,BioRep :=libname%like%"Biological_Replicate"]
CD[,TechRep:=libname%like%"Technical_Replicate"]
CD[,OtherCondition:=libname%like%"Tenth"]
CD[,OtherCondition:=OtherCondition|libname%like%"Golgi"]

setnames(CD,c("wellKey","libname","Stim","Time","BioRep","TechRep","OtherCondition"))
cdat           <- as.data.frame(CD)
rownames(cdat) <- libname

sca_alex_all <- FromMatrix('SingleCellAssay', t(tpm_matrix), cData= cdat, fData=fdat)


#saveRDS(sca_alex_all,file="data/sca_alex_allSCA.rds")

#Filter for Cluster Disruption cells
CDFILTER <- !(exprs(sca_alex_all[,"LYZ1"])<6|exprs(sca_alex_all[,"SERPINB6B"])>4)
CDFILTER <- rownames(CDFILTER)[CDFILTER]
sca_alex_all <- sca_alex_all[CDFILTER,]

# EE <- 2^(exprs(sca_alex_all))-1 #exponentiate for thresholding

# #Threshold
# thresholded         <- thresholdSCRNACountMatrix(EE,nbins=100,min_per_bin=100,bin_by="median")
# sca_alex_all        <- addlayer(sca_alex_all, 'atpm')
# layer(sca_alex_all) <- 'atpm'
# exprs(sca_alex_all) <- thresholded$counts_threshold

#annotate with ngeneson after filtering
ngeneson            <- apply(exprs(sca_alex_all),1,function(x)mean(x>0))
CD                  <- cData(sca_alex_all)
CD$ngeneson         <- ngeneson
CD$cngeneson        <- CD$ngeneson-mean(ngeneson)
rownames(CD)        <- CD$libname
cData(sca_alex_all) <- CD

# 10% expression filter
sca_alex_all_global <- sca_alex_all[cData(sca_alex_all)$Stim%in%c("PAM","PIC","LPS","Unstimulated")&cData(sca_alex_all)$BioRep==FALSE&cData(sca_alex_all)$TechRep==FALSE&cData(sca_alex_all)$OtherCondition==FALSE,]
gene_filter         <- names( which( apply( exprs( sca_alex_all_global ), 2, function(x)mean(x>0)>0.1 ) ) )
sca_alex            <- sca_alex_all_global[,gene_filter]

#AllSuppLibs <- data.table:::melt.data.table(AllSuppLibs)
#setnames( AllSuppLibs, c( "GENES", "libname", "tpm" ) )
#sca_alex_all <- SingleCellAssay( AllSuppLibs,idvars="libname", primerid="GENES", measurement="tpm" )

# filter_by_stim <- function(x){
#     levs <- c("PIC","PAM","LPS")
#     tofilter <- NULL
#     for(i in levs){
#         filt <- subset(x,Stim%in%i&!OtherCondition&!TechRep&!BioRep)
#         tofilter <- c(tofilter,names(which(apply(exprs(filt),2,function(x)mean(x>0)<=0.1))))
#     }

#     setdiff(featureData(x)@data$primerid,unique(tofilter))
# }
# gene_filter_by_stim <- filter_by_stim(sca_alex_all)

# #saveRDS(sca_alex_all,file="data/sca_alex_allProcessed.rds")

# sca_alex_all_by_stim <- sca_alex_all[,gene_filter_by_stim]
# 
#saveRDS(sca_alex_all_by_stim,file="data/sca_alex_allProcessed_Filtered_byStim.rds")
#saveRDS(sca_alex_all_global,file="data/AlexProcessed_Filtered_globally.rds")

# li                    <- unlist(as.list(org.Mm.egSYMBOL2EG))
# uniqGene              <- data.table(gene_id=li, tk=names(li))
# uniqGene[,primerid:=toupper(tk)]
# setkey(uniqGene, primerid)
# uniqGeneSca           <- uniqGene[fData(sca_alex)$primerid,]
# uniqGenedf            <- as.data.frame(uniqGeneSca)
# uniqGenedf$tk         <- uniqGenedf$primerid <- NULL
# row.names(uniqGenedf) <- uniqGeneSca$primerid
# sca_alex              <- combine(sca_alex, uniqGenedf)

# go <- as.data.table(as.data.frame(org.Mm.egGO))
# setkey(go, gene_id)
# setkey(uniqGene, gene_id)
# setkey(uniqGeneSca, gene_id)
# scaGO <- go[uniqGeneSca,,nomatch=0]
# scaGO[,tk:=NULL]
# anyGO <- go[uniqGene,,nomatch=0]
