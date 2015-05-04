#'Processed single cell RNA seqdata for our single cell RNA seq paper from our MAIT study.
#'
#'This data set contains the processed data for our fake 
#'assay, with 100 replicated measurements for each of four subjects.
#'@docType data
#'@format This is a SingleCellAssay object
#'\describe{
#'\item{cData}{Cell level feature annotation.}
#' \itemize{
#' \item{wellKey} unique identifier
#' \item{condition} stimulated or unstimulated
#' \item{ncells} number of cells in the well
#' \item{ngeneson} percentage of genes that are turned on
#' \item{cngeneson} centered ngeneson
#' \item{TRAV1} TRAV1
#' \item{TRBV6} TRBV6
#' \item{TRBV4} TRBV4
#' \item{TRBV20} TRBV20
#' \item{alpha} alpha chain
#' \item{beta}  beta chain
#' \item{ac}    alpha chain + stimulation
#' \item{bc}    beta chain + stimulation
#'	}
#'\item{fData}{Gene level feature annotation.}s
#' \itemize{
#' \item{primerid} entrez id
#' \item{entrez} entrez id
#' \item{symbolid} gene symbols
#'	}
#'\item{exprs}{Expression matrix}
#'}
#'@source The data was from Matin
#'@note We used UCSC hg19 as the reference genome and used RSEM (bowtie) to align and calculate the count and the TPM matrix.
#'@author Jingyuan Deng, Masanao Yajima, Andrew McDavid, Greg Finak
#'@name sca_mait
#'@title Create MAIT data
NULL

# library
suppressPackageStartupMessages({
	require(SingleCellAssay)
	require(plyr)
	require(data.table)
	require(biomaRt)
	require(stringr)
	require(reshape2)
	require(org.Hs.eg.db)
	require(Biobase)
	require(abind)
})
# .testPackageVersion <- function(package,minversion){
#     version <- as.numeric(as.character(ldply(strsplit(as.character(packageVersion(package)),"\\."))[2]))
#     if(version < minversion) {
#         stop(paste0("Package ",package,"must be at least version ",minversion,"\n You have version ",version))
#     }else{
#         message("Found ",package," ",version)
#     }
# }

# .testPackageVersion(package="SingleCellAssay",minversion=921)

#The data path should be set appropriately
data_path <- "/shared/silo_researcher/Gottardo_R/SingleCellRNA_Data/MAIT"

#
## UCSC HG19
## generate count and tpm matrix from RSEM results
rsem.dir <- "/shared/silo_researcher/Gottardo_R/jingyuan_working/MAIT/processedData/Rsem_UCSC"

res.files  <- list.files(rsem.dir, pattern = "*.genes.results$", full = T, recursive = T)
res        <- lapply(res.files, function(x) fread(x))
names(res) <- gsub(".genes.results", "", basename(res.files))
res_long   <- ldply(res)
setnames(res_long, 3, "transcript")
count     <- dcast(res_long, transcript ~ .id, value.var="expected_count")
count_mat <- as.matrix(count[ ,-1])
tpm        <- dcast(res_long, transcript ~ .id, value.var="TPM")
tpm_mat    <- as.matrix(tpm[ ,-1])

len_mat            <- dcast(res_long, transcript ~ .id, value.var="length")[,-1]
effective_len_mat  <- dcast(res_long, transcript ~ .id, value.var="effective_length")[,-1]

## BULK
bulk_rsem_dir <- "/shared/silo_researcher/Gottardo_R/jingyuan_working/MAIT_bulk/RSEM"

res.files.bulk  <- list.files(bulk_rsem_dir, pattern = "*.genes.results$", full = T, recursive = T)
res.bulk        <- lapply(res.files.bulk, function(x) fread(x))
names(res.bulk) <- gsub(".genes.results", "", basename(res.files.bulk))
res.bulk_long   <- ldply(res.bulk)
setnames(res.bulk_long, 3, "transcript")
count_bulk      <- dcast(res.bulk_long, transcript ~ .id, value.var="expected_count")
count_mat_bulk    <- as.matrix(count_bulk[ ,-1])
tpm_bulk        <- dcast(res.bulk_long, transcript ~ .id, value.var="TPM")
tpm_mat_bulk    <- as.matrix(tpm_bulk[ ,-1])

len_bulk_mat        <- dcast(res.bulk_long, transcript ~ .id, value.var="length")[,-1]
effective_len_bulk_mat        <- dcast(res.bulk_long, transcript ~ .id, value.var="effective_length")[,-1]

# get symbolid
idsymbol      <-read.delim("/shared/silo_researcher/Gottardo_R/SingleCellRNA_Data/SCRAMDataPackage/inst/extdata/IdToSymbol",stringsAsFactors=FALSE)
idsymdb       <- data.table(idsymbol)
setnames(idsymdb ,c("UCKGID","SYMBOL"))
setkey(idsymdb ,"UCKGID")
tranid         <- sapply( strsplit( tpm[ ,1], ","), function(x) x[1] ) # only need the first transcript
idsymtab       <- data.frame( idsymdb[tranid,])
stopifnot( all(idsymtab[,1]==tranid)) # double check

fd = data.frame(symbolid= idsymtab$SYMBOL,transcrpitid=tpm[ ,1],stringsAsFactors=FALSE)

# filter
#count_mat  <- as.matrix(count[ ,-1])
#count_mat_filt <-count_mat[!is.na(fDat.dt$primerid)&rowMeans(tpm_mat>0)>0.1,]
#tpm_mat_filt <- tpm_mat#[rowMeans(tpm_mat>0)>0.1,]
#fd_filt      <- fd#[rowMeans(tpm_mat>0)>0.1,] 
esetAll      <- abind(tpm=t(log2(tpm_mat+1)),count=t(log2(count_mat+1)), rev.along=0)

## cdata
qc_info      <- read.table(paste(data_path,"mait_quality.txt",sep="/"),header=TRUE,stringsAsFactors = FALSE)

samples      <- colnames(tpm_mat)
condition    <- sapply(strsplit(samples, "-"), function(x) x[3])
ngo          <- rowMeans(esetAll>0)
## tcr
tcr <- read.table('/shared/silo_researcher/Gottardo_R/SingleCellRNA_Data/MAIT/MAIT_TCR.tsv', check.names=FALSE, as.is=TRUE)
alpha <- str_detect(colnames(tcr), fixed('TRAV'))
alphamax <- colnames(tcr)[alpha][apply(tcr[,alpha,drop=FALSE], 1, function(x){
    wm <- which.max(x)
    if(x[wm]>0) wm else NA
}                                      )]
betamax <- colnames(tcr)[!alpha][apply(tcr[,!alpha,drop=FALSE], 1, function(x){
    wm <- which.max(x)
    if(x[wm]>0) wm else NA
}
                                       )]
tcr.dat <- data.frame(tcr[,colMeans(tcr>0)>.1], alpha=factor(alphamax), beta=factor(betamax))
alphatab <- table(tcr.dat$alpha, exclude=NULL)
alphatabnm <- ifelse(table(tcr.dat$alpha, exclude=NULL)>4, names(alphatab), 'other')
levels(tcr.dat$alpha) <- alphatabnm
betatab <- table(tcr.dat$beta, exclude=NULL)
betatabnm <- ifelse(table(tcr.dat$beta, exclude=NULL)>4, names(betatab), 'other')
levels(tcr.dat$beta) <- betatabnm
tcr.dat$beta <- relevel(tcr.dat$beta, 'other')
tcr.dat$alpha <- relevel(tcr.dat$alpha, 'TRAV1')


#cDat <- data.frame(qc_info[samples,], samples, condition)
cdat <- data.frame(wellKey=samples, condition=condition, qc_info[samples,], ncells=1, ngeneson=ngo, cngeneson=ngo-mean(ngo), tcr.dat[rownames(esetAll),], stringsAsFactors=FALSE)
cdat <- with(cdat, cbind(cdat, data.frame(ac=interaction(alpha, condition), bc=interaction(beta, condition))))

# low quality cell filter
cdat$ourfilter =  cdat$pastFastqc=="PASS"& cdat$exonRate >0.3 & cdat$PercentToHuman>0.6 

## get entrezgene ids
entrez <- as.data.frame( org.Hs.egSYMBOL )
setDT(entrez)
setkey(entrez, symbol)
entrez.order  <- entrez[fd$symbolid,,mult='first']
entrez.order[is.na(gene_id),gene_id:=symbol]
stopifnot(all(entrez.order$symbol==colnames(esetAll)))
entrez.order[,primerid:=make.unique(gene_id)]

dimnames(esetAll)[2]      <- list(primerid=entrez.order$primerid)
#dimnames(tpm_mat_filt)[1] <- list(primerid=entrez.order$primerid)
#dimnames(count_mat_filt)[1] <- list(primerid=entrez.order$primerid)

fdat <- data.frame(primerid=entrez.order$primerid, entrez=entrez.order$gene_id, symbolid=entrez.order$symbol, stringsAsFactors=FALSE)
sca <- FromMatrix('SingleCellAssay', esetAll, cdat, fdat)
dimnames(sca)[3] <- dimnames(esetAll)[3]
sca_filt <- sca[,freq(sca)>.1]
cData(sca_filt)$ngeneson <- rowMeans(exprs(sca_filt)>0)
cData(sca_filt)$cngeneson <- cData(sca_filt)$ngeneson- mean(cData(sca_filt)$ngeneson)

sca_mait_all <- primerAverage(sca_filt, 'symbolid')
sca_mait_all@featureData <- sca_filt@featureData[row.names(fData(sca_mait_all)),]
stopifnot(all(fData(sca_mait_all)$primerid==colnames(exprs(sca_mait_all))))

sca_mait <- subset( sca_mait_all, cData(sca_mait_all)$ourfilter == TRUE )
sca_mait_raw <- sca
# BULK
cell_count  <- factor(unlist(lapply(strsplit(colnames(tpm_mat_bulk),"-"),"[",1)),levels=c("25", "50",  "100", "250", "500", "1000" ))
stimulation <- unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(tpm_mat_bulk),"-"),"[",3)),"_"),"[",1))
stimulation <- ifelse(stimulation=="NonStim","Unstim","Stim")

samples_bulk      <- colnames(tpm_mat_bulk)
ngo_bulk          <- colMeans(tpm_mat_bulk>0)

#cDat <- data.frame(qc_info[samples,], samples, condition)
cdat <- data.frame(wellKey=samples_bulk, condition=stimulation, ncells=cell_count, ngeneson=ngo_bulk, cngeneson=ngo_bulk-mean(ngo_bulk), stringsAsFactors=FALSE)

esetBULK      <- abind(tpm=t(log2(tpm_mat_bulk+1)), rev.along=0)
dimnames(esetBULK)[2]      <- list(primerid=entrez.order$primerid)
fdat$in_singlecell <- freq(sca)>.1

sca_bulk_all <- FromMatrix('SingleCellAssay', esetBULK, cdat, fdat)
sca_bulk_filt <- sca_bulk_all[,freq(sca)>.1]
sca_mait_bulk <- primerAverage(sca_bulk_filt, 'symbolid')
cData(sca_mait_bulk)$ngeneson <- rowMeans(exprs(sca_mait_bulk)>0)
cData(sca_mait_bulk)$cngeneson <- cData(sca_mait_bulk)$ngeneson- mean(cData(sca_mait_bulk)$ngeneson)

esetBULK2      <- abind(tpm=t(log2(tpm_mat_bulk+1)), count=t(log2(count_mat_bulk+1)), rev.along=0)
dimnames(esetBULK2)[2]      <- list(primerid=entrez.order$primerid)

sca_mait_bulk_raw <- FromMatrix('SingleCellAssay', esetBULK2, cdat, fdat)

# sca_mait<-addlayer(sca_mait,"count")

# layer(sca_mait) <- 'count'
# exprs(sca_mait) <-  t(count_mat_filt)
# layer(sca_mait) <- 'tpm'
# exprs(sca_mait) <-   t(log2(tpm_mat_filt+1))

#saveRDS(sca_mait,"sca_mait.rds")

#dimnames(sca_mait)[3] <- list(layer=c('count', 'tpm'))
