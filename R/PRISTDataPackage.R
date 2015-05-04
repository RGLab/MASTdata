#'PRISTDataPackage for PRIST paper
#'
#'This package contains data used for PRIST paper.
#'@docType package
#'@author Jingyuan Deng, Greg Finak, Andrew McDavid, Masaano Yajima
#'@name PRISTDataPackage
#'@title Data package for Single Cell RNA sequence paper
#'@aliases PRISTDataPackage-package
NULL
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
