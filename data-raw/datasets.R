#'PRISTDataPackage for PRIST paper
#'
#'This package contains data used for PRIST paper.
#'@docType package
#'@author Jingyuan Deng, Greg Finak, Andrew McDavid, Masaano Yajima
#'@name PRISTDataPackage
#'@title Data package for Single Cell RNA sequence paper
#'@aliases PRISTDataPackage-package
NULL

sys.source('ProcessAlex.R',envir=topenv())
sys.source('ProcessMAIT.R',envir=topenv())
#sys.source('ProcessMAIT_bulk.R',envir=topenv())
#sys.source('ProcessMonocle.R',envir=topenv())
#sys.source('ProcessCMV.R',envir=topenv())
#keepDataObjects(c("sca_monocle","sca_alex","sca_mait","sca_mait_raw","sca_mait_bulk","sca_mait_bulk_raw","sca_cmv"))
keepDataObjects(c("sca_alex","sca_mait"))
