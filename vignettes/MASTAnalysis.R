## ----eval=FALSE----------------------------------------------------------
#  # library(devtools)
#  # install_git("http://github.com/RGLab/MASTDataPackage")
#  # install_git('http://github.com/RGLab/MAST')

## ----echo=FALSE, error=FALSE---------------------------------------------
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <-
function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
require(grid)

# Make a list from the ... arguments and plotlist
plots <- c(list(...), plotlist)

numPlots = length(plots)

# If layout is NULL, then use 'cols' to determine layout
if (is.null(layout)) {
# Make the panel
# ncol: Number of columns of plots
# nrow: Number of rows needed, calculated from # of cols
layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
ncol = cols, nrow = ceiling(numPlots / cols))
}

if (numPlots == 1) {
print(plots[[1]])

} else {
# Set up the page
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

# Make each plot, in the correct location
for (i in 1:numPlots) {
# Get the i,j matrix positions of the regions that contain this subplot
matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

print(plots[[i]], vp = viewport(
layout.pos.row = matchidx$row,
layout.pos.col = matchidx$col
))
}
}
}
suppressPackageStartupMessages({
library(MIMOSA)
library(plyr)
library(reshape)
library("ggplot2")
library("MASTDataPackage")
library("MAST")
library("GSEABase")
library("limma")
library(grid)
library("reshape2")
library(RColorBrewer)
library("data.table")
library("knitr")
library(psych)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(flowClust)
library(COMPASS)
library(scales)
library(grid)
library(stringr)
library(gridExtra)
library(class)
options(mc.cores = detectCores() - 1) #or set smaller if you find yourself running out of memory
})
LARGE_MEMORY_CORES <- 7 #for large memory operations
knitr::opts_chunk$set(list(
echo = FALSE,eval = TRUE,message = FALSE,error = FALSE,warning = FALSE,cache = TRUE))
plotheme <-
theme(
plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line =
element_line(colour = "black")
)


## ---- cache=FALSE--------------------------------------------------------
opts_chunk$set(tidy=TRUE, cache=TRUE, messages=TRUE)

## ---- dependson="data", cache=FALSE,results='hide'-----------------------
# thresholding
tt <- thresholdSCRNACountMatrix(2^exprs(sca) - 1, nbins = 20, min_per_bin = 30)

if (thres == "adapt") {
    exprs(sca) <- tt$counts_threshold
} else if (thres == "fixed") {
    mat <- exprs(sca)
    mat[mat < 1/log(2)] <- 0
    exprs(sca) <- mat
}

expressed_genes <- colMeans(exprs(sca) > 0) > percent_expressed
sca <- sca[, expressed_genes]
data <- exprs(sca)
cd <- cData(sca)
fd <- fData(sca)

## Data in long format for visualization
dt <- data.table(gene_id = fd$symbolid, t(data))
dt_long <- data.table(melt(dt))
dt_long[, `:=`(condition, ifelse(grepl("-Stim", variable), "Stim", "Unstim"))]

## ----background_correlation,cache=FALSE,dependson="residuals;data;nongeneson;zlm",results='hide'----
library(stringr)
unstim<-colnames(residuals.nong)%like%"Unstim"

S1<-cor(t(residuals.nong[,unstim]))
S1.s<-cor(t(residuals.nong[,!unstim]))

S2<-cor(t(residuals[,unstim]))
S2.s<-cor(t(residuals[,!unstim]))

extractCor<-function(x,genes,class,model){
  S<-x[genes,genes]
  S<-S[upper.tri(S,diag=FALSE)]
  gn<-matrix(rep(genes,length(genes)),ncol=length(genes),byrow=TRUE)
  gn1<-gn[upper.tri(gn,diag=FALSE)]
  gn<-matrix(rep(genes,length(genes)),ncol=length(genes),byrow=FALSE)
  gn2<-gn[upper.tri(gn,diag=FALSE)]
  data.table(rho=S,gene1=gn1,gene2=gn2,class,model)
}

C<-extractCor(S1,rownames(S1),class = "no_ngeneson",model = "no_ngeneson")
C.s<-extractCor(S1.s,rownames(S1.s),class = "no_ngeneson",model = "no_ngeneson")
C2<-extractCor(S2,rownames(S2),class = "ngeneson",model = "ngeneson")
C2.s<-extractCor(S2.s,rownames(S2.s),class = "ngeneson",model = "ngeneson")

C3<-rbind(C,C2)
C3.s<-rbind(C.s,C2.s)

stopifnot(all.equal(C[,list(gene1, gene2)],C2[,list(gene1, gene2)]))
stopifnot(all.equal(C.s[,list(gene1, gene2)],C2.s[,list(gene1, gene2)]))
Cboth <- rbind(data.table("Rho (background)"=C[,rho], "Rho (CDR)"=C2[,rho], Condition='Non-stimulated'),
               data.table("Rho (background)"=C.s[,rho], "Rho (CDR)"=C2.s[,rho], Condition='Stimulated'))
Cboth <- Cboth[sample(nrow(Cboth), 1e5),]



C3 = C3[,delrho := rho[class == "no_ngeneson"] - rho[class == "ngeneson"],list(gene1,gene2)]
C3.s[,delrho := rho[class == "no_ngeneson"] - rho[class == "ngeneson"],list(gene1,gene2)]
setkey(C3,delrho)
setkey(C3.s,delrho)

C3[,bin:=cut(delrho,breaks = 1000)]
C3 = C3[!is.na(bin)]
C3sumry=C3[,.N,bin]
C3sumry[,cdf:=cumsum(N)/sum(N)]
C3sumry[,bin_c:=mean(as.numeric(str_split_fixed(pattern=",",gsub("]","",gsub("\\(","",bin)),2))),bin]

# p = ggplot(C3sumry) + geom_line(aes(x = bin_c,y = cdf,group = 1)) + geom_vline(x =
#                                                                                  0) + theme_linedraw() + theme(axis.text.x = element_text(angle = 90,hjust =
#                                                                                  1),aspect.ratio = 1) + scale_x_continuous("Difference in correlation (without CDR - with CDR)",breaks =
#                                                                                  seq(-1,1,l = 21)) + scale_y_continuous("fraction of gene pairs")+ggtitle("Unstimulated")
#  

C3.s[,bin:=cut(delrho,breaks = 1000)]
C3.s = C3.s[!is.na(bin)]
C3s.sumry=C3.s[,.N,bin]
C3s.sumry[,cdf:=cumsum(N)/sum(N)]
C3s.sumry[,bin_c:=mean(as.numeric(str_split_fixed(pattern=",",gsub("]","",gsub("\\(","",bin)),2))),bin]
       
p=ggplot(data.table(rbind(C3s.sumry,C3sumry),Stim=rep(c("Stim","Unstim"),c(nrow(C3s.sumry),nrow(C3sumry)))))+geom_line(aes(x=bin_c,y=cdf,col=Stim))+ geom_vline(x =
                                                                                 0) + theme_linedraw() + theme(axis.text.x = element_text(angle = 90,hjust =
                                                                                 1),aspect.ratio = 1) + scale_x_continuous("Difference in correlation (without CDR - with CDR)",breaks =
                                                                                 seq(-1,1,l = 21)) + scale_y_continuous("fraction of gene pairs")
# 
# p.s = ggplot(C3s.sumry) + geom_line(aes(x = bin_c,y = cdf,group = 1)) + geom_vline(x =
#                                                                                  0) + theme_linedraw() + theme(axis.text.x = element_text(angle = 90,hjust =
#                                                                                  1),aspect.ratio = 1) + scale_x_continuous("Difference in correlation (without CDR - with CDR)",breaks =
#                                                                                  seq(-1,1,l = 21)) + scale_y_continuous("fraction of gene pairs")+ggtitle("Stimulated")
#  
#                                                                          
# p<-ggplot(C3)+aes(x=delrho)+stat_ecdf()+theme_linedraw()+ggtitle("Unstimulated")+theme(aspect.ratio=1)
# p.s<-ggplot(C3.s)+aes(x=delrho)+stat_ecdf()+theme_linedraw()+ggtitle("Stimulated")+theme(aspect.ratio=1)

library(gridExtra)

colors<-c(brewer.pal(name = "PiYG",n=3)[c(1)],"#5729EE")
layout<-matrix(1,ncol=1)
pdf("../inst/extdata/output/Supplementary_Figure_10.pdf",width=8,height=4)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=1,ncol=1)))
print(p+plotheme+scale_color_discrete("Condition",labels=c("Stimulated","Non-stimulated")),vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
#print(p.s+plotheme+scale_color_discrete("Model",labels=c("CDR adjusted","No CDR")),vp=viewport(layout.pos.row = 1,layout.pos.col = 2))
dev.off()

cond_de_genes_ng<-res_gene_hurdle[adj<0.01&abs(logFC)>FCTHRESHOLD,symbolid]
cond_de_genes_nong<-res_gene_hurdle_nong[adj<0.01&abs(logFC)>FCTHRESHOLD,symbolid]
foo<-data.table(S2[cond_de_genes_ng,cond_de_genes_ng][upper.tri(S2[cond_de_genes_ng,cond_de_genes_ng],diag=FALSE)])
ng.signif.cors.unstim<-nrow(foo)*foo[,mean(p.adjust(r.test(n=sum(unstim),r12=V1)$p,"fdr")<0.01)]
foo<-data.table(S2.s[cond_de_genes_ng,cond_de_genes_ng][upper.tri(S2.s[cond_de_genes_ng,cond_de_genes_ng],diag=FALSE)])
ng.signif.cors.stim<-nrow(foo)*foo[,mean(p.adjust(r.test(n=sum(!unstim),r12=V1)$p,"fdr")<0.01)]

foo<-data.table(S1[cond_de_genes_nong,cond_de_genes_nong][upper.tri(S1[cond_de_genes_nong,cond_de_genes_nong],diag=FALSE)])
nong.signif.cors.unstim<-nrow(foo)*foo[,mean(p.adjust(r.test(n=sum(unstim),r12=V1)$p,"fdr")<0.01,na.rm=TRUE)]
foo<-data.table(S1.s[cond_de_genes_nong,cond_de_genes_nong][upper.tri(S1.s[cond_de_genes_nong,cond_de_genes_nong],diag=FALSE)])
nong.signif.cors.stim<-nrow(foo)*foo[,mean(p.adjust(r.test(n=sum(!unstim),r12=V1)$p,"fdr")<0.01)]


C<-extractCor(S1[cond_de_genes_nong,cond_de_genes_nong],rownames(S1[cond_de_genes_nong,cond_de_genes_nong]),class = "no_ngeneson",model = "no_ngeneson")
C.s<-extractCor(S1.s[cond_de_genes_nong,cond_de_genes_nong],rownames(S1.s[cond_de_genes_nong,cond_de_genes_nong]),class = "no_ngeneson",model = "no_ngeneson")
C2<-extractCor(S2[cond_de_genes_ng,cond_de_genes_ng],rownames(S2[cond_de_genes_ng,cond_de_genes_ng]),class = "ngeneson",model = "ngeneson")
C2.s<-extractCor(S2.s[cond_de_genes_ng,cond_de_genes_ng],rownames(S2.s[cond_de_genes_ng,cond_de_genes_ng]),class = "ngeneson",model = "ngeneson")


# #################
highly.correlated.genes.ng.unstim<-unique(unlist(C2[data.table(S2[cond_de_genes_ng,cond_de_genes_ng][upper.tri(S2[cond_de_genes_ng,cond_de_genes_ng],diag=FALSE)])[,which((p.adjust(r.test(n=sum(unstim),r12=V1)$p,"fdr")<0.01))],list(gene1,gene2)],use.names=FALSE))
highly.correlated.genes.ng.stim<-unique(unlist(C2.s[data.table(S2.s[cond_de_genes_ng,cond_de_genes_ng][upper.tri(S2.s[cond_de_genes_ng,cond_de_genes_ng],diag=FALSE)])[,which((p.adjust(r.test(n=sum(!unstim),r12=V1)$p,"fdr")<0.01))],list(gene1,gene2)],use.names=FALSE))
highly.correlated.genes.nong.unstim<-unique(unlist(C[data.table(S1[cond_de_genes_nong,cond_de_genes_nong][upper.tri(S1[cond_de_genes_nong,cond_de_genes_nong],diag=FALSE)])[,which((p.adjust(r.test(n=sum(unstim),r12=V1)$p,"fdr")<0.01))],list(gene1,gene2)],use.names=FALSE))
highly.correlated.genes.nong.stim<-unique(unlist(C.s[data.table(S1.s[cond_de_genes_nong,cond_de_genes_nong][upper.tri(S1.s[cond_de_genes_nong,cond_de_genes_nong],diag=FALSE)])[,which((p.adjust(r.test(n=sum(!unstim),r12=V1)$p,"fdr")<0.01))],list(gene1,gene2)],use.names=FALSE))

sum(highly.correlated.genes.ng.stim%in%c(mait_signature_genes_up,mait_signature_genes_down))
sum(highly.correlated.genes.nong.stim%in%c(mait_signature_genes_up,mait_signature_genes_down))

sum(highly.correlated.genes.ng.unstim%in%c(mait_signature_genes_up,mait_signature_genes_down))
sum(highly.correlated.genes.nong.unstim%in%c(mait_signature_genes_up,mait_signature_genes_down))

