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

## ----filter_outlying_cell,results='hide'---------------------------------
sca <- subset(sca,nGeneOn > 4000)
eid <-
select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys = fData(sca_mait)$entrez,keytype ="GENEID",columns = c("GENEID","TXNAME"))
ueid <- unique(na.omit(eid)$GENEID)
sca <- sca[,fData(sca)$entrez %in% ueid]

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

## ----zlm, dependson="data",results='hide'--------------------------------
options(mc.cores=LARGE_MEMORY_CORES)
# ZLM (ridge regression for continuous, get standardized deviance residuals)
cond<-factor(cData(sca)$condition)
cond<-relevel(cond,"Unstim")
cData(sca)$condition<-cond
zlm <- zlm.SingleCellAssay(~condition + cngeneson, sca, method = "bayesglm", 
    ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"))
zlm.nongeneson <- zlm.SingleCellAssay(~condition , sca, method = "bayesglm", 
    ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"))
zlm.notrt <- zlm.SingleCellAssay(~cngeneson , sca, method = "bayesglm", 
    ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"))
zlm.null <- zlm.SingleCellAssay(~1 , sca, method = "bayesglm", 
    ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"))
# LRT
lrt <- lrTest(zlm, "condition")
lrt.nongeneson<-lrTest(zlm.nongeneson,"condition")

# Still need to merge on logFC
res_gene <- data.table(melt(lrt))
res_gene[,primerid:=gsub("\\.\\d+","",primerid)]
res_gene <- merge(res_gene, fd, by="primerid")
res_gene_hurdle <- res_gene[metric=="Pr(>Chisq)" & test.type=="hurdle"]
res_gene_hurdle[,adj:=p.adjust(value,"fdr")]
lfc<-getLogFC(zlm)[contrast=="conditionStim"]
setkey(lfc,primerid)
setkey(res_gene_hurdle,primerid)
res_gene_hurdle<-merge(lfc,res_gene_hurdle)

# merge on logFC
res_gene_nong <- data.table(melt(lrt.nongeneson))
res_gene_nong[,primerid:=gsub("\\.\\d+","",primerid)]
res_gene_nong <- merge(res_gene_nong, fd, by="primerid")
res_gene_hurdle_nong <- res_gene_nong[metric=="Pr(>Chisq)" & test.type=="hurdle"]
res_gene_hurdle_nong[,adj:=p.adjust(value,"fdr")]
lfc_nong<-getLogFC(zlm.nongeneson)[contrast=="conditionStim"]
setkey(lfc_nong,primerid)
setkey(res_gene_hurdle_nong,primerid)
res_gene_hurdle_nong<-merge(lfc_nong,res_gene_hurdle_nong)

nrow(res_gene_hurdle_nong[adj<0.01])
nrow(res_gene_hurdle[adj<0.01])

## ----supp_fig_percent_deviance,results='hide'----------------------------
library(scales)

panelb=rbind(data.table(melt((zlm.notrt@deviance-zlm@deviance)/zlm.null@deviance),model="trt"),data.table(melt((zlm.nongeneson@deviance-zlm@deviance)/zlm.null@deviance),model="cdr"))
setnames(panelb,c("gene","component","value","model"))
panelb=dcast(panelb,...~model)

#Shalek data for this SFig
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
cData(filtered)$Stim <- factor(cData(filtered)$Stim,levels=c("Unstimulated","LPS","PAM","PIC"))
cData(filtered)$Time <- factor(cData(filtered)$Time)
filtered_nobaseline <- subset(filtered,!Stim%in%c("Unstimulated")) #drop Time 0
cData(filtered_nobaseline)$Stim <- factor(cData(filtered_nobaseline)$Stim)
cData(filtered_nobaseline)$Time <- factor(cData(filtered_nobaseline)$Time)
ng <- zlm.SingleCellAssay(~cngeneson+Stim/Time,sca=filtered_nobaseline,method="ridge",ebayes=TRUE)
nong <- zlm.SingleCellAssay(~Stim/Time,sca=filtered_nobaseline,method="ridge",ebayes=TRUE)
notrt <- zlm.SingleCellAssay(~cngeneson,sca=filtered_nobaseline,method="ridge",ebayes=TRUE)
nullmodel = zlm.SingleCellAssay(~1,sca=filtered_nobaseline,method="ridge",ebayes=TRUE)


panelc=rbind(data.table(melt((notrt@deviance-ng@deviance)/nullmodel@deviance),model="trt"),data.table(melt((nong@deviance-ng@deviance)/nullmodel@deviance),model="cdr"))
setnames(panelc,c("gene","component","value","model"))
panelc=dcast(panelc,...~model)
#ggplot(panelc)+aes(x=trt,y=cdr,col=component)+geom_point()

toplot = rbind(cbind(panelc,dataset="mdc"),cbind(panelb,dataset="MAIT"))

p3 = ggplot(toplot)+geom_hex()+aes(x=cdr,y=trt)+facet_grid(component~dataset)+theme_linedraw()+scale_color_discrete("Component",labels=c("Continouos","Discrete"))+scale_x_continuous("Proportion of deviance explained by CDR")+scale_y_continuous("Proportion of deviance explained by treatment")+scale_fill_distiller(type = "div",palette = 2,space = "Lab",trans="log",breaks=c(1,10,50,100,250,500))
pdf(file = "../inst/extdata/output/Supplementary_Figure_2.pdf",width = 8,height =
      5)
      print(p3)
      dev.off()
      print(p3)

## ---- dependson="data;zlm",fig.width=8,fig.height=8----------------------
# Let's look at the top 40
genes_to_plot <- res_gene_hurdle[abs(logFC)>FCTHRESHOLD&adj<0.01][order(-abs(logFC),adj)][1:40,symbolid]
ggplot(dt_long[gene_id%in%genes_to_plot][,gene_id:=factor(gene_id,levels=genes_to_plot)],aes(x=condition, y=value,color=condition))+geom_violin()+geom_jitter()+facet_wrap(~gene_id)+ggtitle("Top 40 DE Genes in Activated MAIT Cells")

## ---- dependson="data;zlm",fig.width=15,fig.height=10--------------------
fd<-fData(sca)
fd$primerid<-fd$symbolid
colnames(sca)<-fd$symbolid
sca@featureData@data<-fd
FCTHRESHOLD<-0.58
genes_to_plot <- res_gene_hurdle[abs(logFC)>FCTHRESHOLD&adj<0.01][order(-abs(logFC),adj)][1:100,symbolid]
COMPASS::pheatmap(exprs(sca[,genes_to_plot]),row_annotation=cData(sca)[,"condition",drop=FALSE],main="Top 100 DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)))

## ----boots, dependson="data;zlm"-----------------------------------------
#bootstrap
cache(boots <- bootVcov1(zlm, 50))

## ----gsea, dependson="data;zlm;boots",results='hide'---------------------
module_file <- list.files("../inst/extdata/", pattern = module, full.names = TRUE)
gene_set <- getGmt(module_file)
gene_ids <- geneIds(gene_set)
gene_ids=gene_ids[!names(gene_ids)%like%"TBA"&!names(gene_ids)%like%"B cell"]
saveRDS(gene_ids,file="../inst/extdata/MAIT_Modules_symbols.rds")
sets_indices <- ids2indices(gene_ids, fd$symbolid)
# Only keep modules with at least min_gene_in_module
sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]
cache(gsea <- gseaAfterBoot(zlm, boots, sets_indices, CoefficientHypothesis("conditionStim")))
t_stat <- melt(calcZ(gsea, testType = "normal"))
t_stat_wide <- dcast(t_stat, set ~ comp + metric)
t_stat <- data.frame(calcZ(gsea, testType = "normal", combined = "stouffer"))
t_stat$set <- rownames(t_stat)
t_stat_comb <- data.table(merge(t_stat_wide, t_stat, by = "set"))
setorder(t_stat_comb,P)
t_stat_comb[,adj:=p.adjust(P,"fdr")]

##Filter by effect size: discrete odds ratio of 1
t_stat_comb=t_stat_comb[set%in%names(which(abs(gsea[,"disc","stat","test"]-gsea[,"disc","stat","null"])>log(2.5) ))|set%in%names(which(abs(gsea[,"cont","stat","test"]-gsea[,"cont","stat","null"])>log2(1.5)))]

## ---- eval=TRUE, dependson="data;zlm;boots;gsea"-------------------------
for (i in c(1:10,12,13)) {
    genes_to_plot <- unlist(gene_ids[names(gene_ids) == as.character(t_stat_comb$set[i])])
    gp <- ggplot(dt_long[gene_id %in% genes_to_plot], aes(x = condition, y = value, 
        color = condition)) + geom_violin() + geom_jitter() + facet_wrap(~gene_id,scale="free") + 
        labs(title = paste0(t_stat_comb$set[i], ", Z=", round(t_stat_comb$Z[i], 
            2), ", P=", round(t_stat_comb$P[i], 2)))
    print(gp)
}

## ----scores_fit, dependson="data;gsea;zlm",fig.width=10,fig.height=10----
options("mc.cores" = LARGE_MEMORY_CORES)

score_4_hook <- function(x) {
  if (all(x@fitted)) {
    class(x@fitC) <- c("glm","lm")
    class(x@fitD) <- c("bayesglm","glm","lm")
    wh <- !colnames(x@modelMatrix) %like% "cngeneson"
    wh2 <- !names(coef(x@fitC)) %like% "cngeneson"
    fc <-
      x@modelMatrix[,!wh,drop = FALSE] %*% coef(x@fitC)[!wh2,drop = FALSE] #continuous CDR effect
    fd <-
      arm::invlogit(x@modelMatrix[,c("(Intercept)","cngeneson"),drop = FALSE] %*%
                      coef(x@fitD)[c("(Intercept)","cngeneson"),drop = FALSE]) #discrete CDR effect
    R <-
      matrix((x@response - fc * fd),nrow = 1) #residuals corrected for ngeneson in the continuous part
    colnames(R) <- names(residuals(x@fitD))
    R[x@response == 0] <- 0 - fd
    R
  }
}

zlm.scores <-
  zlm.SingleCellAssay(
    ~condition + cngeneson, sca, method = "bayesglm",
    ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"),hook =
      score_4_hook
  )

## ----scores, dependson="data;gsea;scores_fit;zlm",fig.width=10,fig.height=10----
scores <- do.call(rbind,zlm.scores@hookOut)
rownames(scores) <- names(zlm.scores@hookOut)

genes_to_plot <- sapply(1:nrow(t_stat_comb),function(i) {
unlist(gene_ids[names(gene_ids) == as.character(t_stat_comb$set[i])])
})
names(genes_to_plot) <- t_stat_comb$set

##Include the DE MAIT signature
mait_signature_genes_up = res_gene_hurdle[abs(logFC) > FCTHRESHOLD &
adj < 0.01][order(-abs(logFC),adj)][1:50,][logFC > 0,symbolid]
mait_signature_genes_down = res_gene_hurdle[abs(logFC) > FCTHRESHOLD &
adj < 0.01][order(-abs(logFC),adj)][1:50,][logFC < 0,symbolid]

genes_to_plot$"MAIT Activation [Up]" <- mait_signature_genes_up
genes_to_plot$"MAIT Activation [Down]" <- mait_signature_genes_down

genes_to_plot <- ldply(lapply(genes_to_plot,function(x)
data.frame(x)))
colnames(genes_to_plot) <- c("set","gene")
scores <- (melt(scores))

colnames(scores) <- c("gene","cell","score")

scores <- (merge(genes_to_plot,scores,by = "gene",all.y = TRUE))
condition <- cData(sca)[,"condition",drop = FALSE]
condition$cell <- rownames(condition)
scores <- (merge(scores,condition,by = "cell",all.x = TRUE))
scores <- na.omit(scores)
scores <- data.table(scores)
scores = scores[,set := factor(set)]
#t_stat_comb = t_stat_comb[,adj := p.adjust(P)]
scores$set <-
sapply(as.character(scores$set),function(x)
Kmisc::wrap(x,width = 25))
sets <-
sapply(t_stat_comb[order(-abs(t_stat_comb$Z),t_stat_comb$adj),set][c(1:9)],function(x)
Kmisc::wrap(x,25))
#sets <- gsub("\n$","",gsub(".\\(M.+$","",sets[!sets %like% "TBA"]))
gsea_scores <-
ggplot(scores[set %in% sets,list(score =
mean(score)),list(cell,set,condition)][,set := factor(set,levels = sets)]) +  geom_violin(show_guide =
FALSE) + geom_jitter(show_guide = FALSE) + aes(x = condition,y = score,col =
condition) + facet_wrap( ~set) + theme_linedraw() + theme(
plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line =
element_line(colour = "black"),legend.position = "bottom"
) + scale_color_manual("",values = c(brewer.pal(name = "PiYG",n = 3)[c(1)],"#5729EE"))

## ----residuals,dependson="data;zlm",include=FALSE,eval=TRUE,results='hide'----
zlm.resid <- zlm.SingleCellAssay(~condition + cngeneson, sca, method = "bayesglm", 
    ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"),hook=deviance_residuals_hook)

residuals<-do.call(rbind,zlm.resid@hookOut)
rownames(residuals)<-fData(sca)[,"symbolid"]
colnames(residuals)<-rownames(cData(sca))

genes_to_plot <- res_gene_hurdle[abs(logFC)>FCTHRESHOLD&adj<0.01][order(-abs(logFC),value)][1:20,symbolid]

#correlation amongst cells
C.s<-cor(residuals[genes_to_plot,cData(sca)[,"condition"]%in%"Stim"])
C.u<-cor(residuals[genes_to_plot,cData(sca)[,"condition"]%in%"Unstim"])
C<-cor(residuals[genes_to_plot,])

COMPASS::pheatmap(C,row_annotation=cData(sca)[,"condition",drop=FALSE],show_rownames=TRUE,show_colnames=FALSE,main="Correlation of residuals across all cells using top 20 DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))

## ----nongeneson,dependson="data;zlm",include=FALSE,eval=TRUE-------------
options("mc.cores"=LARGE_MEMORY_CORES)
zlm.resid.nong <- zlm.SingleCellAssay(~condition , sca, method = "bayesglm", 
    ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"),hook=deviance_residuals_hook)
lrt.nong <- lrTest(zlm.resid.nong, "condition")

# Still need to merge on logFC
res_gene_nong <- data.table(melt(lrt.nong))
res_gene_nong <- merge(res_gene_nong, fd, by="primerid")
res_gene_hurdle_nong <- res_gene_nong[metric=="Pr(>Chisq)" & test.type=="hurdle"]
res_gene_hurdle_nong[,adj:=p.adjust(`value`,"fdr")]
lfc_nong<-getLogFC(zlm.resid.nong)[contrast=="conditionStim"]
setkey(lfc_nong,primerid)
setkey(res_gene_hurdle_nong,primerid)
res_gene_hurdle_nong<-merge(lfc_nong,res_gene_hurdle_nong)


residuals.nong<-do.call(rbind,zlm.resid.nong@hookOut)
rownames(residuals.nong)<-fData(sca)[,"symbolid"]
colnames(residuals.nong)<-rownames(cData(sca))

genes_to_plot <- res_gene_hurdle_nong[abs(logFC)>FCTHRESHOLD&adj<0.01][order(-abs(logFC),value)][1:20,symbolid]
nrow(res_gene_hurdle_nong[abs(logFC)>FCTHRESHOLD&adj<0.01])
#correlation amongst cells
C<-cor(residuals.nong[genes_to_plot,])

COMPASS::pheatmap(C,row_annotation=cData(sca)[,"condition",drop=FALSE],show_rownames=FALSE,show_colnames=FALSE,clustering_method="ward",main="Correlation of residuals across all cells using top 20 DE genes (without ngeneson)",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))

## ----residualCor,dependson="data;residuals;zlm",eval=TRUE,include=FALSE,echo=FALSE----
residuals<-do.call(rbind,zlm.resid@hookOut)
rownames(residuals)<-fData(sca)[,"symbolid"]
colnames(residuals)<-rownames(cData(sca))
genes_to_plot <- res_gene_hurdle[abs(logFC)>FCTHRESHOLD&adj<0.01][order(-abs(logFC),value)][1:50,symbolid]
CC.s<-cor(t(residuals)[cData(sca)[,"condition"]%in%"Stim",genes_to_plot])
CC.u<-cor(t(residuals)[cData(sca)[,"condition"]%in%"Unstim",genes_to_plot])
CC.s<-data.table(melt(CC.s))
CC.u<-data.table(melt(CC.u))
setnames(CC.s,c("row","col","cor.s"))
setnames(CC.u,c("row","col","cor.u"))
CC.s<-CC.s[!is.na(cor)]
CC.u<-CC.u[!is.na(cor)]
setkey(CC.s,row,col)
setkey(CC.u,row,col)
CC.su<-merge(CC.s,CC.u)
CC.su[,cor.d:=cor.s-cor.u]
setorder(CC.su,cor.d)

toplot<-dcast(CC.su,row~col,value.var ="cor.s")[,-c(1)]
rownames(toplot)<-colnames(toplot)
o<-hclust(as.dist(1-cor(toplot)),method="ward")$order


toplot.ns<-dcast(CC.su,row~col,value.var ="cor.u")[,-c(1)]
rownames(toplot.ns)<-colnames(toplot.ns)
COMPASS::pheatmap(toplot.ns[o,o],col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21),cluster_rows=FALSE,cluster_cols=FALSE,legend=FALSE,fontsize = 7)
nostim<-grid.grab()
COMPASS::pheatmap(toplot[o,o],col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21),cluster_rows=FALSE,cluster_cols=FALSE,fontsize = 7)
stim<-grid.grab()
cutree(hclust(as.dist(1-cor(toplot)),method="ward"),3)
genes<-c("IFNG","GZMB","JUND","CD69","FOS","GZMA","TXN")

## ----nongresidcor,dependson="data;nongeneson;zlm",eval=FALSE,include=FALSE,echo=FALSE----
#  # genes_to_plot <- res_gene_hurdle_nong[abs(logFC)>FCTHRESHOLD&adj<0.01][order(-abs(logFC),value)][1:50,symbolid]
#  # CC.s<-cor(t(residuals.nong)[cData(sca)[,"condition"]%in%"Stim",genes_to_plot])
#  # CC.u<-cor(t(residuals.nong)[cData(sca)[,"condition"]%in%"Unstim",genes_to_plot])
#  # CC.s<-data.table(melt(CC.s))
#  # CC.u<-data.table(melt(CC.u))
#  # setnames(CC.s,c("row","col","cor.s"))
#  # setnames(CC.u,c("row","col","cor.u"))
#  # CC.s<-CC.s[!is.na(cor)]
#  # CC.u<-CC.u[!is.na(cor)]
#  # setkey(CC.s,row,col)
#  # setkey(CC.u,row,col)
#  # CC.su<-merge(CC.s,CC.u)
#  # CC.su[,cor.d:=cor.s-cor.u]
#  # setorder(CC.su,cor.d)
#  #
#  #
#  # toplot<-dcast(CC.su,row~col,value.var ="cor.u")[,-c(1)]
#  # rownames(toplot)<-colnames(toplot)
#  # COMPASS::pheatmap(toplot,main="Correlation among top DE genes in unstimulated cells (nongneson)",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))
#  #
#  # toplot<-dcast(CC.su,row~col,value.var ="cor.s")[,-c(1)]
#  # rownames(toplot)<-colnames(toplot)
#  # COMPASS::pheatmap(toplot,main="Correlation among top DE genes in stimulated cells (no ngeneson)",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))

## ----pca,dependson="data;residuals;",fig.width=10,include=FALSE----------
library(grid)
PCbiplot <- function(PC, x="PC1", y="PC2", colors=c('black', 'black', 'red', 'red'),arrow_colors="red",point_colors="black",arrow_alpha=0.75,gene_subset=NULL) {
    # PC being a prcomp object
    data <- data.frame(obsnames=row.names(PC$x), PC$x)
    plot <- ggplot(data, aes_string(x=x, y=y)) +geom_point(color=point_colors)#+ geom_text(alpha=.8, size=3, aes(label=obsnames), color=point_colors)
    plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2, color=colors[2])
    datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
    mult <- min(
        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
    datapc <- transform(datapc,
            v1 = .7 * mult * (get(x)),
            v2 = .7 * mult * (get(y))
            )
    plot <- plot + coord_equal() + geom_text(data=datapc[gene_subset,], aes(x=v1, y=v2, label=varnames), size = 4, vjust=1, color=arrow_colors,alpha=arrow_alpha)
    plot <- plot + geom_segment(data=datapc[gene_subset,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=arrow_alpha, color=arrow_colors)
    plot
}
residuals<-do.call(rbind,zlm.resid@hookOut)
rownames(residuals)<-fData(sca)[,"symbolid"]
colnames(residuals)<-rownames(cData(sca))
point_colors<-factor(cData(sca)[,"condition"],labels=c(Unstim=brewer.pal(name="PiYG",n=1)[1],Stim="#5729EE"))

genes_to_plot <- res_gene_hurdle[abs(logFC)>FCTHRESHOLD&adj<0.01][order(-abs(logFC),value)][1:50,symbolid]
pcaplot<-PCbiplot(prcomp(t(residuals[genes_to_plot,])),point_colors=point_colors,arrow_colors="blue",arrow_alpha=0.75,gene_subset=genes)+theme_bw()

## ----network,dependson="data;zlm;residuals",fig.width=10,eval=FALSE,echo=FALSE,include=FALSE,eval=FALSE----
#  # library(Matrix)
#  # library(igraph)
#  # library(glasso)
#  # residuals<-do.call(rbind,zlm.resid@hookOut)
#  # rownames(residuals)<-fData(sca)[,"symbolid"]
#  # colnames(residuals)<-rownames(cData(sca))
#  #
#  # getNetwork<-function(residuals=NULL,genes=NULL,rho=0.35,main=NULL,plot=TRUE){
#  #   g<-glasso(cor(t(residuals[genes,]),method="kendall"),rho=rho)
#  #   colnames(g$w)<-genes
#  #   rownames(g$w)<-genes
#  #   g<-graph.adjacency(g$w,mode="undirected",diag=FALSE,weighted=TRUE)
#  #   sort(igraph:::degree(g))
#  #   G<-subgraph(g,which(igraph:::degree(g)>0))
#  #   set.seed(10)
#  #   lay<-layout.fruchterman.reingold(G,area=vcount(G)^2.3,repulserad=vcount(G)^2)
#  #   if(plot){
#  #     plot(G,layout=lay,edge.width=3,vertex.size=18,vertex.label.cex=0.5,main=main)
#  #   }
#  #   G
#  # }
#  #
#  # genes_to_plot <- res_gene_hurdle[,adj:=p.adjust(value,"fdr")][abs(logFC)>3&adj<0.01][order(-abs(logFC),adj)][,symbolid]
#  # G<-getNetwork(residuals[,cData(sca)[,"condition"]%in%c("Stim","Unstim")],genes_to_plot,0.35,main="Residuals network in MAIT")
#  #
#  #
#  # COMPASS::pheatmap(exprs(sca[,vertex.attributes(G)$name]),row_annotation=cData(sca)[,"condition",FALSE],main="Genes from residuals network in MAIT (with ngeneson)",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))
#  # COMPASS::pheatmap(cor(t(residuals[vertex.attributes(G)$name,])),main="gene-gene correlation from network estimation (with ngeneson)",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))
#  # gene_colors<-c("red","blue")[as.numeric(factor(sign(res_gene_hurdle[symbolid%in%vertex.attributes(G)$name,][order(symbolid),logFC])))]
#  # PCbiplot(prcomp((residuals[sort(vertex.attributes(G)$name),])),arrow_color=NA,point_colors=gene_colors)+ggtitle("PCA of residuals of genes from network estimate with ngeneson")

## ----network_nong, dependson="data;zlm;nongeneson",fig.width=10,eval=FALSE,echo=FALSE,include=FALSE,eval=FALSE----
#  # residuals.nong<-do.call(rbind,zlm.resid.nong@hookOut)
#  # rownames(residuals.nong)<-fData(sca)[,"symbolid"]
#  # colnames(residuals.nong)<-rownames(cData(sca))
#  # genes_to_plot <- res_gene_hurdle_nong[,adj:=p.adjust(value,"fdr")][abs(logFC)>3&adj<0.01][order(-abs(logFC),adj)][,symbolid]
#  # G<-getNetwork(residuals.nong[,cData(sca)[,"condition"]%in%c("Stim","Unstim")],genes_to_plot,0.35,main="Residuals network in MAIT (no ngeneson correction)")
#  #
#  #
#  # gene_colors<-c("red","blue")[as.numeric(factor(sign(res_gene_hurdle_nong[symbolid%in%vertex.attributes(G)$name,][order(symbolid),logFC])))]
#  #
#  # COMPASS::pheatmap(exprs(sca[,vertex.attributes(G)$name]),row_annotation=cData(sca)[,"condition",FALSE],main="Genes from residuals network in MAIT (without ngeneson)",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))
#  # COMPASS::pheatmap(cor(t(residuals.nong[vertex.attributes(G)$name,])),main="Gene-gene correlation from network estimation (without ngeneson)",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)),breaks = seq(-1,1,l=21))
#  # PCbiplot(prcomp((residuals.nong[vertex.attributes(G)$name,])),arrow_color=NA,point_color=gene_colors)+ggtitle("PCA of residuals of genes from network estimate without ngeneson")

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

