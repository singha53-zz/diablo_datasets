#################################################
#
# asthmaDatasets.R
# Author: Amrit Singh
# Date: April 01, 2016
#
#################################################
WhereAmI <- "/Users/asingh/Dropbox/Manuscript/diablo/analyses/casestudy2_asthma/data/"

# load libraries
library(annotate)
library("hugene10sttranscriptcluster.db")
library(org.Hs.eg.db)
library(KEGG.db)
library(dplyr)
library(mixOmics)

# import datasets
genes <- read.csv(paste0(WhereAmI, "data/genRMA.csv"), header=TRUE, row.names=1)
genSym <- getSYMBOL(rownames(genes), "hugene10sttranscriptcluster")
xCell <- cbind(genSym=genSym, genes[,-1])
xCell <- xCell[!is.na(genSym), ]
write_tsv(xCell, paste0(WhereAmI, "data/xCell_genes.tsv"))

## import xCell results
xCell_pval <- read.delim(paste0(WhereAmI, "data/xCell_xCell_genes_xCell_1918013018.pvals.txt"), header=TRUE, row.names=1)
plot(as.numeric(as.matrix(xCell_pval)))
gplots::heatmap.2(as.matrix(xCell_pval), trace = "none")
xCell <- read.delim(paste0(WhereAmI, "data/xCell_xCell_genes_xCell_1918013018.txt"), header=TRUE, row.names=1)
colnames(xCell) <- colnames(genes)[-1]

demo <- read.csv(paste0(WhereAmI, "data/cbc.csv"), header=TRUE, row.names=1)[, 1:2]
cells <- read.csv(paste0(WhereAmI, "data/cbc.csv"), header=TRUE, row.names=1)[, -c(1:2)]
metabolites <- read.csv(paste0(WhereAmI, "data/metabolomics.csv"), header=TRUE, row.names=1)
dim(genes); dim(demo); dim(cells); dim(metabolites); 
all(rownames(demo) == colnames(genes)[-1])
all(rownames(demo) == rownames(cells))
all(rownames(demo) == colnames(metabolites[,-c(1:10)]))

## correspondence between cells and inferred cells
xCell <- t(xCell) %>% as.data.frame
result <- pls(X = cells, Y = xCell)
cim(result)
par(mfrow = c(3, 3))
plot(cells$Erythrocytes..x10.12.L. ~ xCell$Erythrocytes )
plot(cells$Platelets..x10.9.L. ~ xCell$Platelets)
plot(cells$Relative.Neutrophils ~ xCell$Neutrophils)
plot(cells$Relative.Lymphocytes ~ xCell$`CD4+ T-cells`)
plot(cells$Relative.Monocytes ~ xCell$Monocytes)
plot(cells$Relative.Eosinophils ~ xCell$Eosinophils)
plot(cells$Relative.Basophils ~ xCell$Basophils)

## combine both datasets
#cells <- cbind(cells, xCell)  # 90 variables

## datasets to use
genExp <- genes[, -1]
metExp <- metabolites[,-c(1:10)]

## Convert gene expression dataset to KEGG modules dataset
eg <- getEG(rownames(genExp), "hugene10sttranscriptcluster")
kegg <- org.Hs.egPATH2EG
mapped <- mappedkeys(kegg)
kegg2 <- as.list(kegg[mapped])

kegg.affyid <- lapply(kegg2, function(x) rownames(genExp)[!is.na(match(eg, x))])
genesDat <- genExp[unlist(kegg.affyid), ]
dim(genesDat)
names(kegg.affyid)

moduleLabels <- list()
for(i in 1:length(kegg.affyid)){
  moduleLabels[[i]] <- rep(names(kegg.affyid)[i], length(kegg.affyid[[i]]))
}

moduleColors <- unlist(moduleLabels)
xx <- as.list(KEGGPATHID2NAME)
unlist(xx)
modNames <- as.character(unlist(xx)[moduleColors])

genDat <- as.data.frame(genesDat)
genDat$Mod <- modNames

modComp1 <- genDat %>% group_by(Mod) %>% dplyr::do(x = mixOmics::pca(as.matrix(t(.[, -ncol(.), drop = FALSE])), ncomp = 1)$x)
genMEs0 = do.call(cbind, modComp1$x)
colnames(genMEs0) <- modComp1$Mod

## metabolite modules
metDat <- as.data.frame(metExp)
metDat$Mod <- metabolites$SUB_PATHWAY
modComp1 <- metDat %>% group_by(Mod) %>% do(x = pca(t(.[, -ncol(.)]), ncomp = 1)$x)
metMEs0 = do.call(cbind, modComp1$x)
colnames(metMEs0) <- modComp1$Mod

## extract within-var matrix
dim(cells)    # 28 90
dim(genMEs0)  #28 229
dim(metMEs0)  #28 60
all(rownames(cells) == rownames(genMEs0))
all(rownames(cells) == rownames(metMEs0))

rownames(cells)
sample = substr(rownames(cells), 1, 7)
time = factor(substr(rownames(cells), 9, 12), level = c('pre', 'post'))
design = data.frame(sample, time)

gene.module=genMEs0
metabolite.module = metMEs0

save(cells=cells, gene.module=gene.module, metabolite.module=metabolite.module, metExp=metExp, genExp=genExp, 
metabolites=metabolites, genDat=genDat, demo=demo, file=paste0(WhereAmI, "asthmaDatasets.RDATA"))

# for mixOmics website
# --------------------
# data from March
load("/Users/klecao/Dropbox/Manuscript/mixOmics.org:DIABLO/1-Data & Preprocessing/asthma/data/asthmaDatasets.RDATA")

dim(cells)    # 28 23
dim(genMEs0)  #28 229
dim(metMEs0)  #28 60

rownames(cells)
sample = substr(rownames(cells), 1, 7)
time = factor(substr(rownames(cells), 9, 12), level = c('pre', 'post'))
design = data.frame(sample, time)

gene.module=genMEs0
metabolite.module = metMEs0


save(cells = cells, gene.module, metabolite.module, design = design, file = "asthma.mixDIABLO.RDATA")

# rm and test
rm(list = ls())
load('asthma.mixDIABLO.RDATA')

dim(cells)
dim(gene.module)
dim(metabolite.module)
summary(design$time)
