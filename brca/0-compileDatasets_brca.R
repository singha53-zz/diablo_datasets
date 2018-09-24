#################################################
#
# compileDatasets_brca.R
# Date: April 01, 2016
# Author: Amrit Singh
#
#################################################
WhereAmI <- "/Volumes/Data/mac files Dec 05 2016 update/Oz.TRI/TRI/Breast Cancer/NewAnalysis/data/"

## Import datasets
###----------------------------
#
# 1) Clinical Data
#
###----------------------------
# Clinical dataset
clinical <- as.data.frame(t(read.delim(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2015110100.0.0/BRCA.clin.merged.txt"))))
colnames(clinical) <- as.character(as.matrix(clinical[1,]))
clinDat <- clinical[-1,]
rownames(clinDat) <- gsub("-", ".", toupper(clinDat$patient.patient_id))
length(rownames(clinDat)); length(unique(rownames(clinDat))); dim(clinDat);

###----------------------------
#
# Genes
#
###----------------------------
## gene expression
genExp <- read.delim(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"))
genExp2 <- genExp[-c(1:30), ]
genNames <- unlist(lapply(strsplit(as.character(genExp2[, 1]), "|", fixed=TRUE), function(i) i[1]))
genNames[which(genNames == "SLC35E2")[2]] <- "SLC35E2.rep"
rownames(genExp2) <- genNames
genEset <- genExp2[, -1]
dim(genEset)
dim(genEset)
table(colnames(genEset))
table(unlist(lapply(strsplit(colnames(genEset), "\\."), function(i) i[4])))  # 1098 01A

genEset2 <- genEset[, grep("01A", colnames(genEset))]
dim(genEset); dim(genEset2);
colnames(genEset2) <- unlist(lapply(strsplit(colnames(genEset2), "\\."), function(i) i[3]))

# number of common subjects between clinical and gene expression
length(intersect(rownames(clinDat), colnames(genEset2)))   ## 1080

###----------------------------
#
# Proteins
#
###----------------------------
prot <- read.delim(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.Merge_protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.Level_3.2015110100.0.0/BRCA.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt"))[-1,]
rownames(prot) <- prot[, "Sample.REF"]
prot1 <- prot[, -1]
dim(prot1)   #  142 x 410

prot2 <- matrix(0, nrow=nrow(prot1), ncol=ncol(prot1))
rownames(prot2) <- rownames(prot1)
colnames(prot2) <- colnames(prot1)
for(i in 1:nrow(prot2)){
  prot2[i, ] <- as.numeric(as.matrix(prot1[i, ]))
}

prot3 <- prot2[, grep("01A", colnames(prot2))]
dim(prot3)  # 142 proteins x 403 samples
colnames(prot3) <- unlist(lapply(strsplit(colnames(prot3), "\\."), function(x) x[3]))

# number of common subjects between clinical and gene expression
length(Reduce(intersect, list(rownames(clinDat), colnames(genEset2), colnames(prot3))))   ## 400

###----------------------------
#
# microRNA
#
###----------------------------
### miRNA
## HiSeq
mirnaHiseq <- read.delim(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2015110100.0.0/BRCA.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt"))
dim(mirnaHiseq)
mirnaHiseq2 <- mirnaHiseq[-1, seq(2,ncol(mirnaHiseq), by=3)]
rownames(mirnaHiseq2) <- mirnaHiseq[-1, 1]
dim(mirnaHiseq2)
length(unique(colnames(mirnaHiseq2)))

## GA
mirnaGA <- read.delim(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.Merge_mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2015110100.0.0/BRCA.mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt"))
dim(mirnaGA)
mirnaGA2 <- mirnaGA[-1, seq(2,ncol(mirnaGA), by=3)]
rownames(mirnaGA2) <- mirnaGA[-1, 1]
dim(mirnaGA2)

all(rownames(mirnaHiseq2) == rownames(mirnaGA2))
mirna <- cbind(mirnaHiseq2, mirnaGA2)
dim(mirna)
mirna2 <- apply(mirna, 2, as.numeric)
rownames(mirna2) <- rownames(mirna)
mirna3 <- mirna2[, grep("01A", colnames(mirna2))]
dim(mirna3)
colnames(mirna3) <- unlist(lapply(strsplit(colnames(mirna3), "\\."), function(x) x[3]))

# number of common subjects between clinical and gene, proteins and miRNA
length(Reduce(intersect, list(rownames(clinDat), colnames(genEset2), colnames(prot3), colnames(mirna3))))   ## 388

###----------------------------
#
# CpG dataset
#
###----------------------------
## Illumina 27
meth.27 <- read.delim(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0/BRCA.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt"))
meth.27.new <- meth.27[, c(1, 3, 4, 5, grep("Beta_value", as.character(as.matrix(meth.27[1,]))))]
meth.27.new2 <- meth.27.new[-1, ]
colnames(meth.27.new2) <- colnames(meth.27.new)
colnames(meth.27.new2)[1:4] <- c("Composite Element REF", "Gene_Symbol", "Chromosome", "Genomic_Coordinate")

cg.genSym <- meth.27.new2[, 1:2]
other <- apply(meth.27.new2[, -c(1:2)], 2, as.numeric)
cg.genSym.other <- cbind.data.frame(cg.genSym, other)
#saveRDS(cg.genSym.other, paste0(WhereAmI, "firebrowse/meth.27.rds"))

## Illumina 450K
## convert the large methylation datafile in smaller files using the following linux command
## split -l 100000 BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt tcga.meth-
fileNames <- list.files(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0/splitMethDatasets"),
                        full.names = TRUE)
## import the first file separately since it contains the header
meth.450.aa <- read.delim(fileNames[1])
ind <- c(1, 3, 4, 5, grep("Beta_value", as.character(as.matrix(meth.450.aa[1,]))))
meth.450.new <- meth.450.aa[, ind]
meth.450.new2 <- meth.450.new[-1, ]
colnames(meth.450.new2) <- colnames(meth.450.new)
colnames(meth.450.new2)[1:4] <- c("Composite Element REF", "Gene_Symbol", "Chromosome", "Genomic_Coordinate")
cg.genSym <- meth.450.new2[, 1:2]
other <- apply(meth.450.new2[, -c(1:2)], 2, as.numeric)
cg.genSym.other <- cbind.data.frame(cg.genSym, other)

## import the rest of the files
methList <- list()
for(i in 2 : length(fileNames)){
  meth.450 <- read.delim(fileNames[i], header = FALSE)
  meth.450 <- meth.450[, ind]
  colnames(meth.450) <- colnames(cg.genSym.other)
  methList[[i]] <- meth.450
}

methList <- rbind(cg.genSym.other, do.call(rbind, methList))
#saveRDS(methList, paste0(WhereAmI, "firebrowse/meth.450.rds"))

## combine meth.27 and meth.450
meth.27 <- readRDS(paste0(WhereAmI, "firebrowse/meth.27.rds"))
dim(meth.27); # 27578 x 347
dim(methList); # 485577 x 889

length(meth.27$`Composite Element REF`)
length(unique(as.character(meth.27$`Composite Element REF`)))
length(methList$`Composite Element REF`)
length(unique(as.character(methList$`Composite Element REF`)))

rownames(meth.27) <- as.character(meth.27$`Composite Element REF`)
rownames(methList) <- as.character(methList$`Composite Element REF`)
comSubj <- intersect(rownames(meth.27), rownames(methList))
length(comSubj);  # 25978

meth <- cbind(meth.27[comSubj, ], methList[comSubj, ])
dim(meth) # 25978  1236
saveRDS(meth, paste0(WhereAmI, "firebrowse/combined.meth.rds"))
colnames(meth)[1:5]

## methylation dataset
meth0 <- readRDS(paste0(WhereAmI, "firebrowse/combined.meth.rds"))
meth <- meth0[, -c(348:351)]  ## remove duplicate cpg ids, gene symbols, chromosome, genomic coordinates (4 columns)
dim(meth)
ncol(meth) - 4  # number of samples (removed the first 4 columns)
meth2 <- meth[, grep("01A", colnames(meth))]
dim(meth2)
colnames(meth2) <- unlist(lapply(strsplit(colnames(meth2), "\\."), function(x) x[3]))


# number of common subjects between clinical and gene, proteins and miRNA and cpg 
length(Reduce(intersect, list(rownames(clinDat), colnames(genEset2), colnames(prot3), colnames(mirna3), colnames(meth2))))   ## 387


###----------------------------
#
# PAM50 labels
#
###----------------------------
pam50 <- read.delim(paste0(WhereAmI, "BRCA.1182_pam50scores.txt"), row.names = 1)
dim(pam50)
pam50.new <- pam50[grep("01A", rownames(pam50)), ]
dim(pam50.new)
rownames(pam50.new) <- unlist(lapply(strsplit(rownames(pam50.new), "-"), function(x) x[3]))

# number of common subjects between clinical and gene, proteins and miRNA and cpg 
length(Reduce(intersect, list(rownames(clinDat), colnames(genEset2), colnames(prot3), 
                              colnames(mirna3), colnames(meth2), rownames(pam50.new))))   ## 387

comSubjects <- Reduce(intersect, list(rownames(clinDat), colnames(genEset2), colnames(prot3), 
                                      colnames(mirna3), colnames(meth2), rownames(pam50.new)))
###----------------------------
#
# 1) Training dataset (limiting dataset proteomics)
#
###----------------------------
clinTrain <- clinDat[comSubjects, ]
mrnaTrain <- genEset2[, comSubjects]
mirnaTrain <- mirna3[, comSubjects]
protTrain <- prot3[, comSubjects]
methTrain <- meth2[, comSubjects]
pam50Train <- pam50.new[comSubjects ,]
table(pam50Train$Call)

protAnnotation <- read.delim(paste0(WhereAmI, "firebrowse/gdac.broadinstitute.org_BRCA.RPPA_AnnotateWithGene.Level_3.2015110100.0.0/BRCA.antibody_annotation.txt"))
rownames(protAnnotation) <- protAnnotation[, "Composite.Element.REF"]
methAnnotation <- meth[, c("Composite Element REF", "Gene_Symbol", "Chromosome", "Genomic_Coordinate")]

###----------------------------
#
# 2) Test dataset
#
###----------------------------
comSubjects <- setdiff(Reduce(intersect, list(rownames(clinDat), colnames(genEset2),  
                                      colnames(mirna3), colnames(meth2), rownames(pam50.new))), rownames(clinTrain))
length(comSubjects)
clinTest <- clinDat[comSubjects, ]
mrnaTest <- genEset2[, comSubjects]
mirnaTest <- mirna3[, comSubjects]
methTest <- meth2[, comSubjects]
pam50Test <- pam50.new[comSubjects ,]
table(pam50Test$Call)

save(clinTrain=clinTrain, mrnaTrain=mrnaTrain, mirnaTrain=mirnaTrain, protTrain=protTrain, methTrain=methTrain, pam50Train=pam50Train,
     clinTest=clinTest, mrnaTest=mrnaTest, mirnaTest=mirnaTest, methTest=methTest, pam50Test=pam50Test,
     protAnnotation=protAnnotation, methAnnotation=methAnnotation,
     file = paste0(WhereAmI, "trainTestDatasets.RDATA"))


