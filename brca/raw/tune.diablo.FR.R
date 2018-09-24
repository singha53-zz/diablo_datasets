# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


library(mixOmics)
load("/Users/florian/Dropbox (Personal)/Manuscript/mixOmics.org:DIABLO/1-Data & Preprocessing/breastCancer/data/trainTestDatasetsNormalized.RDATA")

Y <- droplevels(pam50Train0$Call)
names(Y) <- rownames(pam50Train0)
X <- list(mRNA = mrnaTrain0, miRNA = mirnaTrain0, CpGs = methTrain0, Proteins = protTrain0)
all(names(Y) == rownames(X[[1]]))
all(names(Y) == rownames(X[[2]]))
all(names(Y) == rownames(X[[3]]))
all(names(Y) == rownames(X[[4]]))
dim(X[[1]]); dim(X[[2]]); dim(X[[3]]); dim(X[[4]]);

# check dimension
lapply(X, dim)

#outcome
summary(Y)

ncomp = 3

design <- matrix(c(0, 1, 1, 1,
1, 0, 0, 0,
1, 0, 0, 0,
1, 0, 0, 0), nrow = 4, ncol = 4, dimnames = list(names(X), names(X)))
design

set.seed(123)
test.keepX = list ("mRNA" = c(1:10,15,seq(20,100,10)), "miRNA" = c(1:10,15,seq(20,100,10)),"CpGs" = c(1:10,15,seq(20,100,10)),"Proteins" = c(1:10,15,seq(20,100,10)))
system.time(tune.block.splsda(X = X, Y = Y, ncomp = ncomp, test.keepX = test.keepX, design = design, constraint = TRUE))
list.keepX = tune.TCGA.constraint$choice.keepX
list.keepX

save(tune.TCGA.constraint,list.keepX, file = 'result-TCGA-constraint-diablo.RData')


set.seed(123)
test.keepX = list ("mRNA" = c(1:10,15,seq(20,100,10)), "miRNA" = c(1:10,15,seq(20,100,10)),"CpGs" = c(1:10,15,seq(20,100,10)),"Proteins" = c(1:10,15,seq(20,100,10)))
tune.TCGA = tune.block.splsda(X = X, Y = Y, ncomp = ncomp, test.keepX = test.keepX, design = design, constraint = FALSE)
list.keepX = tune.TCGA$choice.keepX
list.keepX

save(tune.TCGA,list.keepX, file = 'result-TCGA-diablo.RData')




