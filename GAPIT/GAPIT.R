#!/home/agkanogiannis/bin/Rscript

#
# GAPIT.R
# Script part of CassavaGWASCIAT project
#
# Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

if (!requireNamespace("BiocManager",quietly=T)) install.packages("BiocManager")
BiocManager::install(update=F)
library(dplyr)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler)
library(scatterplot3d)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

file.relatedness2 <- "Cassava_RAD_CIAT_GWAS.snps.relatedness2.txt"

# Convert relatedness2 from vcftools to kinship csv, if it has not been already
if(!file.exists("GAPIT.Kin.relatedness2.csv") & file.exists(file.relatedness2)){
  relatedness2 <- read.table(file.relatedness2, head=T,colClasses=c(NA,NA,"NULL","NULL","NULL","NULL",NA))
  colnames(relatedness2) <- c("INDV1","INDV2","kinhsip")
  relatedness2.samples.id <- unique(relatedness2$INDV1)
  relatedness2.kinship <- data.frame(matrix(ncol=length(relatedness2.samples.id)+1, nrow=length(relatedness2.samples.id)))
  colnames(relatedness2.kinship) <- c("taxa",relatedness2.samples.id)
  class(relatedness2.kinship[,1]) <- "character"
  for(i in c(2:ncol(relatedness2.kinship))) relatedness2.kinship[,i] <- as.numeric(relatedness2.kinship[,i])
  relatedness2.kinship[,1] <- relatedness2.samples.id
  for (s1 in relatedness2.samples.id){
    for (s2 in relatedness2.samples.id){
      kinship_s1_s2 <- relatedness2[relatedness2$INDV1==s1 & relatedness2$INDV2==s2,"kinhsip"]
      relatedness2.kinship[relatedness2.kinship$taxa==s1,s2] <- kinship_s1_s2
    }
  }
  write.table(relatedness2.kinship,file="GAPIT.Kin.relatedness2.csv",sep=",", quote=F,row.names=F,col.names=F)
}
  
# Get kinship
#myKI <- read.csv("GAPIT.Kin.relatedness2.csv", header=F)  
#myKI <- read.csv("GAPIT.Kin.VanRaden.csv", header=F)  

# Get structure coeffs and remove last (to make them independend)
myCV <- read.table("Cassava_RAD_CIAT_GWAS.snps.qmatrix.txt", head=F)
colnames(myCV) <- c("taxa","G1","G2","G3","G4","G5","G6","G7")
#myCV <- myCV[,-ncol(myCV)]

# Get phenotypes
myY <- read.table("BLUPs.csv", head=T,sep="\t",stringsAsFactors=F)
myY <- myY %>% rename("taxa"="Genotype")

# Get Genotypes from hapmap.
#initial = read.table("Cassava_RAD_CIAT_GWAS.snps.hmp.txt",nrows=100,comment.char="")
#classes = sapply(initial, class)
myG <- read.table("Cassava_RAD_CIAT_GWAS.snps.hmp.txt",header=F, comment.char="")
myG[,3] <- as.numeric(myG[,3])
myG[,4] <- as.numeric(myG[,4])
# Convert hapmap to numerical
#convertGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
#myGD <- convertGAPIT$GD
#myGM <- convertGAPIT$GM
# Get numerical genotypes GD and GM
#myGD <- read.table("Cassava_RAD_CIAT_GWAS.snps.numerical.txt", head=T)
#myGM <- read.table("Cassava_RAD_CIAT_GWAS.snps.map.txt", head=T)

# Find common samples.id in Y and G
samples.Y.id <- unique(myY$taxa)
samples.G.id <- unique(as.character(myG[1,c(12:ncol(myG))]))
samples.id <- samples.Y.id[samples.Y.id %in% samples.G.id]

# Subset Y and CV
myY.sub <- subset(myY, taxa %in% samples.id)
myCV.sub <- subset(myCV, taxa %in% samples.id)


# Run GAPIT
myGAPIT <- GAPIT(
  Y=myY.sub[,1:5],
  #GD=myGD,
  #GM=myGM,
  G=myG,
  #KI=myKI,
  CV=myCV.sub,
  PCA.total=0,
  #PCA.total=6,
  kinship.algorithm="VanRaden",
  model=c("MLM","MLMM","FarmCPU","SUPER"),
  #model=c("GLM", "MLM", "CMLM", "SUPER","MLMM", "FarmCPU", "Blink"),
  Multiple_analysis=T,
  Model.selection=F,
  Major.allele.zero=T
  #disable compression
  #group.from=length(samples.id),
  #group.to=length(samples.id),
  #group.by=1
)
