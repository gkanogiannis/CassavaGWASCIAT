#!/home/agkanogiannis/bin/Rscript

#
# get_BLUPs.R
# Script part of CassavaWFMapCIAT project
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

if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(update=F)
if (!require(gdsfmt)) BiocManager::install(c("gdsfmt")) 
library(gdsfmt)
if (!require(SNPRelate)) BiocManager::install(c("SNPRelate"))
library(SNPRelate)
library(dplyr)
library(lme4)
library(lmerTest)

pheno.all <- data.frame(Genotype=character(),Year=character(),Repetition_Field=character(),stringsAsFactors=F)
trait.names <- list("DM","HCN","PPD","CT")
for(trait in trait.names){
  if(trait!="PPD") pheno.Trait <- read.table(paste0("pheno.",trait,".tsv"),sep="\t",header=T,colClasses=c("character","character","character","numeric"),stringsAsFactors=F)
  else pheno.Trait <- read.table(paste0("pheno.",trait,".tsv"),sep="\t",header=T,colClasses=c("character","character","character","character","numeric"),stringsAsFactors=F)
  pheno.Trait <- aggregate(as.formula(paste(trait, "~ Genotype+Year+Repetition_Field")),pheno.Trait, mean)
  pheno.all <- merge(pheno.all,pheno.Trait,all=T,by=c("Genotype","Year","Repetition_Field"))
}
write.table(pheno.all,file="pheno.all.tsv",sep="\t",quote=F,row.names=F)

BLUP.all <- data.frame(Genotype=character(),stringsAsFactors=F)
pdf(file = "BLUP.plots.pdf",paper="usr",width=11, height=8.5)
for(trait in trait.names){
  layout(matrix(c(1,2,3,4,5,5), 2, 3, byrow=T))
  #layout.show(5)
  pheno.Trait <- pheno.all[!is.na(pheno.all[trait]),c("Genotype","Year","Repetition_Field",trait)]
  pheno.Trait$Genotype <- as.factor(pheno.Trait$Genotype)
  pheno.Trait$Year <- as.factor(pheno.Trait$Year)
  pheno.Trait$Repetition_Field <- as.factor(pheno.Trait$Repetition_Field)
  pheno.Trait <- pheno.Trait %>% rename("Trait"=trait)
  #pheno.Trait <- unique(pheno.Trait)
  #temporary until CT data is fixed
  #pheno.Trait <- pheno.Trait %>% distinct(Genotype, Year, Repetition_Field, .keep_all = TRUE)
  
  #par(cex.axis=0.5)
  #par(cex.lab=0.75)
  hist(pheno.Trait$Trait, col="gold",xlab=trait, main=paste0(trait," histogram"))
  boxplot(pheno.Trait$Trait~pheno.Trait$Year, xlab="Year", ylab=trait, main=paste0(trait," by Year"), col="pink",)
  boxplot(pheno.Trait$Trait~pheno.Trait$Repetition_Field, xlab="Location", ylab=trait, main=paste0(trait," by Location"), col="pink")
  
  if(length(levels(pheno.Trait$Repetition_Field))>1)
    model <- lmer(Trait~
                         (1|Genotype)
                        +(Year)
                        +(Repetition_Field)
                        #+(1|Year:Repetition_Field)
                        +(1|Genotype:Year)
                        +(1|Genotype:Repetition_Field)
                        #+(1|Genotype:Repetition_Field:Year)
                      ,data=pheno.Trait)
  else
    model <- lmer(Trait~
                         (1|Genotype)
                        +(Year)
                        #+(1|Genotype:Year)
                      ,data=pheno.Trait)
  print(summary(model))
  BLUP <- ranef(model)
  BLUP.genotype <- data.frame(Genotype=rownames(BLUP$Genotype),Trait=BLUP$Genotype,stringsAsFactors=F)  
  colnames(BLUP.genotype) <- c("Genotype",trait)
  write.table(BLUP.genotype,file=paste0("BLUPs.",trait,".csv"),sep="\t",quote=F,row.names=F,col.names=T)
  BLUP.all <- merge(BLUP.all,BLUP.genotype,all=T,by=c("Genotype"))
  
  ## Creating plots with the BLUPs
  # Create a histogram with the BLUP for each genotype
  hist(BLUP.genotype[,2], col="brown",xlab=paste0(trait," BLUP"),main=paste0(trait," BLUP histogram"))
  ## Compare BLUP to line averages on a scatterplot
  lmean = tapply(pheno.Trait$Trait, pheno.Trait$Genotype, na.rm=T, mean)
  par(pty="s")
  plot(BLUP.genotype[,2], lmean, col="blue", 
       xlab=paste0(trait," BLUP"),
       ylab=paste0(trait," genotypes means"))
  
  #print(ls_means(BLUP.model, test.effs=NULL, method.grad='simple'))
}
dev.off()
write.table(BLUP.all,file=paste0("BLUPs.csv"),sep="\t",quote=F,row.names=F,col.names=T)
