#good tut https://ycl6.github.io/16S-Demo/4_picrust2_tutorial.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")
library("data.table")   # Also requires R.utils to read gz and bz2 files
library("phyloseq")
#BiocManager::install("ALDEx2")
library("ALDEx2")
library("dplyr")
library("R.utils")
setwd("/projectnb/davieslab/jfifer/Fouha/PICRUST2/all.rare.trim/")

p2EC = as.data.frame(fread("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread("EC_pathways_out/path_abun_unstrat.tsv.gz"))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)

#load the phyloseq objects 
load("16S_PostTrim&Rarefy.RData")



###comps within  edge center and sed###
p2EC_edge = p2EC[,sample_names(onlyedgeps.rare.trim)]
p2EC_center = p2EC[,sample_names(onlycenterps.rare.trim)]
p2EC_sed= p2EC[,sample_names(onlysedps.rare.trim)]

p2KO_edge = p2KO[,sample_names(onlyedgeps.rare.trim)]
p2KO_center = p2KO[,sample_names(onlycenterps.rare.trim)]
p2KO_sed= p2KO[,sample_names(onlysedps.rare.trim)]

p2PW_edge = p2PW[,sample_names(onlyedgeps.rare.trim)]
p2PW_center = p2PW[,sample_names(onlycenterps.rare.trim)]
p2PW_sed= p2PW[,sample_names(onlysedps.rare.trim)]


##comps between edge, center and sed
p2EC_all = p2EC[,sample_names(ps.all.rare.trim)]
p2KO_all = p2KO[,sample_names(ps.all.rare.trim)]
p2PW_all= p2PW[,sample_names(ps.all.rare.trim)]



##First EC##
#If the test is ‘glm’, then effect should be FALSE. The ‘glm’ option evaluates the data as a one-way ANOVA using the glm and Kruskal-Wallace test
#If the test is ‘t’, then effect should be set to TRUE. The ‘t’ option evaluates the data as a two-factor experiment using both the Welch’s t and the Wilcoxon rank tests
#All tests include a Benjamini-Hochberg correction of the raw P values

test=sample_data(ps.all.rare.trim)$position
g=as.data.frame(test)

#make model matrices changing the reference group
mmC=model.matrix(~test,g)
colnames(mmC)
mmE=cbind(mmC[,2],mmC[,1],mmC[,3])
colnames(mmE)=c("(Intercept)", "testCenter","testOther")
#center as ref
#useMC refers to multicore not montecarlo
#first try, denom=iqlr gives subscript out of bounds not sure why
clr.EC.C=aldex.clr(p2EC_all,mmC, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C.glm_EC=aldex.glm(clr.EC.C, verbose=TRUE)

#edge as ref
clr.EC.E=aldex.clr(p2EC_all,mmE, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E.glm_EC=aldex.glm(clr.EC.E, verbose=TRUE)

##Next KO##
#center as ref
clr.KO.C=aldex.clr(p2KO_all,mmC, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C.glm_KO=aldex.glm(clr.KO.C, verbose=TRUE)

#edge as ref
clr.KO.E=aldex.clr(p2KO_all,mmE, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E.glm_KO=aldex.glm(clr.KO.E, verbose=TRUE)

##Next pathways##
#center as ref
clr.PW.C=aldex.clr(p2PW_all,mmC, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C.glm_PW=aldex.glm(clr.PW.C, verbose=TRUE)




#######################################################################################################


#edge as ref

clr.PW.E=aldex.clr(p2PW_all,mmE, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E.glm_PW=aldex.glm(clr.PW.E, verbose=TRUE)

#save.image("aldexoutput_EdgeCenSed.RData")

############## Now within edge comparing sites ############
test=sample_data(onlyedgeps.rare.trim)$site
g=as.data.frame(test)

#make model matrices changing the reference group
mm2=model.matrix(~test,g)
colnames(mm2)
mm3=cbind(mm2[,2],mm2[,1],mm2[,3])
colnames(mm3)=c("(Intercept)", "test2","test4")
#EC
clr.EC.E2=aldex.clr(p2EC_edge,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E2.glm_EC=aldex.glm(clr.EC.E2, verbose=TRUE)

clr.EC.E3=aldex.clr(p2EC_edge,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E3.glm_EC=aldex.glm(clr.EC.E3, verbose=TRUE)

#none sig 
#KO
clr.KO.E2=aldex.clr(p2KO_edge,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E2.glm_KO=aldex.glm(clr.KO.E2, verbose=TRUE)

clr.KO.E3=aldex.clr(p2KO_edge,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E3.glm_KO=aldex.glm(clr.KO.E3, verbose=TRUE)
#none sig
#pathways
clr.PW.E2=aldex.clr(p2PW_edge,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E2.glm_PW=aldex.glm(clr.PW.E2, verbose=TRUE)

clr.PW.E3=aldex.clr(p2PW_edge,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
E3.glm_PW=aldex.glm(clr.PW.E3, verbose=TRUE)
#none significant 

#save.image("aldexoutput_Edge.RData")
############## Now within center comparing sites ############
test=sample_data(onlycenterps.rare.trim)$site
g=as.data.frame(test)

#make model matrices changing the reference group
mm2=model.matrix(~test,g)
colnames(mm2)
mm3=cbind(mm2[,2],mm2[,1],mm2[,3])
colnames(mm3)=c("(Intercept)", "test2","test4")
#EC
clr.EC.C2=aldex.clr(p2EC_center,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C2.glm_EC=aldex.glm(clr.EC.C2, verbose=TRUE)

clr.EC.C3=aldex.clr(p2EC_center,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C3.glm_EC=aldex.glm(clr.EC.C3, verbose=TRUE)
#KO
clr.KO.C2=aldex.clr(p2KO_center,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C2.glm_KO=aldex.glm(clr.KO.C2, verbose=TRUE)

clr.KO.C3=aldex.clr(p2KO_center,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C3.glm_KO=aldex.glm(clr.KO.C3, verbose=TRUE)
#pathways
clr.PW.C2=aldex.clr(p2PW_center,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C2.glm_PW=aldex.glm(clr.PW.C2, verbose=TRUE)

clr.PW.C3=aldex.clr(p2PW_center,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
C3.glm_PW=aldex.glm(clr.PW.C3, verbose=TRUE)

#save.image("aldexoutput_Center.RData")

################ Now within sed comparing sites #############

test=sample_data(onlysedps.rare.trim)$site
g=as.data.frame(test)

#make model matrices changing the reference group
mm2=model.matrix(~test,g)
colnames(mm2)
mm3=cbind(mm2[,2],mm2[,1],mm2[,3])
colnames(mm3)=c("(Intercept)", "test2","test4")
#EC
clr.EC.S2=aldex.clr(p2EC_sed,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
S2.glm_EC=aldex.glm(clr.EC.S2, verbose=TRUE)

#Sigs! different between 2 and 3!
subset(S2.glm_EC,S2.glm_EC[,14] < .1)  

clr.EC.S3=aldex.clr(p2EC_sed,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
S3.glm_EC=aldex.glm(clr.EC.S3, verbose=TRUE)

subset(S3.glm_EC,S3.glm_EC[,15] < .1)  


#KO
clr.KO.S2=aldex.clr(p2KO_sed,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
S2.glm_KO=aldex.glm(clr.KO.S2, verbose=TRUE)

subset(S2.glm_KO,S2.glm_KO[,15] < .1)  


clr.KO.S3=aldex.clr(p2KO_sed,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
S3.glm_KO=aldex.glm(clr.KO.S3, verbose=TRUE)

subset(S3.glm_KO,S3.glm_KO[,14] < .1)  
#none sig 

#pathways
clr.PW.S2=aldex.clr(p2PW_sed,mm2, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
S2.glm_PW=aldex.glm(clr.PW.S2, verbose=TRUE)

clr.PW.S3=aldex.clr(p2PW_sed,mm3, mc.samples = 500, denom="all", verbose=FALSE, useMC=FALSE)
S3.glm_PW=aldex.glm(clr.PW.S3, verbose=TRUE)

#none sig

save.image("aldexoutput_all.RData")






# fileConn<-file("aldex.sh")
# writeLines(c("#!/bin/bash -l","#$ -V","#$ -cwd","#$ -N aldex","#$ -l h_rt=12:00:00",
#              "#$ -M james.e.fifer@gmail.com","#$ -m be","module load R","Rscript ./aldex2.R"), fileConn)
# close(fileConn)



