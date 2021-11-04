install.packages("indicspecies")
library(indicspecies)
setwd("~/BOSTON/Davies/Fouha Microbiome/dada2")

load("~/BOSTON/Davies/Fouha Microbiome/dada2/16S_PostTrim&Rarefy.RData")
#Figure out most abundant ASVs for each group
#edge
edgeASVsums=data.frame(taxa_sums(onlyedgeps.rare.trim))
edgeASVsums$percent=edgeASVsums$taxa_sums.onlyedgeps.rare.trim./sum(edgeASVsums$taxa_sums.onlyedgeps.rare.trim.)*100
#highest ASV is like 3% of total

#center
centerASVsums=data.frame(taxa_sums(onlycenterps.rare.trim))
centerASVsums$percent=centerASVsums$taxa_sums.onlycenterps.rare.trim./sum(centerASVsums$taxa_sums.onlycenterps.rare.trim.)*100

allTaxa <- taxa_names(onlycenterps.rare.trim)
prunedTaxa <- prune_taxa(allTaxa,onlycenterps.rare.trim )
#ASV 7 10.3% Kingdom Bacteria
#ASV 33 8.2% Kingdom Bacteria
#ASV 10 8.2% Xenococcaceae
#ASV 2 8% Endozoicomonas

#Family level
#glom doesnt work if there are families that have different listed kingdoms for example, so intead of dealing with that
#can just look at mean relative abundances for each Family 



#quick check with relative abundances
all.rare.rel <- transform_sample_counts(ps.all.rare.trim, function(x) x / sum(x))
melt.all.rare.rel<- psmelt(all.rare.rel)
#plot_bar(all.rare.rel) all 1.0 sanity check
TEST=melt.all.rare.rel %>% group_by(Sample) %>% top_n(1, Abundance)
TEST2=filter(TEST,position!="Other")


##grouping by 
#grouping by fam
library(dplyr)
#MAKE SURE YOU DETACH BEFORE RUNNING MELT
detach(package:Rmisc)
detach(package:plyr)
huh=melt.all.rare.rel %>%
  group_by(Family,Order, Sample, Colony, position,site) %>%
  summarize(Abundance = sum(Abundance))

#plotting top 5 families
library(RColorBrewer)
n <- 59
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector
edgeonly=filter(huh,position=="Edge")
#CHECK MEANs by family
edgeonlyfammeans=aggregate(edgeonly[,7 ], list(edgeonly$Family), mean)
#save for supplementary file
#write.csv(edgeonlyfammeans,file="EdgeFamilyMeanAbund.csv")
top10.edge=edgeonly %>% group_by(Sample) %>% top_n(5, Abundance)
top10.edge <- top10.edge[order(top10.edge$site),]
top10.edge$site<-factor(top10.edge$site,levels=unique(top10.edge$site))
top10.edge$Sample<-factor(top10.edge$Sample,levels=unique(top10.edge$Sample))
top10.edge$Family<-factor(top10.edge$Family,levels=unique(top10.edge$Family))

#center
centeronly=filter(huh,position=="Center")
centeronlyfammeans=aggregate(centeronly[,7 ], list(centeronly$Family), mean)
write.csv(centeronlyfammeans,file="CenterFamilyMeanAbund.csv")

top10.center=centeronly %>% group_by(Sample) %>% top_n(3, Abundance)
top10.center <- top10.center[order(top10.center$site),]
#top10.center$site<-factor(top10.center$site,levels=unique(top10.center$site))
top10.center$Sample<-factor(top10.center$Sample,levels=unique(top10.center$Sample))
top10.center$Family<-factor(top10.center$Family,levels=unique(top10.center$Family))

dd <- union(top10.edge$Family,top10.center$Family)
library(RColorBrewer)
n <- 59
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector
length(col_vector)<-length(dd)
colours <- setNames(as.character(col_vector), 
                    dd)

ok.edge=as.character(unique(top10.edge$Family))
ok.center=as.character(unique(top10.center$Family))


col.edge=colours[ok.edge]
col.cen=colours[ok.center]


p <- ggplot(data=top10.edge, aes(x=Sample, y=Abundance,fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=col.edge,na.value="black")+
  theme_bw()+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #theme_cowplot()+
  
  guides(fill=guide_legend(nrow=5))+xlab(NULL)


p1 <- ggplot(data=top10.center, aes(x=Sample, y=Abundance,fill=Family))
p1 + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=col.cen,na.value="black")+
  theme_bw()+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #theme_cowplot()+

  guides(fill=guide_legend(nrow=5))+xlab(NULL)

#sed

sedonly=filter(huh,position=="Other")
sedonlyfammeans=aggregate(sedonly[,7 ], list(sedonly$Family), mean)
write.csv(sedonlyfammeans,file="SedFamilyMeanAbund.csv")


top10.sed=sedonly %>% group_by(Sample) %>% top_n(5, Abundance)
top10.sed <- top10.sed[order(top10.sed$site),]
#top10.sed$site<-factor(top10.sed$site,levels=unique(top10.sed$site))
top10.sed$Sample<-factor(top10.sed$Sample,levels=unique(top10.sed$Sample))

p <- ggplot(data=top10.sed, aes(x=Sample, y=Abundance, fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=col_vector,na.value="black")+
  theme_bw()+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  
  guides(fill=guide_legend(nrow=5))+xlab(NULL)


#indicator start
#normal 
#create community data matrix from phyloseq object

cdm1=as(otu_table(ps.all.trim), "matrix")
if(taxa_are_rows(ps.all.trim)){cdm1 <- t(cdm1)}
indval = multipatt(cdm1, new.new.new.meta$site, control = how(nperm=999))
summary(indval)


#edge only 
cdm_edge=as(otu_table(onlyedgepsall.trim), "matrix")
if(taxa_are_rows(onlyedgepsall.trim)){cdm_edge <- t(cdm_edge)}
edge.meta=subset(new.new.new.meta, position=="Edge")
indvale = multipatt(cdm_edge, edge.meta$site, control = how(nperm=999))
summary(indvale)

#center only
cdm_center=as(otu_table(onlycenterpsall.trim), "matrix")
if(taxa_are_rows(onlycenterpsall.trim)){cdm_center <- t(cdm_center)}
center.meta=subset(new.new.new.meta, position=="Center")
indvalc = multipatt(cdm_center, center.meta$site, control = how(nperm=999))
summary(indvalc)

#sed only
cdm_sed=as(otu_table(onlysedpsall.trim), "matrix")
if(taxa_are_rows(onlysedpsall.trim)){cdm_sed <- t(cdm_sed)}
sed.meta=subset(new.new.new.meta, position=="Other")
indvals = multipatt(cdm_sed, sed.meta$site, control = how(nperm=999))
summary(indvals)

##Checking taxa
edgetax=taxa_names(onlyedgepsall.trim)
edgetaxind=row.names(subset(indvale[["sign"]],p.value <0.05))
edgeind.ps=prune_taxa(edgetaxind, onlyedgepsall.trim)
plot_bar(edgeind.ps,x="site",fill="Family")+
  facet_wrap(~Family,scales="free")

###rarefied samples ###

#create community data matrix from phyloseq object
library(indicspecies)
cdm.all.rare=as(otu_table(ps.all.rare.trim), "matrix")
if(taxa_are_rows(ps.all.rare.trim)){cdm.all.rare <- t(cdm.all.rare)}
indval.rare = multipatt(cdm.all.rare, new.new.new.meta$site, control = how(nperm=999))
summary(indval.rare)

#edge only 
cdm_edge=as(otu_table(onlyedgeps.rare.trim), "matrix")
if(taxa_are_rows(onlyedgeps.rare.trim)){cdm_edge <- t(cdm_edge)}
edge.meta=subset(new.new.new.meta, position=="Edge")
indvale.rare = multipatt(cdm_edge, edge.meta$site, control = how(nperm=999))
summary(indvale.rare)

#center only
cdm_center=as(otu_table(onlycenterps.rare.trim), "matrix")
if(taxa_are_rows(onlycenterps.rare.trim)){cdm_center <- t(cdm_center)}
center.meta=subset(new.new.new.meta, position=="Center")
indvalc.rare = multipatt(cdm_center, center.meta$site, control = how(nperm=999))
summary(indvalc.rare)

#sed only
cdm_sed=as(otu_table(onlysedps.rare.trim), "matrix")
if(taxa_are_rows(onlysedps.rare.trim)){cdm_sed <- t(cdm_sed)}
sed.meta=subset(new.new.new.meta, position=="Other")
indvals.rare = multipatt(cdm_sed, sed.meta$site, control = how(nperm=999))
summary(indvals.rare)

#EdgeOnly Genus glom
ps.edge.gen <- tax_glom(onlyedgeps.rare.trim, "Genus")
cdm_edge_gen=as(otu_table(ps.edge.gen), "matrix")
if(taxa_are_rows(ps.edge.gen)){cdm_edge_gen <- t(cdm_edge_gen)}
edge.meta=subset(new.new.new.meta, position=="Edge")
indvale.gen.rare = multipatt(cdm_edge_gen, edge.meta$site, control = how(nperm=999))
summary(indvale.gen.rare)

#multitest correction, actually according to  Caceres et al (2010) multiple test correction is not necessary. "When reporting the results of indicator species analysis for several species, users should be aware of multiple testing issues. Let us say that we conduct indicator species analysis with ?? 0.05 on 200 species without correcting for multiple testing and obtain 10 signifi cant results. If we say that there are 10 indicator species, this 'experiment-wise' statement will probably be wrong, because under the null hypothesis of no association the expected number of signifi cant results with 200 species is ?? 200 10. Corrections for multiple testing are advisable in this case. Th ese procedures (e.g. see the 'p.adjust' function of R) modify the p-values in order to keep the probability of fi nding, among all the statistical tests, at least one signifi cant result at the chosen signifi cance level. After the correction, we can report the number of signifi - cant indicators more safely. If, however, we are interested in reporting that a given species is an indicator, we do not need any correction because we are not making any experimentwise statement. If the test is exact, the probability of type I error will be equal to ?? in that case."
#extracting sigs (no mulitple test correction)

#specifying sigs 
indval.rare.out <- data.frame(indval.rare[["sign"]])
indval.rare.out <- indval.rare.out[complete.cases(indval.rare.out),]
indval.rare.out.sig <- subset(indval.rare.out,p.value <= 0.05)

goodtaxa <- rownames(indval.rare.out.sig)
allTaxa <- taxa_names(ps.all.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
all.rare.rel <- transform_sample_counts(ps.all.rare.trim, function(x) x / sum(x))
ind.rare.rel <- prune_taxa(allTaxa, all.rare.rel)
plot_bar(ind.rare.rel,x="site",fill="Genus")+
  facet_wrap(~Class,scales="free")

#edge
indvale.rare.out <- data.frame(indvale.rare[["sign"]])
indvale.rare.out <- indvale.rare.out[complete.cases(indvale.rare.out),]
indvale.rare.out.sig <- subset(indvale.rare.out,p.value <= 0.05)

goodtaxa <- rownames(indvale.rare.out.sig)
allTaxa <- taxa_names(onlyedgeps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
edge.rare.rel <- transform_sample_counts(onlyedgeps.rare.trim, function(x) x / sum(x))
ind.edge.rare.rel <- prune_taxa(allTaxa, edge.rare.rel)
plot_bar(ind.edge.rare.rel,x="site",fill="Species")+
  facet_wrap(~Genus,scales="free")

#bubble plot
library(Rmisc)
melt.ind.edge.rare.rel<- psmelt(ind.edge.rare.rel)
melt.ind.edge.rare.rel.sum <- summarySE(melt.ind.edge.rare.rel,measurevar="Abundance",groupvars = c("OTU","Genus","site","Family","Order"))

melt.indic.edge= melt(indvale.rare.out.sig)
what=data.frame(c(indvale.rare.out.sig$s.2,indvale.rare.out.sig$s.3,indvale.rare.out.sig$s.4))
names=rep(rownames(indvale.rare.out.sig),3)
what$rownames=names
colnames(what)=c("index","rownames")
newdata <- what[order(what$rownames),]

melt.indic.edge2=cbind(melt.ind.edge.rare.rel.sum,newdata)
melt.indic.edge2$Genus=as.character(melt.indic.edge2$Genus)
#Add family to genus if missing genus 
melt.indic.edge2[58:60,2]="Planococcaceae"
melt.indic.edge2[10:12,2]="Methyloligellaceae"
melt.indic.edge2[7:9,2]="Kiloniellaceae"
melt.indic.edge2[16:18,2]="Spirochaeta"
melt.indic.edge2[28:30,2]="Cyanobium"
melt.indic.edge2$Genus_ASV=paste(melt.indic.edge2$Genus,melt.indic.edge2$OTU)


ggplot(melt.indic.edge2, aes(x = site, y = Genus_ASV)) + 
  geom_point(aes(size = Abundance, fill = site,shape=as.factor(index)), alpha = 0.75)+
 # facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
 scale_fill_manual(values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+
  scale_shape_manual(values=c(1, 21))+
  guides(fill=FALSE,shape=FALSE)+ggtitle("Edge Only")+ylab(NULL)+xlab(NULL)


#Genera that are doing interesting things:
#Acanthopleuribacter- up in 4, low in 2/3, cant find anything on it
#catenococcus, 2,3,4 decreasing gradient, indicator in 2 only
#Few marine studies have recorded Catenococcus spp. [78,79,80]. However, this member of the Vibrionaceae family is described as a pathogen in the seaweed Kappaphycus alvarezii in which infection by Catenococcus thiocyli causes bleaching [79]. Bacterial extracts from sponge Stylotella sp. of Proteobacteria closely related to Catenococcus thiocycli, showed toxicity activity [80]. As nearly all dying clams harbored a typical microbiome mainly composed of Vibrionaceae (Catenococcus spp. or an unclassified genus), we hypothesize that the combined effect of secondary metabolites from Acropora corals and increased water temperature may have weakened the clams' defenses against Vibrio infection.
#^

#center
indvalc.rare.out <- data.frame(indvalc.rare[["sign"]])
indvalc.rare.out <- indvalc.rare.out[complete.cases(indvalc.rare.out),]
indvalc.rare.out.sig <- subset(indvalc.rare.out,p.value <= 0.05)

goodtaxa <- rownames(indvalc.rare.out.sig)
allTaxa <- taxa_names(onlycenterps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
center.rare.rel <- transform_sample_counts(onlycenterps.rare.trim, function(x) x / sum(x))
ind.center.rare.rel <- prune_taxa(allTaxa, center.rare.rel)
plot_bar(ind.center.rare.rel,x="site",fill="Genus")+
  facet_wrap(~Class,scales="free")

#bubble plot
library(Rmisc)
melt.ind.center.rare.rel<- psmelt(ind.center.rare.rel)
melt.ind.center.rare.rel.sum <- summarySE(melt.ind.center.rare.rel,measurevar="Abundance",groupvars = c("OTU","Genus","site","Family","Order"))


what=data.frame(c(indvalc.rare.out.sig$s.2,indvalc.rare.out.sig$s.3,indvalc.rare.out.sig$s.4))
names=rep(rownames(indvalc.rare.out.sig),3)
what$rownames=names
colnames(what)=c("index","rownames")
newdata <- what[order(what$rownames),]

melt.indic.center2=cbind(melt.ind.center.rare.rel.sum,newdata)
melt.indic.center2$Genus=as.character(melt.indic.center2$Genus)
#Add family to genus if missing genus 
melt.indic.center2[13:15,2]="Planococcaceae"
melt.indic.center2[22:24,2]="Francisellaceae"
melt.indic.center2[25:27,2]="Planococcaceae"
melt.indic.center2[28:30,2]="Burkholderiaceae"

melt.indic.center2$Genus_ASV=paste(melt.indic.center2$Genus,melt.indic.center2$OTU)


ggplot(melt.indic.center2, aes(x = site, y = Genus_ASV)) + 
  geom_point(aes(size = Abundance, fill = site,shape=as.factor(index)), alpha = 0.75)+
  # facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
  scale_fill_manual(values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+
  scale_shape_manual(values=c(1, 21))+
  guides(fill=FALSE,shape=FALSE)+ggtitle("Center Only")+ylab(NULL)+xlab(NULL)

#sed
indvals.rare.out <- data.frame(indvals.rare[["sign"]])
indvals.rare.out <- indvals.rare.out[complete.cases(indvals.rare.out),]
indvals.rare.out.sig <- subset(indvals.rare.out,p.value <= 0.05)

goodtaxa <- rownames(indvals.rare.out.sig)
allTaxa <- taxa_names(onlysedps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
sed.rare.rel <- transform_sample_counts(onlysedps.rare.trim, function(x) x / sum(x))
ind.sed.rare.rel <- prune_taxa(allTaxa, sed.rare.rel)
plot_bar(ind.sed.rare.rel,x="site",fill="Genus")+
  facet_wrap(~Genus,scales="free")
###
#Check overlap between the above. 
#center/edge indics
intersect(rownames(indvale.rare.out.sig), rownames(indvalc.rare.out.sig))
#728! indic species in 1 in both!!!
#ASV_728  "Bacteria" "Firmicutes"     "Bacilli"             "Bacillales"              "Planococcaceae"
#https://www.nature.com/articles/ismej201485-Planococcaceae%20and%20Aeromonadaceae%20were%20in%20the%20top%2015%20families%20in%20HH.
#Pyrene polluted bacterial communities switch from Planococcaceae to Pseudomonadaceae  https://www.sciencedirect.com/science/article/pii/S001393512030027X?casa_token=IeUNcUtPjpQAAAAA:x-9Ba0EMkcaUQ8L1-R8O1pZYeUXsNYBa06S4HQjS9aZqsEc_-Jzm8Goh2HAc3iiKjc_mi5r0xw
#Planococcaceae,well represented (2.5%) among the cultured isolates, were absent from the environmental samples file:///C:/Users/james/Downloads/microorganisms-09-00619.pdf
#Way too many different genera with different functions to make any broad strokes assignment of function- https://link.springer.com/referenceworkentry/10.1007%2F978-3-642-30120-9_351
 
goodtaxa = intersect(rownames(indvale.rare.out.sig), rownames(indvalc.rare.out.sig))
allTaxa <- taxa_names(onlyedgeps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
edge.rare.rel <- transform_sample_counts(onlyedgeps.rare.trim, function(x) x / sum(x))
ind.edge.rare.rel <- prune_taxa(allTaxa, edge.rare.rel)
library(cowplot)
plot_bar(ind.edge.rare.rel,x="site",fill="site")+
theme_cowplot()+theme(text = element_text(size = 20),legend.position = "none")+
    scale_fill_manual(name="site",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+xlab(NULL)+
annotate(geom="text", x=2.5, y=0.06, label="ASV_728",size=6)+
  annotate(geom="text", x=2.7, y=0.055, label="Planococcaceae",size=6)
  
allTaxa <- taxa_names(onlycenterps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
center.rare.rel <- transform_sample_counts(onlycenterps.rare.trim, function(x) x / sum(x))
ind.center.rare.rel <- prune_taxa(allTaxa, center.rare.rel)


plot_bar(ind.center.rare.rel,x="site",fill="site")+
  theme_cowplot()+theme(text = element_text(size = 20),legend.position = "none")+
  scale_fill_manual(name="site",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+xlab(NULL)+
  annotate(geom="text", x=2.5, y=0.16, label="ASV_728",size=6)+
  annotate(geom="text", x=2.7, y=0.15, label="Planococcaceae",size=6)

#center/seds indics
intersect(rownames(indvalc.rare.out.sig), rownames(indvals.rare.out.sig))
#0
#edge/seds indics
intersect(rownames(indvale.rare.out.sig), rownames(indvals.rare.out.sig))
#ASV_259
#YESSSS BOTH @ SITE 1!!!!
#"Bacteria" "Proteobacteria" "Deltaproteobacteria" "Bdellovibrionales"       "Bacteriovoracaceae" "Peredibacter"
#The Bdellovibrio and like organisms are a family of unique
#prokaryotes that prey on many Gram-negative bacteria.They are ubiquitous in soil, water, and sewage [8, 21, 22,28] and have also been recovered from animals [11] and
#humans (Edao 2000, Ph.D. thesis, University of Leipzig).
#The group consists of four genera, Bdellovibrio, Bacteriovorax [2], Peredibacter (formerly Bdellovibrio starrii)
#[5], and a newly established genus, Bacteriolyticum
#(formerly Bdellovibrio stolpii) [2, 18]. As the result of a
#recent re-classification, the Bacteriovorax now consist exclusively of the salt-water predators [18] and are referred to as such in this article
#https://link.springer.com/content/pdf/10.1007/s00248-009-9499-7.pdf
goodtaxa=intersect(rownames(indvale.rare.out.sig), rownames(indvals.rare.out.sig))

allTaxa <- taxa_names(onlyedgeps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
edge.rare.rel <- transform_sample_counts(onlyedgeps.rare.trim, function(x) x / sum(x))
ind.edge.rare.rel <- prune_taxa(allTaxa, edge.rare.rel)
library(cowplot)
plot_bar(ind.edge.rare.rel,x="site",fill="site")+
  theme_cowplot()+theme(text = element_text(size = 20),legend.position = "none")+
  scale_fill_manual(name="site",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+xlab(NULL)+
  annotate(geom="text", x=2.5, y=0.009, label="ASV_259",size=6)+
  annotate(geom="text", x=2.59, y=0.008, label="Peredibacter",size=6)

allTaxa <- taxa_names(onlysedps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
sed.rare.rel <- transform_sample_counts(onlysedps.rare.trim, function(x) x / sum(x))
ind.sed.rare.rel <- prune_taxa(allTaxa, sed.rare.rel)


plot_bar(ind.sed.rare.rel,x="site",fill="site")+
  theme_cowplot()+theme(text = element_text(size = 20),legend.position = "none")+
  scale_fill_manual(name="site",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+xlab(NULL)+
  annotate(geom="text", x=2.5, y=0.027, label="ASV_259",size=6)+
  annotate(geom="text", x=2.59, y=0.025, label="Peredibacter",size=6)

save.image(file="16S_indic_analyses.RData")
##Edge vs Center vs Sed

cdm.all.rare=as(otu_table(ps.all.rare.trim), "matrix")
if(taxa_are_rows(ps.all.rare.trim)){cdm.all.rare <- t(cdm.all.rare)}
indval.rare = multipatt(cdm.all.rare, new.new.new.meta$position, control = how(nperm=999))
summary(indval.rare)

indval.rare.out <- data.frame(indval.rare[["sign"]])
indval.rare.out <- indval.rare.out[complete.cases(indval.rare.out),]
indval.rare.out.sig <- subset(indval.rare.out,p.value <= 0.05)

goodtaxa <- rownames(indval.rare.out.sig)
allTaxa <- taxa_names(ps.all.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
all.rare.rel <- transform_sample_counts(ps.all.rare.trim, function(x) x / sum(x))
ind.rare.rel <- prune_taxa(allTaxa, all.rare.rel)
plot_bar(ind.rare.rel,x="position",fill="Genus")+
  facet_wrap(~Class,scales="free")+theme(legend.position = "none")

#center + edge is boring its just coral taxa...but i guess worth putting in, center only& edge only is intersting and edge + sediment is intersting. 
#center+edge
centerplusedge=subset(indval.rare.out.sig,index==4)
goodtaxa <- rownames(centerplusedge)
allTaxa <- taxa_names(ps.all.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
all.rare.rel <- transform_sample_counts(ps.all.rare.trim, function(x) x / sum(x))
ind.rare.rel <- prune_taxa(allTaxa, all.rare.rel)
#sample_data(all.rare.rel)$position=c("Center","Edge","Sed")
library(Rmisc)
melt.ind.centerplusedge.rare.rel<- psmelt(ind.rare.rel)
melt.ind.centerplusedge.rare.rel.sum <- summarySE(melt.ind.centerplusedge.rare.rel,measurevar="Abundance",groupvars = c("OTU","Genus","position","Family","Order"))

melt.ind.centerplusedge.rare.rel.sum$Genus=as.character(melt.ind.centerplusedge.rare.rel.sum$Genus)

melt.ind.centerplusedge.rare.rel.sum[1:3,2]="Xenococcaceae"
melt.ind.centerplusedge.rare.rel.sum[4:6,2]="Francisellaceae"
melt.ind.centerplusedge.rare.rel.sum[19:21,2]="Rickettsiales"
melt.ind.centerplusedge.rare.rel.sum[28:30,2]="Prochlorococcus"
melt.ind.centerplusedge.rare.rel.sum[34:36,2]="Burkholderia"

melt.ind.centerplusedge.rare.rel.sum$Genus_ASV=paste(melt.ind.centerplusedge.rare.rel.sum$Genus,melt.ind.centerplusedge.rare.rel.sum$OTU)

ggplot(melt.ind.centerplusedge.rare.rel.sum, aes(x = position, y = Genus_ASV)) + 
  geom_point(aes(size = Abundance, colour = as.factor(position)), alpha = 0.75)+
  # facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
  scale_colour_manual(values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  #scale_shape_manual(values=c(3))+
  guides(fill=FALSE,shape=FALSE, colour=FALSE)+ggtitle("Indicator in Center & Edge")+ylab(NULL)+xlab(NULL)


# plot_bar(ind.rare.rel,x="position",fill="Genus")+
#   facet_wrap(~Genus,scales="free")+theme(legend.position = "none")

#edge only
edgeonly=subset(indval.rare.out.sig,index==2)
goodtaxa <- rownames(edgeonly)
allTaxa <- taxa_names(ps.all.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
#ps <- tax_glom(ps.all.rare.trim, "Genus")
#all.rare.rel <- transform_sample_counts(ps, function(x) x / sum(x))
#all.rare.rel2=merge_samples(all.rare.rel,"position")
#sample_data(all.rare.rel2)$position=c("Center","Edge","Sed")
ps2 <- transform_sample_counts(ps.all.rare.trim, function(x) x / sum(x))
ind.rare.rel <- prune_taxa(allTaxa, ps2)
#plot_bar(ind.rare.rel,x="position",fill="Genus")+
#  facet_wrap(~Family,scales="free")+theme(legend.position = "none")
melt.ind.edgeonly.rare.rel<- psmelt(ind.rare.rel)
melt.ind.edgeonly.rare.rel.sum <- summarySE(melt.ind.edgeonly.rare.rel,measurevar="Abundance",groupvars = c("OTU","Genus","position","Family","Order","Class"))
melt.ind.edgeonly.rare.rel.sum$Family=as.character(melt.ind.edgeonly.rare.rel.sum$Family)

melt.ind.edgeonly.rare.rel.sum[94:96,4]="Solibacteracea"
melt.ind.edgeonly.rare.rel.sum[199:201,4]="NA"


ggplot(melt.ind.edgeonly.rare.rel.sum, aes(x = position, y = Family)) + 
  geom_point(aes(size = Abundance, colour = as.factor(position)), alpha = 0.75)+
  # facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
  scale_colour_manual(values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  #scale_shape_manual(values=c(3))+
  guides(fill=FALSE,shape=FALSE, colour=FALSE)+ggtitle("Indicator at Edge Only")+ylab(NULL)+xlab(NULL)

#center only
centeronly=subset(indval.rare.out.sig,index==1)

#one asv 851
#Taxonomy Table:     [1 taxa by 7 taxonomic ranks]:Kingdom    Phylum           Class                 Order                   Family            
#ASV_851 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Betaproteobacteriales" "Burkholderiaceae"
#Genus                                        Species
#ASV_851 "Burkholderia-Caballeronia-Paraburkholderia" NA  

goodtaxa <- rownames(centeronly)
goodtaxa
allTaxa <- taxa_names(ps.all.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]

#ps <- tax_glom(ps.all.rare.trim, "Genus")
ps2 <- transform_sample_counts(ps.all.rare.trim, function(x) x / sum(x))
ind.rare.rel <- prune_taxa(allTaxa, ps2)
melt.ind.centeronly.rare.rel<- psmelt(ind.rare.rel)


#edge+sed
edgesed=subset(indval.rare.out.sig,index==6)
goodtaxa <- rownames(edgesed)
allTaxa <- taxa_names(ps.all.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
#ps <- tax_glom(ps.all.rare.trim, "Genus")
#all.rare.rel <- transform_sample_counts(ps, function(x) x / sum(x))
#all.rare.rel2=merge_samples(all.rare.rel,"position")
#sample_data(all.rare.rel2)$position=c("Center","Edge","Sed")
ps2 <- transform_sample_counts(ps.all.rare.trim, function(x) x / sum(x))
ind.rare.rel <- prune_taxa(allTaxa, ps2)
plot_bar(ind.rare.rel,x="position",fill="Genus")+
  facet_wrap(~Order,scales="free")+theme(legend.position = "none")

melt.ind.edgeANDsed.rare.rel<- psmelt(ind.rare.rel)
melt.ind.edgeANDsed.rare.rel.sum <- summarySE(melt.ind.edgeANDsed.rare.rel,measurevar="Abundance",groupvars = c("OTU","Genus","position","Family","Order","Class"))

melt.ind.edgeANDsed.rare.rel.sum$Family=as.character(melt.ind.edgeANDsed.rare.rel.sum$Family)

melt.ind.edgeANDsed.rare.rel.sum[262:264,4]="Deferrisoma"

ggplot(melt.ind.edgeANDsed.rare.rel.sum, aes(x = as.factor(position), y = Family)) + 
  geom_point(aes(size = Abundance, colour = as.factor(position)), alpha = 0.75)+
  # facet_grid(indic_res~.,scales="free",space="free")+
  theme_cowplot()+
  scale_colour_manual(values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  #scale_shape_manual(values=c(3))+
  guides(fill=FALSE,shape=FALSE,colour=FALSE)+ggtitle("Indicator at Edge and Sediment")+ylab(NULL)+xlab(NULL)





