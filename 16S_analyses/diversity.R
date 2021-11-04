setwd("~/BOSTON/Davies/Fouha Microbiome/dada2")
load("~/BOSTON/Davies/Fouha Microbiome/dada2/16S_PostTrim&Rarefy.RData")



#### alpha diversity ####
#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, ?? obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of ?? that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/??) and the Gini-Simpson index (1 ??? ??).
#get rid of outliers tho so use new.new.new.meta
library(phyloseq); packageVersion("phyloseq")
ps.all.rare=phyloseq(otu_table(data.matrix(seqtab.all.rare.4mcmc), taxa_are_rows=F), 
                          sample_data(new.new.new.meta), 
                          tax_table(tax_tab_phy)) 
onlyedgeps.rare=subset_samples(ps.all.rare, position=="Edge")
onlycenterps.rare=subset_samples(ps.all.rare, position=="Center")
onlysedps.rare=subset_samples(ps.all.rare, position=="Other")

#edge
onlyedgeps.south.rare=subset_samples(onlyedgeps.rare, side=="South")
onlyedgeps.north.rare=subset_samples(onlyedgeps.rare, side=="North")

#center
onlycenterps.south.rare=subset_samples(onlycenterps.rare, side=="South")
onlycenterps.north.rare=subset_samples(onlycenterps.rare, side=="North")

#site 2 only
only2ps.rare=subset_samples(ps.all.rare, site=="2")
only3ps.rare=subset_samples(ps.all.rare, site=="3")
only4ps.rare=subset_samples(ps.all.rare, site=="4")



#plotting
library(ggplot2)
#all
plot_richness(ps.all.rare, x="position", measures=c("Shannon", "Observed","InvSimpson"), color="position") + theme_bw()+
  theme(text = element_text(size = 20), legend.position = "none", axis.text.x = element_text(angle = 45, vjust =.55 ))+
  scale_color_manual(name="Position",values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+xlab(NULL)+
  annotate(geom="text",x=Inf, y = Inf,vjust=1.2, hjust=1.8, label="p < 0.001***",size=4)+ geom_boxplot()

  

#edge
plot_richness(onlyedgeps.rare, x="site", measures=c("Shannon", "Observed","InvSimpson"), color="site") + theme_bw()+
  theme(text = element_text(size = 20), legend.position = "none", axis.text.x = element_text(angle = 0, vjust =.55 ))+
  scale_color_manual(name="Position",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+xlab(NULL)+
  annotate(geom="text",x=Inf, y = Inf,vjust=1.2, hjust=2.5, label="p = ",size=4)+
  geom_boxplot() 

plot_richness(onlycenterps.rare, x="site", measures=c("Shannon", "Observed","InvSimpson"), color="site") + theme_bw()+
  theme(text = element_text(size = 20), legend.position = "none", axis.text.x = element_text(angle = 0, vjust =.55 ))+
  scale_color_manual(name="Position",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+xlab(NULL)+
  annotate(geom="text",x=Inf, y = Inf,vjust=1.2, hjust=2.5, label="p = ",size=4)+
  geom_boxplot() 

plot_richness(onlysedps.rare, x="site", measures=c("Shannon", "Observed","InvSimpson"), color="site") + theme_bw()+
  theme(text = element_text(size = 20), legend.position = "none", axis.text.x = element_text(angle = 0, vjust =.55 ))+
  scale_color_manual(name="Position",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+xlab(NULL)+
  annotate(geom="text",x=Inf, y = Inf,vjust=1.2, hjust=2.5, label="p = ",size=4)+
  geom_boxplot() 
#edge subsets

plot_richness(onlyedgeps.north.rare, x="site", measures=c("Shannon", "Simpson","Observed","InvSimpson"), color="site") + theme_bw()
plot_richness(onlyedgeps.south.rare, x="site", measures=c("Shannon", "Simpson","Observed","InvSimpson"), color="site") + theme_bw()
#center subs
plot_richness(onlycenterps.north.rare, x="site", measures=c("Shannon", "Simpson","Observed","InvSimpson"), color="site") + theme_bw()
plot_richness(onlycenterps.south.rare, x="site", measures=c("Shannon", "Simpson","Observed","InvSimpson"), color="site") + theme_bw()

#site specific
plot_richness(only2ps.rare, x="position", measures=c("Shannon", "Simpson","Observed","InvSimpson"), color="position") + theme_bw()
plot_richness(only3ps.rare, x="position", measures=c("Shannon", "Simpson","Observed","InvSimpson"), color="position") + theme_bw()
plot_richness(only4ps.rare, x="position", measures=c("Shannon", "Simpson","Observed","InvSimpson"), color="position") + theme_bw()


#Normality and Heteroskedasticity: 
df.all <- data.frame(estimate_richness(ps.all.rare, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
df.edge <- data.frame(estimate_richness(onlyedgeps.rare, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
df.center <- data.frame(estimate_richness(onlycenterps.rare, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
df.sed= data.frame(estimate_richness(onlysedps.rare, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))

df.all$id <- rownames(df.all)
df.all.div <- cbind(df.all,new.new.new.meta) #add sample data

df.edge$id <- rownames(df.edge)
df.edge.div <- cbind(df.edge,subset(new.new.new.meta, position=="Edge")) #add sample data

df.center$id <- rownames(df.center)
df.center.div <- cbind(df.center,subset(new.new.new.meta, position=="Center")) #add sample data

df.sed$id <- rownames(df.sed)
df.sed.div <- cbind(df.sed,subset(new.new.new.meta, position=="Other")) #add sample data



library(car)

#Testing edge vs sed vs center
shapiro.test(df.all.div$Observed) #nope
shapiro.test(df.all.div$InvSimpson) #nope
hist(df.all.div$Shannon)

shapiro.test(sqrt(df.all.div$Shannon)) #no
hist(df.all.div$Observed)
leveneTest(sqrt(df.all.div$Observed)~position,data=df.all.div) #n
#all data not normal or homoskedastic for all three
#obs sigs with welch's
t.test(df.edge.div$Observed,df.center.div$Observed,var.equal=FALSE,conf.level=0.95)
#data:  df.edge.div$Observed and df.center.div$Observed
#t = 7.4856, df = 22.038, p-value = 1.729e-07
##alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  331.4054 585.3737
#sample estimates:
#  mean of x mean of y 
#487.8696   29.4800 
t.test(df.edge.div$Observed,df.sed.div$Observed,var.equal=FALSE,conf.level=0.95)
#Welch Two Sample t-test

#data:  df.edge.div$Observed and df.sed.div$Observed
#t = -10.452, df = 38.989, p-value = 7.225e-13
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1009.2060  -681.9437
#sample estimates:
#  mean of x mean of y 
#487.8696 1333.4444
t.test(df.center.div$Observed,df.sed.div$Observed,var.equal=FALSE,conf.level=0.95)
#Welch Two Sample t-test

#data:  df.center.div$Observed and df.sed.div$Observed
#t = -24.639, df = 17.039, p-value = 9.168e-15
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1415.602 -1192.326
##sample estimates:
#  mean of x mean of y 
#29.480  1333.444
# SHannon sigs with welch's
t.test(df.edge.div$Shannon,df.center.div$Shannon,var.equal=FALSE,conf.level=0.95)
# Welch Two Sample t-test
# 
# data:  df.edge.div$Shannon and df.center.div$Shannon
# t = 11.686, df = 27.505, p-value = 3.535e-12
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   2.224933 3.171626
# sample estimates:
#   mean of x mean of y 
# 4.716358  2.018078
t.test(df.edge.div$Shannon,df.sed.div$Shannon,var.equal=FALSE,conf.level=0.95)
# data:  df.edge.div$Shannon and df.sed.div$Shannon
# t = -6.0963, df = 31.039, p-value = 9.284e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.9642551 -0.9794828
# sample estimates:
#   mean of x mean of y 
# 4.716358  6.188226
t.test(df.center.div$Shannon,df.sed.div$Shannon,var.equal=FALSE,conf.level=0.95)
# Welch Two Sample t-test
# 
# data:  df.center.div$Shannon and df.sed.div$Shannon
# t = -32.404, df = 36.531, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -4.431016 -3.909281
# sample estimates:
#   mean of x mean of y 
# 2.018078  6.188226
#InvSimpson with welch's
t.test(df.edge.div$Observed,df.center.div$Observed,var.equal=FALSE,conf.level=0.95)
# data:  df.edge.div$Observed and df.center.div$Observed
# t = 7.4856, df = 22.038, p-value = 1.729e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   331.4054 585.3737
# sample estimates:
#   mean of x mean of y 
# 487.8696   29.480
t.test(df.edge.div$Observed,df.sed.div$Observed,var.equal=FALSE,conf.level=0.95) 
# data:  df.edge.div$Observed and df.sed.div$Observed
# t = -10.452, df = 38.989, p-value = 7.225e-13
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1009.2060  -681.9437
# sample estimates:
#   mean of x mean of y 
# 487.8696 1333.4444
t.test(df.center.div$Observed,df.sed.div$Observed,var.equal=FALSE,conf.level=0.95)
#data:  df.center.div$Observed and df.sed.div$Observed
#t = -24.639, df = 17.039, p-value = 9.168e-15


#edge normal
shapiro.test(df.edge.div$Observed) #yeah
shapiro.test(df.edge.div$Shannon) #yes
shapiro.test(log(df.edge.div$InvSimpson)) #yeah
#edge hetero
leveneTest(df.edge.div$Observed~site,data=df.edge.div) #homosked
leveneTest(df.edge.div$Shannon~site,data=df.edge.div) #homosked
leveneTest(log(df.edge.div$InvSimpson)~site,data=df.edge.div) #homosked
#aovs
summary(aov(df.edge.div$Observed~site,data=df.edge.div)) 
#            Df  Sum Sq Mean Sq F value Pr(>F)
# site         2  145755   72878   0.833  0.449
# Residuals   20 1750435   87522 
summary(aov(df.edge.div$Shannon~site,data=df.edge.div)) 
# Df Sum Sq Mean Sq F value Pr(>F)
# site         2   2.90   1.452   0.492  0.619
# Residuals   20  59.05   2.953       
summary(aov(log(df.edge.div$InvSimpson)~site,data=df.edge.div)) 
# Df Sum Sq Mean Sq F value Pr(>F)
# site         2   1.73  0.8652   0.723  0.498
# Residuals   20  23.94  1.1968

#center normal
shapiro.test(df.center.div$Observed) #yeah
shapiro.test(log(df.center.div$Shannon)) #yes
shapiro.test(1/(df.center.div$InvSimpson)) #none

#center hetero
leveneTest(df.center.div$Observed~site,data=df.center.div) #homosked
leveneTest(log(df.center.div$Shannon)~site,data=df.center.div) #homosked
leveneTest(df.center.div$InvSimpson~site,data=df.center.div) #homosked
#aovs
summary(aov(df.center.div$Observed~site,data=df.center.div)) 
#Df Sum Sq Mean Sq F value Pr(>F)
#site         2     85   42.54   0.258  0.774
#Residuals   47   7757  165.04 
summary(aov(log(df.center.div$Shannon)~site,data=df.center.div)) 
#Df Sum Sq Mean Sq F value Pr(>F)
#site         2  2.286  1.1430   1.997  0.147
#Residuals   47 26.896  0.5723

#invSimpson not normal so need to do kruskal
kruskal.test(df.center.div$InvSimpson~site,data=df.center.div)
#Kruskal-Wallis rank sum test
#data:  df.center.div$InvSimpson by site
#Kruskal-Wallis chi-squared = 1.7213, df = 2, p-value = 0.4229


#sed normal
shapiro.test(1/log(df.sed.div$Observed)) #none
shapiro.test(df.sed.div$Shannon) #yes
shapiro.test(df.sed.div$InvSimpson) #yes

#sed hetero
leveneTest(df.sed.div$Observed~site,data=df.sed.div) #homosked
leveneTest(df.sed.div$Shannon~site,data=df.sed.div) #homosked
leveneTest(df.sed.div$InvSimpson~site,data=df.sed.div) #homosked

#aovs observed not normal so need to do kruskal
kruskal.test(df.sed.div$Observed~site,data=df.sed.div)
#data:  df.sed.div$Observed by site
#Kruskal-Wallis chi-squared = 0.24561, df = 2, p-value = 0.8844
summary(aov(log(df.sed.div$Shannon)~site,data=df.sed.div)) 
#Df  Sum Sq  Mean Sq F value Pr(>F)
#site         2 0.00551 0.002754   0.481  0.627
#Residuals   15 0.08581 0.005720
#
summary(aov(df.sed.div$InvSimpson~site,data=df.sed.div)) 
#Df Sum Sq Mean Sq F value Pr(>F)
#site         2  17218    8609   0.568  0.578
#Residuals   15 227441   15163      

###BOXPLOTS for ALPHA DIVERSITY


###Beta Diversity 

#### pcoa plots ####
library(stats)
library(MCMC.OTU)
library(vegan)
library(cowplot)
#transform to relative abundance
ps.trim.rel <- transform_sample_counts(ps.all.trim, function(x) x / sum(x))
seq.trim.rel <- data.frame(otu_table(ps.trim.rel))

#Adonis
#All
#rel
sampdata=data.frame(as(sample_data(ps.trim.rel), "matrix"))
adonis(seq.trim.rel ~ position, data=sampdata, permutations = 1000, method="bray")
# Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
# position   2   11.9767  5.9884  89.666 0.67082 0.000999 ***
#   Residuals 88    5.8771  0.0668         0.32918             
# Total     90   17.8538                 1.00000          
#non-rel
sampdata=data.frame(as(sample_data(ps.all.trim), "matrix"))
adonis(seqtab.all.trim[-c(1)] ~ position, data=sampdata, permutations = 1000, method="bray")
#pairwise

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
sampdata=data.frame(as(sample_data(ps.trim.rel), "matrix"))
pairwise.adonis(seq.trim.rel, factors=sampdata$position,  perm = 1000, sim.method="bray")


sampdata=data.frame(as(sample_data(ps.all.trim), "matrix"))
pairwise.adonis(seqtab.all.trim[-c(1)], factors=sampdata$position, perm = 1000, sim.method="bray")
?pairwise.adonis
###### BETA FUN #####
#Comparing position a at site x with position b at site x
#sTOLEN FROM https://rdrr.io/github/jeffkimbrel/jakR/src/R/plotDistances.R

library(dplyr)
library(tidyr)
library(reshape2)

beta.fun.meta=new.new.new.meta
beta.fun.meta$side.site.pos=paste(beta.fun.meta$site,beta.fun.meta$position)
beta.fun.meta$site.pos=paste(beta.fun.meta$site,beta.fun.meta$position)
ROWNAMES=row.names(beta.fun.meta)
ITS2Meta= read.csv("~/BOSTON/Davies/Fouha Microbiome/ITS2_code/ITS2.meta.csv")[,c(7,9)]
beta.fun.meta=merge(ITS2Meta,beta.fun.meta, by="Site",all.y=TRUE)
rownames(beta.fun.meta)=ROWNAMES


ps.all.beta.fun=phyloseq(otu_table(data.matrix(seqtab.all.rare.4mcmc), taxa_are_rows=F), 
                     sample_data(beta.fun.meta), 
                     tax_table(tax_tab_phy)) 
ps.all.beta.fun.rel <- transform_sample_counts(ps.all.beta.fun, function(x) x / sum(x))
#dont load vegan b4 this
test2=distance(ps.all.beta.fun.rel, method="bray", type="sample")
class(test2)
#install.packages("usedist")
library(usedist)

test3=dist_setNames(test2, new.new.new.meta$sample)
wu.m = melt(as.matrix(test3))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(ps.all.beta.fun.rel) %>%
  select(sample, site.pos) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data

colnames(sd) = c("Var1", "Type1")


wu.sd = left_join(wu.m, sd, by= c('Var1' = 'Var1'))

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")


# plot
library(ggplot2)
p = ggplot(wu.sd, aes(x = Type2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
  scale_color_identity()+
  facet_wrap(~ Type1, scales = "free_x") +
 theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))# +
  #ggtitle(paste0("Distance Metric = ", m)) +
  #ylab(m) +
  #xlab(d)
p
# return

#Summary: no obvious trends here, except that edges are closer to sediment than centers are. 
#An edge colony from site x is not more similar to a sediment colony from site x compared to a sediment colony from site y though

# wu.sd$samesite <- ifelse(grepl("3", wu.sd$Type1) & grepl("3",wu.sd$Type2), "yes3",
#                          ifelse(grepl("2", wu.sd$Type1) & grepl("2",wu.sd$Type2), "yes2", ))



###Beta fun within colonies


?merge

test3=dist_setNames(test2, new.new.new.meta$sample)
wu.m = melt(as.matrix(test3))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(ps.all.beta.fun.rel) %>%
  select(sample, site.pos, Colony.x) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data

colnames(sd) = c("Var1", "Type1","Colony")


wu.sd = left_join(wu.m, sd, by= c('Var1' = 'Var1'))
length(unique(sd$Colony))

colnames(sd) = c("Var2", "Type2","Colony")
wu.sd = left_join(wu.sd, sd, by = "Var2")

mteq <- wu.sd[wu.sd$Colony.x==wu.sd$Colony.y, ]
names(mteq)[names(mteq) == 'Var1'] <- 'sample'
merged.mteq=na.omit(merge(mteq,beta.fun.meta,by.x="sample", by.y="sample"))

#adding switching variable
ITS2Switching= read.csv("~/BOSTON/Davies/Fouha Microbiome/ITS2_code/ITS2_switchmeta.csv")[-c(36,6,24,18,40),]
ITS2Meta= read.csv("~/BOSTON/Davies/Fouha Microbiome/ITS2_code/ITS2.meta.csv")
ITS2Switching_2=merge(ITS2Switching, ITS2Meta, by="X")
merged.mteq.ITS2=merge(merged.mteq,ITS2Switching_2, by="Site" )
# paired t-test
Same=subset(merged.mteq.ITS2, Unique=="Switched")
Switched=subset(merged.mteq.ITS2, Unique=="Same")
t.test(Same$value.x,Switched$value.y, var.equal=TRUE,conf.level=0.95)

# Two Sample t-test
# 
# data:  Same$value.x and Switched$value.y
# t = -1.2871, df = 36, p-value = 0.2063
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.10959065  0.02449662
# sample estimates:
#   mean of x mean of y 
# 0.8416435 0.8841905 


#all
ggplot(merged.mteq.ITS2, aes(x=Unique,y=value.x))+geom_boxplot() + theme_bw()+
  theme(text = element_text(size = 20), legend.position = "none")+
 # scale_color_manual(name="Position",values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  xlab(NULL)+ylab(label="Bray-Curtis")+
  annotate(geom="text",x=Inf, y = Inf,vjust=1.5, hjust=1.6, label="p = 0.21",size=4)


# #edge/center only
# 
# what=merged.mteq[!grepl('Other',merged.mteq$Type1),]
# what2=what[!grepl('Other',what$Type2),]
# ggplot(what2, aes(x=factor(Colony.x.x),y=value))+geom_point()
# 
# #center/sed only
# what=merged.mteq[!grepl('Edge',merged.mteq$Type1),]
# what2=what[!grepl('Edge',what$Type2),]
# ggplot(what2, aes(x=site,y=value))+geom_point()
# 
# #edge/sed only
# 
# what=merged.mteq[!grepl('Center',merged.mteq$Type1),]
# what2=what[!grepl('Center',what$Type2),]
# ggplot(what2, aes(x=site,y=value))+geom_point()


###Beta Fun Over###

#edge


#transform to relative abundance
onlyedgepsall.trim.rel <- transform_sample_counts(onlyedgepsall.trim, function(x) x / sum(x))
seq.onlyedgepsall.trim.rel <- data.frame(otu_table(onlyedgepsall.trim.rel))
#rel
sampdata=data.frame(as(sample_data(onlyedgepsall.trim.rel), "matrix"))
adonis(seq.onlyedgepsall.trim.rel ~ site, data=sampdata, permutations = 1000, method="bray")
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# site       2   0.15636 0.078181 0.69843 0.06528 0.7333
# Residuals 20   2.23876 0.111938         0.93472       
# Total     22   2.39512                  1.00000
#non-rel
seqtab.edge.trim=data.frame(otu_table(onlyedgepsall.trim))
sampdata=data.frame(as(sample_data(onlyedgepsall.trim), "matrix"))
adonis(seqtab.edge.trim ~ site, data=sampdata, permutations = 1000, method="bray")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#site       2   0.21204 0.10602  0.7679 0.07131 0.6783
#Residuals 20   2.76126 0.13806         0.92869       
#Total     22   2.97330                 1.00000   

#center
onlycenterpsall.trim.rel <- transform_sample_counts(onlycenterpsall.trim, function(x) x / sum(x))
seq.onlycenterpsall.trim.rel <- data.frame(otu_table(onlycenterpsall.trim.rel))
#rel
sampdata=data.frame(as(sample_data(onlycenterpsall.trim.rel), "matrix"))
adonis(seq.onlycenterpsall.trim.rel ~ site, data=sampdata, permutations = 1000, method="bray")
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#site       2   0.06279 0.031395  1.0718 0.04362 0.3516
#Residuals 47   1.37674 0.029292         0.95638       
#Total     49   1.43953                  1.00000 
#non-rel
seqtab.center.trim=data.frame(otu_table(onlycenterpsall.trim))
sampdata=data.frame(as(sample_data(onlycenterpsall.trim), "matrix"))
adonis(seqtab.center.trim ~ site, data=sampdata, permutations = 1000, method="bray")
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#site       2    0.2443 0.122149  1.4702 0.05888 0.1548
#Residuals 47    3.9049 0.083084         0.94112       
#Total     49    4.1492                  1.00000 


#sed
onlysedpsall.trim.rel <- transform_sample_counts(onlysedpsall.trim, function(x) x / sum(x))
seq.onlysedpsall.trim.rel <- data.frame(otu_table(onlysedpsall.trim.rel))
#rel
sampdata=data.frame(as(sample_data(onlysedpsall.trim.rel), "matrix"))
adonis(seq.onlysedpsall.trim.rel ~ site, data=sampdata, permutations = 1000, method="bray")
#Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
#site       2   0.45116 0.22558  2.1264 0.22089 0.03097 *
#  Residuals 15   1.59127 0.10608         0.77911          
#Total     17   2.04242                 1.00000          
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

sampdata=data.frame(as(sample_data(onlysedpsall.trim), "matrix"))
pairwise.adonis(seq.onlysedpsall.trim.rel[-c(1)], factors=sampdata$site, perm = 1000, sim.method="bray")
#pairs Df  SumsOfSqs   F.Model         R2     p.value p.adjusted sig
#1 3 vs 2  1 0.32598751 3.0708372 0.23493806 0.022977023 0.06893107    
#2 3 vs 4  1 0.08287068 0.6864206 0.06423297 0.637362637 1.00000000    
#3 2 vs 4  1 0.26758123 2.9273881 0.22644854 0.006993007 0.02097902  

#sed all 6
#onlysedpsall.trim.rel@sam_data$side.site.pos=paste(onlysedpsall.trim.rel@sam_data$position,onlysedpsall.trim.rel@sam_data$side.site)
ord <- ordinate(onlysedpsall.trim.rel, "PCoA", "bray")
plot_ordination(onlysedpsall.trim.rel, ord,color="side.site", shape="side.site")+
  stat_ellipse()

gg.pcoa.site.rare <- 
  plot_ordination(onlysedpsall.trim.rel, ord,color="side.site")+
  stat_ellipse(size=1)+
  theme_cowplot()+theme(text = element_text(size = 20))+
  #scale_color_manual(name="Position",values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  #scale_shape_manual(name="Position",values=c(21,16,16))+
  xlab("Axis 1 (68.1%)")+
  ylab("Axis 2 (7%)")+
  annotate(geom="text", x=-0.25, y=0.65, label="padj < 0.001***",size=6) #+
#ggtitle("Rarefied")
gg.pcoa.site.rare


#not rel
seqtab.sed.trim=data.frame(otu_table(onlysedpsall.trim))
sampdata=data.frame(as(sample_data(onlysedpsall.trim), "matrix"))
adonis(seqtab.sed.trim ~ site, data=sampdata, permutations = 1000, method="bray")
# Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
# site       2   0.46786 0.23393  2.0988 0.21865 0.02697 *
#   Residuals 15   1.67185 0.11146         0.78135          
# Total     17   2.13971                 1.00000          
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ord <- ordinate(ps.trim.rel, "PCoA", "bray")
plot_ordination(ps.trim.rel, ord,color="position", shape="position")+
  stat_ellipse()

gg.pcoa.site.rare <- 
  plot_ordination(ps.trim.rel, ord,color="position")+
  stat_ellipse(size=1)+
  theme_cowplot()+theme(text = element_text(size = 20))+
  scale_color_manual(name="Position",values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  #scale_shape_manual(name="Position",values=c(21,16,16))+
  xlab("Axis 1 (68.1%)")+
  ylab("Axis 2 (7%)")+
  annotate(geom="text", x=-0.25, y=0.65, label="padj < 0.001***",size=6) #+
  #ggtitle("Rarefied")
gg.pcoa.site.rare

#All accross 6 sites
ps.trim.rel@sam_data$side.site.pos=paste(ps.trim.rel@sam_data$position,ps.trim.rel@sam_data$side.site)

ord <- ordinate(ps.trim.rel, "PCoA", "bray")
plot_ordination(ps.trim.rel, ord,color="side.site.pos", shape="position")+
  stat_ellipse()

gg.pcoa.site.rare <- 
  plot_ordination(ps.trim.rel, ord,color="side.site.pos",shape="position")+
  stat_ellipse(size=1)+
  theme_cowplot()+theme(text = element_text(size = 20))+
  #scale_color_manual(name="Position",values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  #scale_shape_manual(name="Position",values=c(21,16,16))+
  xlab("Axis 1 (68.1%)")+
  ylab("Axis 2 (7%)")+
  annotate(geom="text", x=-0.25, y=0.65, label="padj < 0.001***",size=6) #+
#ggtitle("Rarefied")
gg.pcoa.site.rare

#NO SED
ps.trim.rel.nosed=subset_samples(ps.trim.rel, position!="Other")
ord <- ordinate(ps.trim.rel.nosed, "PCoA", "bray")


gg.pcoa.site.rare <- 
  plot_ordination(ps.trim.rel.nosed, ord,color="side.site.pos",shape="position")+
  stat_ellipse(size=1)+
  theme_cowplot()+theme(text = element_text(size = 20))+
  #scale_color_manual(name="Position",values=c("#ff9955ff","#aa4400ff","#4d4d4dff"))+
  #scale_shape_manual(name="Position",values=c(21,16,16))+
  xlab("Axis 1 (68.1%)")+
  ylab("Axis 2 (7%)")+
  annotate(geom="text", x=-0.25, y=0.65, label="padj < 0.001***",size=6) #+
#ggtitle("Rarefied")
gg.pcoa.site.rare



ord <- ordinate(onlyedgepsall.trim.rel, "PCoA", "bray")
plot_ordination(onlyedgepsall.trim.rel, ord,color="site", shape="site")+
  stat_ellipse()

gg.pcoa.site.rare <- plot_ordination(onlyedgepsall.trim.rel, ord,color="site")+
  stat_ellipse(size=1)+
  theme_cowplot()+theme(text = element_text(size = 20))+
  scale_color_manual(name="Site",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+
 # scale_shape_manual(name="Position",values=c(21,16,16))+
  xlab("Axis 1 (45.1%)")+
  ylab("Axis 2 (12%)")+
  annotate(geom="text", x=-0.48, y=0.65, label="p = 0.733",size=6) +
ggtitle("Edge")
gg.pcoa.site.rare


ord <- ordinate(onlycenterpsall.trim.rel, "PCoA", "bray")
plot_ordination(onlycenterpsall.trim.rel, ord,color="site", shape="site")+
  stat_ellipse()

gg.pcoa.site.rare <- plot_ordination(onlycenterpsall.trim.rel, ord,color="site")+
  stat_ellipse(size=1)+
  theme_cowplot()+theme(text = element_text(size = 20))+
  scale_color_manual(name="Site",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+
  # scale_shape_manual(name="Position",values=c(21,16,16))+
  xlab("Axis 1 (58.2%)")+
  ylab("Axis 2 (15.4%)")+
  annotate(geom="text", x=-.05, y=0.25, label="p = 0.351",size=6) +
  ggtitle("Center")
gg.pcoa.site.rare

ord <- ordinate(onlysedpsall.trim.rel, "PCoA", "bray")
plot_ordination(onlysedpsall.trim.rel, ord,color="site", shape="site")+
  stat_ellipse()

gg.pcoa.site.rare <- plot_ordination(onlysedpsall.trim.rel, ord,color="site")+
  stat_ellipse(size=1)+
  theme_cowplot()+theme(text = element_text(size = 20))+
  scale_color_manual(name="Site",values=c("#4d4d4dff","#87cddeff","#0055d4ff"))+
  # scale_shape_manual(name="Position",values=c(21,16,16))+
  xlab("Axis 1 (41.6%)")+
  ylab("Axis 2 (14.2%)")+
  annotate(geom="text", x=-0.4, y=0.65, label="2vs3 padj = 0.069",size=6) +
  annotate(geom="text", x=-0.4, y=0.58, label="3vs4 padj = 1.0",size=6) +
  annotate(geom="text", x=-0.4, y=0.51, label="2vs4 padj < 0.02*",size=6) +
  
  ggtitle("Sediment")
gg.pcoa.site.rare

