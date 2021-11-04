
load("~/BOSTON/Davies/Fouha Microbiome/ITS2_code/PostTrimRarefy_ITS.RData")
setwd("~/BOSTON/Davies/Fouha Microbiome/ITS2_code/")

#write.csv(new.meta, file="ITS2.meta.csv")
#
library(reshape2)
library(RColorBrewer)
library(ggplot2)
#courtney Klepac's symportal data
SymPortalDataK=read.table(file="184_20211015_03_DBV_20211018T045550.profiles.relative.abund_and_meta.txt", skip=6, sep="\t",header=T)[-c(28,29,32,59,60),-1]

library(stringr)
SymPortalDataK$label = str_extract(SymPortalDataK$X, "s|t")
SymPortalDataK$var = str_extract(SymPortalDataK$X, "LV|MV|HV")
SymPortalDataK$colony=str_sub(SymPortalDataK$X, 1,3)
row.names(SymPortalDataK)=SymPortalDataK$X
new.melt=melt(SymPortalDataK,id=c('X','label','var','colony'))
SymPortalDataK_1=subset(new.melt,value != 0 )


test=subset(SymPortalDataK_1, duplicated(SymPortalDataK_1$colony) | duplicated(SymPortalDataK_1$colony, fromLast = TRUE))

dd <- union(SymPortalDataK_1$variable,test$variable)
library(RColorBrewer)
n <- 59
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector
length(col_vector)<-length(dd)
colours <- setNames(as.character(col_vector), 
                    dd)

ok.all=as.character(unique(SymPortalDataK_1$variable))
ok.test=as.character(unique(test$variable))


col.all=colours[ok.all]
col.test=colours[ok.test]

SymPortalDataK_1$var_f = factor(SymPortalDataK_1$var, levels=c('LV','MV','HV'))

ggplot(SymPortalDataK_1, aes(x=as.factor(X), y=as.numeric(value), fill=variable))+ geom_bar(aes(), stat="identity") +
  scale_fill_manual(values=col.all,na.value="black")+
  theme_bw()+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  guides(fill=guide_legend(nrow=5))+
  scale_y_continuous(limits = c(0, 1), breaks=c(.25,.5,.75,1))+xlab(NULL)+ylab(NULL)+facet_grid(~var_f, 
                                                                                                scales = "free_x", # Let the x axis vary across facets.
                                                                                                space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                                                                                switch = "x")      # Move the facet labels to the bottom.
p1
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(p1)
p3=ggplot(test[test$label %in% c('t'),], aes(x=as.factor(colony), y=as.numeric(value), fill=variable))+ geom_bar(aes(), stat="identity") +
  scale_fill_manual(values=col.test,na.value="black")+
  theme_bw()+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  guides(fill=guide_legend(nrow=5))+
  scale_y_continuous(limits = c(0, 1), breaks=c(.25,.5,.75,1))+xlab(NULL)+ylab(NULL)+facet_grid(~label*var, 
                                                                                                scales = "free_x", # Let the x axis vary across facets.
                                                                                                space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                                                                                switch = "x")      # Move the facet labels to the bottom.
p3
p2=ggplot(test[test$label %in% c('s'),], aes(x=as.factor(colony), y=as.numeric(value), fill=variable))+ geom_bar(aes(), stat="identity") +
  scale_fill_manual(values=col.test,na.value="black")+
  theme_bw()+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  guides(fill=guide_legend(nrow=5))+
  scale_y_continuous(limits = c(0, 1),breaks=c(.25,.5,.75,1))+xlab(NULL)+ylab(NULL)+facet_grid(~label*var, 
                                                                                               scales = "free_x", # Let the x axis vary across facets.
                                                                                               space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                                                                               switch = "x")      # Move the facet labels to the bottom.
library(gridExtra)

grid.arrange(p3,p2,nrow=2)

#chi-square, top vs side
chisq.test(SymPortalDataK_1$variable, SymPortalDataK_1$label, correct=FALSE)
#Pearson's Chi-squared test

#data:  SymPortalDataK_1$variable and SymPortalDataK_1$label
#X-squared = 5.5824, df = 6, p-value = 0.4716

#top
SymPortalDataK_1_t=subset(SymPortalDataK_1, label=="t")
chisq.test(SymPortalDataK_1_t$variable, SymPortalDataK_1_t$var, correct=FALSE)
#Pearson's Chi-squared test

#data:  SymPortalDataK_1_t$variable and SymPortalDataK_1_t$var
#X-squared = 15.867, df = 10, p-value = 0.1035

#side
SymPortalDataK_1_s=subset(SymPortalDataK_1, label=="s")

chisq.test(SymPortalDataK_1_s$variable, SymPortalDataK_1_s$var, correct=FALSE)
#Pearson's Chi-squared test

#data:  SymPortalDataK_1_s$variable and SymPortalDataK_1_s$var
#X-squared = 9.2857, df = 8, p-value = 0.3188

########################my symportal data
SymPortalData=read.table(file="177_20210916_all_03_DBV_20210916T225814.profiles.relative.abund_and_meta.txt", skip=6, sep="\t",header=T)[-c(72,73),-1]
row.names(SymPortalData)=SymPortalData$X
new1<-transform(merge(SymPortalData,new.meta,by=0), row.names=Row.names, Row.names=NULL)[,-c(15,17,18,19,20,21)]

new.melt=melt(new1,id=c('X','site','position','side','Colony','side.site','sample'))

#Abundance plots

n <- 59
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector

new.melt <- new.melt[order(new.melt$site),]
new.melt <- new.melt[order(new.melt$position),]

#top10.center$site<-factor(top10.center$site,levels=unique(top10.center$site))
new.melt$sample<-factor(new.melt$sample,levels=unique(new.melt$sample))

#colonies that have both edge/center samples 
new.melt1=subset(new.melt,value != 0 )
test=subset(new.melt1, duplicated(new.melt1$Colony) | duplicated(new.melt1$Colony, fromLast = TRUE))

dd <- union(new.melt$variable,test$variable)
library(RColorBrewer)
n <- 59
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector
length(col_vector)<-length(dd)
colours <- setNames(as.character(col_vector), 
                    dd)

ok.all=as.character(unique(new.melt$variable))
ok.test=as.character(unique(test$variable))


col.all=colours[ok.all]
col.test=colours[ok.test]



p <- ggplot(data=subset(new.melt,value != 0), aes(x=as.factor(Colony), y=as.numeric(value), fill=variable))
p + geom_bar(aes(), stat="identity") +
  scale_fill_manual(values=col.all,na.value="black")+
  theme_bw()+
  theme(legend.position="right",legend.text = element_text( size=10),legend.key.size = unit(.4, 'cm'),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  scale_y_continuous(limits = c(0, 1), breaks=c(.25,.5,.75,1))+xlab(NULL)+ylab(NULL)+facet_grid(~site*position, 
                                                                                                scales = "free_x", # Let the x axis vary across facets.
                                                                                                space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                                                                                switch = "x")      # Move the facet labels to the bottom.

#subbing by only colonies that have both edge and center
new.melt1=subset(new.melt,value != 0 )
test=subset(new.melt1, duplicated(new.melt1$Colony) | duplicated(new.melt1$Colony, fromLast = TRUE))
 hmm=test %>%
  group_by(Colony, variable) %>%
  mutate(Unique = if_else(duplicated(variable) == TRUE|duplicated(variable, fromLast=TRUE), "Same", "Switched"))
  #mutate(code=ifelse(row_number()%%2==0,'A','a'))
write.csv(hmm,"ITS2_switchmeta.csv")
 
p1= ggplot(test[test$position %in% c('Center'),], aes(x=as.factor(Colony), y=as.numeric(value), fill=variable))+ geom_bar(aes(), stat="identity") +
  scale_fill_manual(values=col.test,na.value="black")+
  theme_bw()+
  theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  guides(fill=guide_legend(nrow=5))+
  scale_y_continuous(limits = c(0, 1), breaks=c(.25,.5,.75,1))+xlab(NULL)+ylab(NULL)+facet_grid(~position*site, 
                                                                        scales = "free_x", # Let the x axis vary across facets.
                                                                        space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                                                        switch = "x")      # Move the facet labels to the bottom.

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(p1)
p3=ggplot(test[test$position %in% c('Center'),], aes(x=as.factor(Colony), y=as.numeric(value), fill=variable))+ geom_bar(aes(), stat="identity") +
  scale_fill_manual(values=col.test,na.value="black")+
  theme_bw()+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  guides(fill=guide_legend(nrow=5))+
  scale_y_continuous(limits = c(0, 1), breaks=c(.25,.5,.75,1))+xlab(NULL)+ylab(NULL)+facet_grid(~position*site, 
                                                                                                scales = "free_x", # Let the x axis vary across facets.
                                                                                                space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                                                                                switch = "x")      # Move the facet labels to the bottom.

p2= ggplot(test[test$position %in% c('Edge'),], aes(x=as.factor(Colony), y=as.numeric(value), fill=variable))+ geom_bar(aes(), stat="identity") +
  scale_fill_manual(values=col.test,na.value="black")+
  theme_bw()+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme_cowplot()+
  guides(fill=guide_legend(nrow=5))+
  scale_y_continuous(limits = c(0, 1),breaks=c(.25,.5,.75,1))+xlab(NULL)+ylab(NULL)+facet_grid(~position*site, 
                                                                        scales = "free_x", # Let the x axis vary across facets.
                                                                        space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                                                                        switch = "x")      # Move the facet labels to the bottom.


grid.arrange(p3,p2,nrow=2)

?grid.arrange

##Chi-square test##
#Should do the same thing but within center/edge accross sites



input=test[,c(3,8)]
new.input=droplevels(input)
new.input
new.input=ctable <- table(new.input)
new.input
CTable <- data.frame(matrix(ctable, nrow = dim(ctable)[1], dimnames = dimnames(ctable)))
droplevels(CTable)
chisq.test(input$variable, input$position, correct=FALSE)

# Pearson's Chi-squared test
# 
# data:  input$variable and input$position
# X-squared = 7.7298, df = 8, p-value = 0.4603

## Subset edge all, look accross site
edge.new.melt1=subset(new.melt1, position=="Edge")

input=edge.new.melt1[,c(2,8)]

chisq.test(input$variable, input$site, correct=FALSE)

#Pearson's Chi-squared test

#data:  input$variable and input$site
#X-squared = 15.754, df = 16, p-value = 0.4703
center.new.melt1=subset(new.melt1, position=="Center")

input=center.new.melt1[,c(2,8)]

chisq.test(input$variable, input$site, correct=FALSE)
#Pearson's Chi-squared test

#data:  input$variable and input$site
#X-squared = 26.682, df = 22, p-value = 0.2235

