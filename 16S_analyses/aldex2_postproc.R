#####Below is what I run after getting aldex results ######
load("aldexoutput_all.RData")
#sig different between center and edge
Sig.PW=subset(C.glm_PW,C.glm_PW[,14]<0.05)
# #A ton are sig, filtering out the pathways with actual names and including only BH <0.05 edge/center comp
# Sig.PW.clean=subset(C.glm_PW[!grepl("^PWY|^P1|^P2|^P4", rownames(C.glm_PW)),], C.glm_PW[!grepl("^PWY|^P1|^P2|^P4", rownames(C.glm_PW)),][,14]<.05)
# #Those that are sig different between sed and edge
Sig.PW.SE=subset(E.glm_PW, E.glm_PW[,15]<0.05)
# #make list of those that are sig different between edge/center AND edge/sed 

Sig.PW.SE.names=rownames(Sig.PW.SE)
Sig.PW.EC.And_SE=Sig.PW[(rownames(Sig.PW) %in% Sig.PW.SE.names),]

#the coefficient for center/edge comp (i.e. model.testEdgeEstimate) is 
#greater than the coefficient for center/sed comp 
#(i.e. model.testOtherEstimate) the coefficient means Edge is X times greater 
#than the OVERALL mean and thus if it is larger than the Other coefficient 
#Edge is greater than Sed.

#Of those sig different btwn edge/center AND edge/sed the following are 
#greatest in Edge then Sed then Center
Sig.PW.Edge_Sed_Center=subset(Sig.PW.EC.And_SE, Sig.PW.EC.And_SE$`model.testEdge Estimate`>Sig.PW.EC.And_SE$`model.testOther Estimate`& Sig.PW.EC.And_SE$`model.testEdge Estimate`>0)
#greatest in Sed then Edge then center 
Sig.PW.Sed_Edge_Center=subset(Sig.PW.EC.And_SE, Sig.PW.EC.And_SE$`model.testEdge Estimate`<Sig.PW.EC.And_SE$`model.testOther Estimate`& Sig.PW.EC.And_SE$`model.testEdge Estimate`>0)
#greatest in Center then Edge then Sed
Sig.PW.Center_Edge_Sed=subset(Sig.PW.EC.And_SE, Sig.PW.EC.And_SE$`model.testEdge Estimate`>Sig.PW.EC.And_SE$`model.testOther Estimate`& Sig.PW.EC.And_SE$`model.testEdge Estimate`<0)
#greatest in Center then Sed then Edge
Sig.PW.Center_Sed_Edge=subset(Sig.PW.EC.And_SE, Sig.PW.EC.And_SE$`model.testEdge Estimate`<Sig.PW.EC.And_SE$`model.testOther Estimate`& Sig.PW.EC.And_SE$`model.testEdge Estimate`<0)

pathnames=getFeatureNames(clr.PW.C)
reads=getReads(clr.PW.C)
sampnames=getSampleIDs(clr.PW.C)
props= lapply( reads, function(x){ x/sum(x, na.rm=TRUE)} )
df <- data.frame(matrix(unlist(props), nrow=length(props), byrow=TRUE))
row.names(df)=sampnames
colnames(df)=pathnames

#Plotting, eventually will need to change row annotation to MYC pathway name         
Group=Sig.PW.Edge_Sed_Center
goodnames=rownames(Group)
df.goodsig <- df[, goodnames]
df.goodsig.2=cbind(sample_data(ps.all.rare.trim),df.goodsig)
#write.csv(df.goodsig.2, file="Sig.PW.EC.SE.clean.csv")
# df.goodsig.2=read.csv(file="Sig.PW.EC.SE.clean.csv",row.names = 1)
# 
# #create z scores for each column
df.goodsig.2[,-c(1:13)] <- scale(df.goodsig.2[,-c(1:13)])
df.goodsig.2= df.goodsig.2[order(df.goodsig.2$position),]
df.goodsig.3 = t(data.matrix(df.goodsig.2[,-c(1:13)]))
library(pheatmap)
pheatmap(df.goodsig.3, show_rownames=T, labels_row= row.names(df.goodsig.3),
         show_colnams = T, angle_col= "90",
         cellheight = 4, cellwidth = 4, fontsize_row=2, cluster_cols=T, scale='row', fontsize_col=2,
         cluster_rows = T, color=colour, labels_col = df.goodsig.2$position,
         treeheight_col = 10, treeheight_row = 10, legend=F,
         width=15, height=25)

Group=Sig.PW.Sed_Edge_Center
goodnames=rownames(Group)
df.goodsig <- df[, goodnames]
df.goodsig.2=cbind(sample_data(ps.all.rare.trim),df.goodsig)
#write.csv(df.goodsig.2, file="Sig.PW.EC.SE.clean.csv")
# df.goodsig.2=read.csv(file="Sig.PW.EC.SE.clean.csv",row.names = 1)
# 
# #create z scores for each column
df.goodsig.2[,-c(1:13)] <- scale(df.goodsig.2[,-c(1:13)])
df.goodsig.2= df.goodsig.2[order(df.goodsig.2$position),]
df.goodsig.3 = t(data.matrix(df.goodsig.2[,-c(1:13)]))
library(pheatmap)
pheatmap(df.goodsig.3, show_rownames=T, labels_row= row.names(df.goodsig.3),
         show_colnams = T, angle_col= "90",
         cellheight = 3, cellwidth = 3, fontsize_row=3, cluster_cols=T, scale='row', fontsize_col=3,
         cluster_rows = T, color=colour, labels_col = df.goodsig.2$position,
         treeheight_col = 10, treeheight_row = 10, legend=F,
         width=20, height=25)

Group=Sig.PW.Center_Edge_Sed
goodnames=rownames(Group)
df.goodsig <- df[, goodnames]
df.goodsig.2=cbind(sample_data(ps.all.rare.trim),df.goodsig)
#write.csv(df.goodsig.2, file="Sig.PW.EC.SE.clean.csv")
# df.goodsig.2=read.csv(file="Sig.PW.EC.SE.clean.csv",row.names = 1)
# 
# #create z scores for each column
df.goodsig.2[,-c(1:13)] <- scale(df.goodsig.2[,-c(1:13)])
df.goodsig.2= df.goodsig.2[order(df.goodsig.2$position),]
df.goodsig.3 = t(data.matrix(df.goodsig.2[,-c(1:13)]))
library(pheatmap)
colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(df.goodsig.3, show_rownames=T, labels_row= row.names(df.goodsig.3),
         show_colnams = T, angle_col= "90",
         cellheight = 3, cellwidth = 3, fontsize_row=3, cluster_cols=T, scale='row', fontsize_col=3,
         cluster_rows = T, color=colour, labels_col = df.goodsig.2$position,
         treeheight_col = 10, treeheight_row = 10, legend=F,
         width=20, height=25)


# heatmap_plot= function(Group){
#         Group=Sig.PW.Edge_Sed_Center
#         pathnames=getFeatureNames(clr.PW.C)
#         reads=getReads(clr.PW.C)
#         sampnames=getSampleIDs(clr.PW.C)
#         props= lapply( reads, function(x){ x/sum(x, na.rm=TRUE)} )
#         df <- data.frame(matrix(unlist(props), nrow=length(props), byrow=TRUE))
#         row.names(df)=sampnames
#         colnames(df)=pathnames
#         goodnames=rownames(Group)
#         df.goodsig <- df[, goodnames]
#         df.goodsig.2=cbind(sample_data(ps.all.rare.trim),df.goodsig)
#         #write.csv(df.goodsig.2, file="Sig.PW.EC.SE.clean.csv")
#         # df.goodsig.2=read.csv(file="Sig.PW.EC.SE.clean.csv",row.names = 1)
#         # 
#         # #create z scores for each column
#         df.goodsig.2[,-c(1:13)] <- scale(df.goodsig.2[,-c(1:13)])
#         df.goodsig.2= df.goodsig.2[order(df.goodsig.2$position),]
#         df.goodsig.3 = t(data.matrix(df.goodsig.2[,-c(1:13)]))
#         
#         # #heatmap
#         library(gplots)
#         colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#         colnames(df.goodsig.3)
#         heatmap.2(,col = colour, Rowv = TRUE, Colv = TRUE, scale = "row", 
#                   dendrogram = "both",
#                   trace = "none", labRow = row.names(df.goodsig.3),cexRow=.3, labCol = df.goodsig.2$position,cexCol =.3,
#                   margin = c(3,6))
# }
# heatmap_plot(Sig.PW.Edge_Sed_Center)
# heatmap_plot(Sig.PW.Sed_Edge_Center)
# heatmap_plot(Sig.PW.Center_Edge_Sed)
# heatmap_plot(Sig.PW.Center_Sed_Edge)


# 
# 
# 
# #make list of those that are sig different between edge/center AND NOT edge/sed
Sig.PW.EC.Not_SE=Sig.PW[!(rownames(Sig.PW) %in% Sig.PW.SE.names),]
Group=Sig.PW.EC.Not_SE
goodnames=rownames(Group)
df.goodsig <- df[, goodnames]
df.goodsig.2=cbind(sample_data(ps.all.rare.trim),df.goodsig)
#write.csv(df.goodsig.2, file="Sig.PW.EC.SE.clean.csv")
# df.goodsig.2=read.csv(file="Sig.PW.EC.SE.clean.csv",row.names = 1)
# 
# #create z scores for each column
df.goodsig.2[,-c(1:13)] <- scale(df.goodsig.2[,-c(1:13)])
df.goodsig.2= df.goodsig.2[order(df.goodsig.2$position),]
df.goodsig.3 = t(data.matrix(df.goodsig.2[,-c(1:13)]))
library(pheatmap)
pheatmap(df.goodsig.3, show_rownames=T, labels_row= row.names(df.goodsig.3),
         show_colnams = T, angle_col= "90",
         cellheight = 3, cellwidth = 3, fontsize_row=3, cluster_cols=T, scale='row', fontsize_col=3,
         cluster_rows = T, color=colour, labels_col = df.goodsig.2$position,
         treeheight_col = 10, treeheight_row = 10, legend=F,
         width=20, height=25)


# 
# #plot
# pathnames=getFeatureNames(clr.PW.C)
# reads=getReads(clr.PW.C)
# sampnames=getSampleIDs(clr.PW.C)
# props= lapply( reads, function(x){ x/sum(x, na.rm=TRUE)} )
# df <- data.frame(matrix(unlist(props), nrow=length(props), byrow=TRUE))
# row.names(df)=sampnames
# colnames(df)=pathnames
# 
# goodnames=rownames(Sig.PW.EC.Not_SE)
# df.goodsig <- df[, goodnames]
# df.goodsig.2=cbind(sample_data(ps.all.rare.trim),df.goodsig)
# write.csv(df.goodsig.2, file="Sig.PW.EC.Not_SE.csv")
# df.goodsig.2=read.csv(file="Sig.PW.EC.Not_SE.csv",row.names = 1)
# 
# #create z scores for each column
# df.goodsig.2[,13:16] <- scale(df.goodsig.2[,13:16])
# df.goodsig.2= df.goodsig.2[order(df.goodsig.2$position),]
# df.goodsig.3 = t(data.matrix(df.goodsig.2[,13:16]))
# 
# #heatmap
# library(gplots)
# colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
# colnames(df.goodsig.3)
# heatmap.2(df.goodsig.3, col = colour, Rowv = TRUE, Colv = TRUE, scale = "row", 
#           dendrogram = "both",
#           trace = "none", labRow = row.names(df.goodsig.3),cexRow=.3, labCol = df.goodsig.2$position,cexCol =.3,
#           margin = c(3,6))
# 
# 
# 
# #plotting everything, first extract counts for relevant pathways and transform to proportions
# pathnames=getFeatureNames(clr.PW.C)
# reads=getReads(clr.PW.C)
# sampnames=getSampleIDs(clr.PW.C)
# props= lapply( reads, function(x){ x/sum(x, na.rm=TRUE)} )
# df <- data.frame(matrix(unlist(props), nrow=length(props), byrow=TRUE))
# row.names(df)=sampnames
# colnames(df)=pathnames
# 
# goodnames=rownames(Sig.PW.clean)
# df.goodsig <- df[, goodnames]
# df.goodsig.2=cbind(sample_data(ps.all.rare.trim),df.goodsig)
# write.csv(df.goodsig.2, file="picrust.sig.prop.csv")
# df.goodsig.2=read.csv(file="picrust.sig.prop.csv",row.names = 1)
# 
# #create z scores for each column
# df.goodsig.2[,13:100] <- scale(df.goodsig.2[,13:100])
# df.goodsig.2= df.goodsig.2[order(df.goodsig.2$position),]
# df.goodsig.3 = t(data.matrix(df.goodsig.2[,13:100]))
# 
# #heatmap
# library(gplots)
# colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
# colnames(df.goodsig.3)
# heatmap.2(df.goodsig.3, col = colour, Rowv = TRUE, Colv = TRUE, scale = "row", 
#           dendrogram = "both",
#           trace = "none", labRow = row.names(df.goodsig.3),cexRow=.3, labCol = df.goodsig.2$position,cexCol =.3,
#           margin = c(3,6))
####Plottting only the PWs that are sig diff in center/edge but not edge/sed


#re-run clr object with center vs edge comparison and then filter out those that were significant from here?
#NAH
# ps_EdCen=subset_samples(ps.all.rare.trim, position=="Center"|position=="Edge")
# p2PW_EdCen= p2PW[,sample_names(ps_EdCen)]
# 
# set.seed(12345)
# system.time({
#   aldex2_PW_EdCen = aldex(p2PW_EdCen, sample_data(ps_EdCen)$position, mc.samples = 500, test = "t", 
#                      effect = TRUE, denom = "iqlr", verbose = TRUE)
# })

##Take pvals from glm or t.test, I dont think it really matters. Use effect size for heatmap. Filter
#out pathways so it is only the coolest pathways. 
#If there are too many, could do heatmap with subset and table with rest. 
#If I do a hm should do for each individual sample though, alternatively could just box plot
#those clr numbers and group according to edge/center. 
#shoudl compare edge, sed and center and use the clr numbers. 
