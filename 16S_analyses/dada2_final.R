#https://benjjneb.github.io/dada2/tutorial.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10")

library(dada2); packageVersion("dada2")
library(ShortRead)
library(Biostrings)
library(phyloseq); packageVersion("phyloseq")
library(tibble)

setwd("~/BOSTON/Davies/Fouha Microbiome/dada2")
#this was all done in VM, because windows doesn't like cutadapt
path <- "/media/jamesf/TOSHIBA_EXT/R_stuff/Fouha/16S"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_NoIll_16S.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_NoIll_16S.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
library(stringr)
sample.names=str_remove(sample.names, "FouhaFinalJF_")
sample.names

FWD <- "GTGYCAGCMGCCGCGGTA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)



fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
#filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)



primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/home/jamesf/.local/bin/cutadapt"  
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
####Output
# === Second read: Adapter 4 ===
#   
#   Sequence: TACCGCGGCKGCTGRCAC; Type: regular 3'; Length: 18; Trimmed: 94 times.
# 
# No. of allowed errors:
# 0-9 bp: 0; 10-18 bp: 1
# 
# Bases preceding removed adapters:
#   A: 4.3%
#   C: 1.1%
#   G: 5.3%
#   T: 89.4%
#   none/other: 0.0%
# WARNING:
#     The adapter is preceded by "T" extremely often.
#     The provided adapter sequence could be incomplete at its 3' end.
# 
# Overview of removed sequences
# length	count	expect	max.err	error counts
# 3	18	865.0	0	18
# 4	5	216.2	0	5
# 17	45	0.0	1	31 14
# 18	2	0.0	1	2
# 19	1	0.0	1	1
# 22	2	0.0	1	2
# 28	1	0.0	1	1
# 30	1	0.0	1	1
# 36	1	0.0	1	1
# 37	13	0.0	1	13
# 40	1	0.0	1	1
# 43	1	0.0	1	1
# 48	1	0.0	1	1
# 73	1	0.0	1	0 1
# 215	1	0.0	1	0 1

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

cutFs <- sort(list.files(path.cut, pattern = "R1_NoIll_16S.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_NoIll_16S.fastq", full.names = TRUE))
cutFs

#checkings forwards
plotQualityProfile(cutFs[2])

#checking revs
plotQualityProfile(cutRs[2])
#for now even though the proportion of reads drops off at ~190, I think the quality of reads drops at ~230 so a 225 cutoff was used 
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
#Initial settings
#out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(225,225),
#                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                     compress=TRUE, multithread=TRUE)
#led to a lot being filtered out:
# > head(track)
# input filtered denoisedF denoisedR merged nonchim
# A10_S10 56669    27466     27341     27355  26992   25542
# A11_S11 44084    22046     20922     20822  20012   18598
# A12_S12 56557    33454     30163     30092  27043   24948
# A1_S1   13646    11583     11483     11516  11260   11260
# A2_S2   48291    18587     17086     17237  15894   15430
# A3_S3   56168    47053     41360     41127  33957   31870

#2nd try:
# out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(225,225),
#                                           maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
#                                           compress=TRUE, multithread=TRUE)
#The parameter you are probably looking for is truncLen. This will remove all reads less than truncLen, and also truncate the reads at truncLen.

#                                       reads.in reads.out
#FouhaFinalJF_A10_S10_R1_NoIll_16S.fastq    56669     28041
#FouhaFinalJF_A11_S11_R1_NoIll_16S.fastq    44084     22506
#FouhaFinalJF_A12_S12_R1_NoIll_16S.fastq    56557     34290
#FouhaFinalJF_A1_S1_R1_NoIll_16S.fastq      13646     12015
#FouhaFinalJF_A2_S2_R1_NoIll_16S.fastq      48291     19141
#FouhaFinalJF_A3_S3_R1_NoIll_16S.fastq      56168     48658
#So maxEE doesn't do much difference, truncLen is the main culprit here and I think its because its removing all reads less than 225, 50%
#of them are ~190 

#per Nicola's suggestion, if I truncLen at 190 I should still maintain that last ~30-40 basepairs that are real on half of them via the reverse reads. 
#3rd try
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(185,185),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)
#reads.in reads.out
#FouhaFinalJF_A10_S10_R1_NoIll_16S.fastq    56669     55786
#FouhaFinalJF_A11_S11_R1_NoIll_16S.fastq    44084     43372
#FouhaFinalJF_A12_S12_R1_NoIll_16S.fastq    56557     55714
#FouhaFinalJF_A1_S1_R1_NoIll_16S.fastq      13646     13031
#FouhaFinalJF_A2_S2_R1_NoIll_16S.fastq      48291     47546
#FouhaFinalJF_A3_S3_R1_NoIll_16S.fastq      56168     55075


errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
#https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling
#run below for pseudopool (don't need to pseudopool, only real need is if processing time was insanely long for pool
#pooling detects more rares- should just set pool=TRUE)
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool="pseudo")
#if I want to look for some specific rares, could set priors. Also explained in above link 

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#96 9795
#^ number of ASVs
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))
#188  189  190  191  192  193  194  195  196  197  198  203  207  209  210  211  212  213  214 
#2    1    4    1    2    8    1    2    2    1    1    1    1    2    1    3   16    1    1 
#217  218  222  225  226  227  228  229  230  232  234  236  238  241  243  244  246  247  248 
#1    2    1    2    1   13    5    1    1    4    2    3    1    2    2    2    4   14    8 
#249  250  251  252  253  254  255  256  257  258  259  260  261  262  263  264  265  266  267 
#4    5   11   26  226 7456  590   55   58   18    5    5    4    4    2    3    2    1    1 
#269  270  271  273  281  286  308  325  327  346  347  349 
#2    1    1    1    1    2    1    1    1    2    1   53
sum(seqtab.nochim)/sum(seqtab)
head(seqtab.nochim)
#98 percent not chimeras
dim(seqtab.nochim)
#96 8663
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)quit
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)
# input filtered denoisedF denoisedR merged nonchim
# A10_S10 56669    55786     55551     55643  55397   53715
# A11_S11 44084    43372     42348     42294  40804   39509
# A12_S12 56557    55714     53147     52805  47367   45372
# A1_S1   13646    13031     12975     12992  12802   12802
# A2_S2   48291    47546     46344     46337  44210   43777
# A3_S3   56168    55075     51206     50614  39906   37848

save.image("E:/R_stuff/Fouha/16S/16S_Fouha_ws_pseudo_pool.RData")

#load RData for everything above
load("E:/R_stuff/Fouha/16S/16S_Fouha_ws.RData")
#or with the pseudopool
load("E:/R_stuff/Fouha/16S/16S_Fouha_ws_pseudo_pool.RData")


#then moved to cluster
#settings i did on cluster
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
save.image("Post_taxa_pseudo.RData")
?assignTaxonomy()
##Phyloseq
#load rdata
setwd("~/BOSTON/Davies/Fouha Microbiome/dada2/")
#load("~/BOSTON/Davies/Fouha Microbiome/dada2/Post_taxa.RData")
#or load pseudopool
load("~/BOSTON/Davies/Fouha Microbiome/dada2/Post_taxa_pseudo.RData")
write.csv(track, file="readtrack.csv")
#load metadata
meta=read.csv(file="~/BOSTON/Davies/Fouha Microbiome/16S_metadata.csv",na.strings=c("NA","NaN", ""))
#fix metadata
meta$side <- ifelse(grepl("N", meta$Site, ignore.case = T), "North", 
                  ifelse(grepl("S", meta$Site, ignore.case = T), "South", "Other"))
meta$site=ifelse(grepl("N1", meta$Site, ignore.case = T), "1", 
                ifelse(grepl("N2", meta$Site, ignore.case = T), "2", 
                       ifelse(grepl("N3", meta$Site, ignore.case = T), "3",
                              ifelse(grepl("N4", meta$Site, ignore.case = T), "4",
                                     ifelse(grepl("S1", meta$Site, ignore.case = T), "1",
                                            ifelse(grepl("S2", meta$Site, ignore.case = T), "2",
                                                   ifelse(grepl("S3", meta$Site, ignore.case = T), "3",
                                                          ifelse(grepl("S4", meta$Site, ignore.case = T), "4","Other"))))))))
meta$position = ifelse(grepl("E", meta$Site, ignore.case = T), "Edge", 
                    ifelse(grepl("C", meta$Site, ignore.case = T), "Center", "Other"))
meta$side.site= paste(meta$side,meta$site)
new.meta <- meta[order(meta$sample),]
names=row.names(seqtab.nochim)
sort.names=sort(names)
row.names(new.meta)=sort.names


library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

##getting rid of tha contaminants 

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
###
#create fasta for all NAs
asv_tax.df=as.data.frame(asv_tax)
nas_asvs <- asv_tax.df[is.na(asv_tax.df$Class),]
nas_asvs_names=row.names(nas_asvs)
nas_indices=which(asv_fasta %in% paste0(">", nas_asvs_names))
want <- sort(c(nas_indices, nas_indices + 1))
nas_asv_fasta <- asv_fasta[want]
#the "wb" is only for windows
output.file <- file("./nas_ASVs.fa", "wb")
write(nas_asv_fasta, file=output.file)
#this following line is only necessary on windows
close(output.file)


#


###decontaim
BiocManager::install("decontam")
library(decontam)

colnames(asv_tab) # our blanks are the 4th sample in this case
vector_for_decontam <- c(rep(FALSE, 3), rep(TRUE, 1), rep(FALSE, 92))
vector_for_decontam

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 16 as contaminants

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
contam_asvs
asv_tax[row.names(asv_tax) %in% contam_asvs, ]

#[1] "ASV_153"  "ASV_155"  "ASV_191"  "ASV_232"  "ASV_647"  "ASV_775"  "ASV_824"  "ASV_868" 
#[9] "ASV_994"  "ASV_1108" "ASV_1273" "ASV_1363" "ASV_1621" "ASV_1767" "ASV_2726" "ASV_2910"
#grep -w -A1 "^>ASV_38\|^>ASV_152\|^>ASV_154\|^>ASV_193\|^>ASV_199\|^>ASV_308\|^>ASV_321\|^>ASV_337\|^>ASV_632\|^>ASV_636\|^>ASV_810\|^>ASV_915\|^>ASV_1129\|^>ASV_1326\|^>ASV_1377\|^>ASV_1706\|^>ASV_1865\|^>ASV_2221\|^>ASV_2632\" ASVs.fa


# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
test_indices=which(asv_fasta %in% paste0(">", nas_asvs_names))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)

getwd()
library("phyloseq")
library("vegan")
#BiocManager::install("DESeq2")
library("DESeq2")
library("ggplot2")
#BiocManager::install("dendextend")
library("dendextend")
#install.packages("tidyr")
library("tidyr")
library("viridis")
#install.packages("reshape")
library("reshape")

count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T,row.names = 1,
                        check.names=F, sep="\t")[ , -c(4)]

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

#sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t")

# and setting the color column to be of type "character", which helps later
#sample_info_tab$color <- as.character(sample_info_tab$color)
#in bash removed all sequences that hit eukaryote
#add lines here

#now create new lists without eukaryotes
euk_contam_asvs=read.table(file="euk_contam_asvs",fill = TRUE, stringsAsFactors = FALSE)
euk_contam_asvs=euk_contam_asvs$V1
euk_contam_asvs
# making new count table
asv_tab_no_euk_contam <- asv_tab[!row.names(count_tab) %in% euk_contam_asvs, ]

# making new taxonomy table
asv_tax_no_euk_contam <- asv_tax[!row.names(tax_tab) %in% euk_contam_asvs, ]


#Could come back and create a new fasta if I wanted to but for now lets just make new counts files

write.table(asv_tab_no_euk_contam, "ASVs_counts-no-euk-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_euk_contam, "ASVs_taxonomy-no-euk-contam.tsv",
            sep="\t", quote=F, col.names=NA)
  
  
new_count_tab <- read.table("ASVs_counts-no-euk-contam.tsv", header=T,row.names = 1,
                        check.names=F, sep="\t")[ , -c(4)]

new_tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-euk-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
  
#need to get rid of the chloroplast 
##removing chloroplast
library(stringr)

new_tax_tab_dat= data.frame(new_tax_tab)
is.chloro= rownames(new_tax_tab_dat[new_tax_tab_dat$Order == "Chloroplast",])
is.chloro=is.chloro[!str_detect(is.chloro,pattern="NA")]

is.mito= rownames(new_tax_tab_dat[new_tax_tab_dat$Family == "Mitochondria",])
is.mito=is.mito[!str_detect(is.mito,pattern="NA")]
is.mito.chloro=append(is.mito, is.chloro)

is.euk=rownames(new_tax_tab_dat[new_tax_tab_dat$Kingdom == "Eukaryota",])
is.euk=is.euk[!str_detect(is.euk,pattern="NA")]
is.mito.chloro=append(is.mito, is.chloro)
is.mito.chloro.euk=append(is.mito.chloro, is.euk)

is.arch=rownames(new_tax_tab_dat[new_tax_tab_dat$Kingdom == "Archaea",])
is.mito.chloro.euk.arch=append(is.mito.chloro.euk, is.arch)
new_count_tab.no.chloro.mito.euk.arch  = new_count_tab[!row.names(new_count_tab)%in%is.mito.chloro.euk.arch,]

#delete the neg control from the sample_info sheet
new.new.meta = new.meta[-1,]
df[match(target, df$name),]
new.new.new.meta=new.new.meta[ (colnames(new_count_tab)), ]

tax_tab_phy <- tax_table(new_tax_tab)
#tax_tab_phy

ps2 <- phyloseq(otu_table(new_count_tab.no.chloro.mito.euk.arch, taxa_are_rows=T), 
               sample_data(new.new.meta), 
               tax_table(tax_tab_phy))  



save.image("16S_PostFiltChloroMitoEukArch.RData")
load("~/BOSTON/Davies/Fouha Microbiome/dada2/16S_PostFiltChloroMitoEukArch.RData")
#add library sizes for each sample 
test=data.frame(sample_sums(ps2))

librarysize <- test[ order(row.names(test)), ]
lib.meta=cbind(librarysize, new.new.meta)
#D7 and B6 also look like outliers via pcoa as well as being the lowest lib size 
#boxplot of lib sizes
library(ggplot2)
ggplot(data = lib.meta, mapping = aes_string(x = "position",y = "librarysize"
                                            )) +
  geom_boxplot(fill = NA) +
  geom_point(size = 1, alpha = 0.3,
             position = position_jitter(width = 0.3)) +
  theme(legend.position="none")

do.call("rbind", with(lib.meta, tapply(librarysize, position,
                                          function(x) unlist(shapiro.test(x)[c("statistic", "p.value")]))))
#install.packages("RVAideMemoire")
library(RVAideMemoire)

pairwise.var.test(lib.meta$librarysize,lib.meta$position)
model=aov(librarysize ~ position, data=lib.meta)
summary(model)
TukeyHSD(model, which = "position")
#center and edge not signiciantly different

####OUTLIER AND TRIMMING UNDERREPRESENTED ASVS#####


ps2 <- phyloseq(otu_table(new_count_tab.no.chloro.mito.euk.arch, taxa_are_rows=T), 
                sample_data(lib.meta), 
                tax_table(tax_tab_phy))  

seqtab.all=data.frame(otu_table(ps2))
write.csv(seqtab.all, file= "seqtab.all.csv")
seqtab.all=read.csv("seqtab.all.csv",row.names = 1)
install.packages("MCMC.OTU")
library(MCMC.OTU)
seqtab.all.4mcmc=data.frame(t(otu_table(ps2)))
seqtab.all.4mcmc=add_column(seqtab.all.4mcmc, sample = row.names(seqtab.all.4mcmc), .after = 0)
seqtab.all.trim=purgeOutliers(seqtab.all.4mcmc,count.columns=2:7654,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)

#[1] "samples with counts below z-score -2.5 :"
#[1] "B6_S18" "C4_S28" "D7_S43" "E4_S52"
#[1] "zscores:"
#B6_S18    C4_S28    D7_S43    E4_S52 
#-3.067871 -2.948323 -3.919400 -2.814662 
#[1] "OTUs passing frequency cutoff  1e-04 : 890"
#[1] "OTUs with counts in 0.02 of samples:"

#FALSE  TRUE 
#66   824 

#rarefying 
library(vegan)
rarecurve(t(otu_table(ps2)), step=50, cex=0.005)

#ps2noout=subset_samples(ps2, sample!="B6"& sample!="D7"& sample!="C4"& sample!="E4")
#rarecurve(t(otu_table(ps2noout)), step=50, cex=0.005)
#ps.rarefied.all=rarefy_even_depth(ps2noout, rngseed=1, sample.size=5000, replace=F)

#edge
rarecurve(t(otu_table(onlyedgepsallsites)), step=50, cex=0.005)

#center
rarecurve(t(otu_table(onlycenterpsallsites)), step=50, cex=0.005)

#sed #for this I think it makes sense to rarefy to minimum 
rarecurve(t(otu_table(onlysedpsallsites)), step=50, cex=0.005)

seqtab.all.rare <- rrarefy(seqtab.all,sample=5000)
rarecurve(seqtab.all,step=100,label=FALSE)
rarecurve(seqtab.all.rare,step=100,label=FALSE)
seqtab.all.rare.4mcmc=data.frame(t(seqtab.all.rare))
seqtab.all.rare.4mcmc=add_column(seqtab.all.rare.4mcmc, sample = row.names(seqtab.all.rare.4mcmc), .after = 0)
#remove outliers 
rownames(seqtab.all.rare.4mcmc)
remove.rows=c("B6_S18","C4_S28","D7_S43","E4_S52")
seqtab.all.rare.4mcmc=seqtab.all.rare.4mcmc[!(row.names(seqtab.all.rare.4mcmc) %in% remove.rows), ]
seqtab.all.rare.trim=purgeOutliers(seqtab.all.rare.4mcmc,count.columns=2:7654,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#[1] "samples with counts below z-score -2.5 :"
#character(0)
#[1] "zscores:"
#named numeric(0)
#[1] "OTUs passing frequency cutoff  1e-04 : 1736"
#[1] "OTUs with counts in 0.02 of samples:"
#FALSE  TRUE 
#188  1548 
#^ Rarefied doesnt find the same 4 outliers- should I just remove beforehand then?YEAH 



new.new.new.meta=new.new.meta[ (rownames(seqtab.all.trim)), ]
ps.all.trim <- phyloseq(otu_table(data.matrix(seqtab.all.trim), taxa_are_rows=F), 
                sample_data(new.new.new.meta), 
                tax_table(tax_tab_phy))  

onlyedgepsall.trim=subset_samples(ps.all.trim, position=="Edge")
onlycenterpsall.trim=subset_samples(ps.all.trim, position=="Center")
onlysedpsall.trim=subset_samples(ps.all.trim, position=="Other")

ps.all.rare.trim=phyloseq(otu_table(data.matrix(seqtab.all.rare.trim), taxa_are_rows=F), 
                          sample_data(new.new.new.meta), 
                          tax_table(tax_tab_phy)) 
onlyedgeps.rare.trim=subset_samples(ps.all.rare.trim, position=="Edge")
onlycenterps.rare.trim=subset_samples(ps.all.rare.trim, position=="Center")
onlysedps.rare.trim=subset_samples(ps.all.rare.trim, position=="Other")

save.image("16S_PostTrim&Rarefy.RData")
setwd("~/BOSTON/Davies/Fouha Microbiome/dada2/")
load("~/BOSTON/Davies/Fouha Microbiome/dada2/16S_PostTrim&Rarefy.RData")
write.csv(new.new.new.meta,"new.new.new.meta.csv")
###Creating input files for picrust###
# making new fasta file
library("Biostrings")

s = readDNAStringSet("ASVs.fa")
all.rare.trim.asvs=colnames(seqtab.all.rare.trim)
length(all.rare.trim.asvs)
selected_sequences <- s[all.rare.trim.asvs]
writeXStringSet(selected_sequences, "all.rare.trim.fa")
#biom format, i dont think neceessary can just make .tsv below 
# BiocManager::install("biomformat")
# library(biomformat)
# seqtab_biom <- t(seqtab.all.rare.trim)
# biom_object <- make_biom(data = seqtab_biom)
# write_biom(biom_object, biom_file = "all.rare.trim.biom")

head(seqtab_biom)


# count table:
asv_tab <- t(seqtab.all.rare.trim)
write.table(asv_tab, "ASVs_counts.all.rare.trim.tsv", sep="\t", quote=F)


