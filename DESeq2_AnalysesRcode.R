#setwd("~/MCPetro_IEGSeq")

#Three input files are required:

###gene_count_matrix+1.csv
###mcpetro_phenodataDESeq.csv
###MTxIPxSpecies_SigCounts.csv

#Load & install required packages
library(readr)
library(ggplot2)
library(viridis)
library("pheatmap")
library("RColorBrewer")
library(Rtsne)

#library(Biobase)
#library("DESeq2")

#Please note: Petro instead of PC used in all code for Petrotilapia chtimba (rock species). MC used for Mchenga conophoros (sand species).

#DESeq2 Analyses
#This is the read count matrix which is the input for DESeq2. 
#This was produced by StringTie, followed by a python script (prepDE.py) provided here http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual. 
#Then a pseudocount of 1 was added to all cells.

countData <- as.matrix(read.csv("gene_count_matrix+1.csv", row.names="gene_id"))

#This is the phenotype data.

colData <- read.csv("mcpetro_phenodataDESeq.csv", sep=",", row.names=1)
colData$sampletype<-as.factor(colData$sampletype)
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

#Test effect of sampletype controlling for species & test:    
dds1 <- DESeqDataSetFromMatrix(countData = countData, 
                               colData = colData, design = ~ batch + species + test + sampletype)

dds1 <- DESeq(dds1,test="LRT", reduced= ~ batch + species + test)
res1 <- results(dds1)
summary(res1)
resultsNames(dds1)

boxplot(log10(assays(dds1)[["cooks"]]), range=0, las=2)
plotDispEsts(dds1)

(res1Ordered <- res1[order(res1$padj), ])

res1$padj[is.na(res1$padj)] <- 100
res1Sig <- as.data.frame(res1[res1$padj < 0.05,])

#save(dds1,res1,res1Sig, colData,file="IPvsInputOnly.RData") # save it for the future

#load("~/IPvsInputOnly.RData")

#write.table(res1Sig, file="IP vs Input Only Significant.csv",sep=",",row.names=T)

#Volcano Plot
res1_table <- as.data.frame(res1)
#write.table(res1_table, file="Model1_InputvIP_DESeq2_output.csv",sep=",",row.names=T)

##Highlight genes that have an absolute fold change > 2 and adjusted p-value < 0.05
res1_table$threshold = as.factor(abs(res1_table$log2FoldChange) > 2 & res1_table$padj < 0.05)

#png("InputvsIPVolcanoPlot.png", units="in", width=5, height=4, res=300)
ggplot(data=res1_table, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  scale_colour_manual(values=c("slategray3","darkorchid2")) +
  theme(legend.position = "none",axis.text.x=element_text(size=14, color='black')
        ,axis.text.y=element_text(size=18, color='black')
        ,text = element_text(size=24, color='black')) +
  geom_hline(yintercept = 1.301, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  xlim(c(-10, 10)) + ylim(c(0, 20)) +
  xlab(expression(log[2]("Fold Change"))) + ylab(expression(log[10]("p-value")))
#dev.off()

#Test two-way interaction of test x samples, controlling for batch, species, and test:    
dds2 <- DESeqDataSetFromMatrix(countData = countData, 
                               colData = colData, design = ~ batch + sampletype + test + species + sampletype:test)

dds2 <- DESeq(dds2, test="LRT", reduced = ~ batch + sampletype + test + species)
res2 <- results(dds2)
summary(res2)
matrix(resultsNames(dds2))

boxplot(log10(assays(dds2)[["cooks"]]), range=0, las=2)
plotDispEsts(dds2)

(res2Ordered <- res2[order(res2$padj), ])

res2$padj[is.na(res2$padj)] <- 100
res2Sig <- as.data.frame(res2[res2$padj < 0.05,])

res2_table <- as.data.frame(res2)
#write.table(res2_table, file="Model2_MTxIP_DESeq2_output.csv",sep=",",row.names=T)

#save(dds2,res2,res2Sig, colData,file="MTxIP.RData") # save it for the future 

#load("~/MTxIP.RData")

#write.table(res2Sig, file="MTxIP Significant.csv",sep=",",row.names=T)

#Volcano Plots

##Highlight genes that have an absolute fold change > 2 and adjusted p-value < 0.05
res2_table$threshold = as.factor(abs(res2_table$log2FoldChange) > 2 & res2_table$padj < 0.05)

#png("MTxIPVolcanoPlot.tiff", units="in", width=5, height=4, res=300)
ggplot(data=res2_table, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  scale_colour_manual(values=c("slategray3","darkorchid2")) +
  theme(legend.position = "none",axis.text.x=element_text(size=14, color='black')
        ,axis.text.y=element_text(size=18, color='black')
        ,text = element_text(size=24, color='black')) +
  geom_hline(yintercept = 1.301, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  xlim(c(-40, 40)) + ylim(c(0, 5)) +
  xlab(expression(log[2]("Fold Change"))) + ylab(expression(log[10]("p-value")))
#dev.off()

#Test for effect of full 3-way interaction genes that show differential fold enrichment between treatment and species:
dds3 <- DESeqDataSetFromMatrix(countData = countData, 
                               colData = colData, design = ~ batch + sampletype + test + species + sampletype:species + test:species + sampletype:test + sampletype:test:species)

dds3 <- DESeq(dds3, test="LRT", reduced= ~ batch + sampletype + test + species + sampletype:species + test:species + sampletype:test)
res3 <- results(dds3)
summary(res3)
matrix(resultsNames(dds3))

boxplot(log10(assays(dds3)[["cooks"]]), range=0, las=2)
plotDispEsts(dds3)

(res3Ordered <- res3[order(res3$padj), ])
res3$padj[is.na(res3$padj)] <- 100
res3Sig <- as.data.frame(res3[res3$padj < 0.05,])

res3_table <- as.data.frame(res3)
#write.table(res3_table, file="Model3_MTxIPxSpecies_DESeq2_output.csv",sep=",",row.names=T)

#save(dds3,res3,res3Sig, colData,file="MTxIPxSpecies.RData") # save it for the future 

#load("~/MTxIPxSpecies.RData")

#write.table(res3, file="MTxIPxSpecies_out.csv",sep=",",row.names=T)
#write.table(res3Sig, file="MTxIPxSpecies Significant.csv",sep=",",row.names=T)

#Volcano Plots

##Highlight genes that have an absolute fold change > 2 and adjusted p-value < 0.05
res3_table$threshold = as.factor(abs(res3_table$log2FoldChange) > 2 & res3_table$padj < 0.05)

#png("MTxIPxSpeciesVolcanoPlot.tiff", units="in", width=5, height=4, res=300)
ggplot(data=res3_table, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  scale_colour_manual(values=c("slategray3","darkorchid2")) +
  theme(legend.position = "none",axis.text.x=element_text(size=14, color='black')
        ,axis.text.y=element_text(size=18, color='black')
        ,text = element_text(size=24, color='black')) +
  geom_hline(yintercept = 1.301, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  xlim(c(-40, 40)) + ylim(c(0, 5)) +
  xlab(expression(log[2]("Fold Change"))) + ylab(expression(log[10]("p-value")))
#dev.off()


##Fold Enrichment (FE) Analyses

#Reload countdata as data frame to allow for manipulation & fold enrichment calculations.
countData <- as.data.frame(read.csv("gene_count_matrix+1.csv", row.names="gene_id"))
colData <- read.csv("mcpetro_phenodataDESeq.csv", sep=",", row.names=1)

attach(countData)

countData$MC10_FE <-(MC10.IP/MC10.Input)
countData$MC21_FE <-(MC21.IP/MC21.Input)
countData$MC22_FE <-(MC22.IP/MC22.Input)
countData$MC9_FE <-(MC9.IP/MC9.Input)
countData$Petro11_FE <-(Petro11.IP/Petro11.Input)
countData$Petro12_FE <-(Petro12.IP/Petro12.Input)
countData$Petro23_FE <-(Petro23.IP/Petro23.Input)
countData$Petro24_FE <-(Petro24.IP/Petro24.Input)

attach(countData)

countData$MCControl_FE<-(MC10_FE+MC22_FE)/2
countData$MCMT_FE<-(MC9_FE+MC21_FE)/2
countData$PetroControl_FE<-(Petro11_FE+Petro24_FE)/2
countData$PetroMT_FE<-(Petro12_FE+Petro23_FE)/2

countData$Control_FE<-(MC10_FE+MC22_FE+Petro11_FE+Petro24_FE)/4
countData$MT_FE<-(MC9_FE+MC21_FE+Petro12_FE+Petro23_FE)/4

attach(countData)

countData$NonPetroMT_FE<-(MCControl_FE+MCMT_FE+PetroControl_FE)/3
countData$NonMCMT_FE<-(MCControl_FE+PetroMT_FE+PetroControl_FE)/3

countData$MTDiff_FE<-(MT_FE-Control_FE)
countData$MTDiff_FE_bin<-ifelse(abs(countData$MTDiff_FE)>2,countData$MTDiff_FE_bin<-1,countData$MTDiff_FE_bin<-0)

countData$MTPetroDiff_FE<-(PetroMT_FE-(MCControl_FE+MCMT_FE+PetroControl_FE)/3)
countData$MTMCDiff_FE<-(MCMT_FE-(MCControl_FE+PetroMT_FE+PetroControl_FE)/3)

countData$nonMTPetroDiff_FE<-(PetroControl_FE-(MCControl_FE+MCMT_FE+PetroMT_FE)/3)
countData$nonMTMCDiff_FE<-(MCControl_FE-(MCMT_FE+PetroMT_FE+PetroControl_FE)/3)

countData$PetroMTvControl_FE<-(PetroMT_FE)/(PetroControl_FE)
countData$MCMTvControl_FE<-(MCMT_FE)/(MCControl_FE)

#load("~/MTxIP.RData")

res2<-as.data.frame(res2)

res2$genename<-rownames(res2)
countData$genename<-rownames(countData)

total2 <- merge(countData,res2,by="genename")

total2$res2sig<-ifelse(total2$padj < 0.05,total2$res2sig<-1,total2$res2sig<-0)
total2$res2sig<-as.factor(total2$res2sig)
total2$res2sig<-factor(total2$res2sig,levels=c(0,1), labels=c("Not significant","p < 0.05"))

total2.ordered <- total2[order(total2$res2sig),] 
label(total2.ordered$res2sig) <- "Significance" 

#tiff("MTvControlFEplot.tiff", units="in", width=6, height=6, res=300)
ggplot(total2.ordered, aes(x=log10(Control_FE), y=log10(MT_FE), shape=res2sig, alpha=res2sig, color=log2FoldChange)) +
  geom_point()+
  xlim(-5, 5)+ylim(-5, 5)+
  scale_shape_manual(values=c(1, 16), name="Significance")+
  scale_alpha_manual(values=c(0.1,1), guide="none")+
  scale_color_viridis()+
  geom_vline(xintercept = 0, linetype="dashed", color = "black")+
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  geom_abline(intercept=0, linetype="dashed", color = "gray")+
  theme_classic()+
  theme(legend.position = c(0.15,0.775), text = element_text(size=12, color='black'))
#dev.off()

total2.ordered <- total2[order(total2$padj),] 

#Individual Gene Count Plots

# MSTRG.39541 (HSD 17-B 8)
genecounts <- plotCounts(dds2, gene = "MSTRG.39541", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("HSD17B8counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("HSD17B8") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

# MSTRG.8258 (GnRH2)
genecounts <- plotCounts(dds2, gene = "MSTRG.8258", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("GNRH2counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("GNRH2") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

# MSTRG.6317 (CACNA1H)
genecounts <- plotCounts(dds2, gene = "MSTRG.6317", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("CACNA1Hcounts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("CACNA1H") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.42512 (PTGFRN) 
genecounts <- plotCounts(dds2, gene = "MSTRG.42512", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("PTGFRNcounts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("PTGFRN") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.42558 (POMC)
genecounts <- plotCounts(dds3, gene = "MSTRG.42558", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("POMCcounts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("POMC") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.40407 (NPY)
genecounts <- plotCounts(dds3, gene = "MSTRG.40407", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("NPYcounts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("NPY") +
  theme(plot.title = element_text(hjust = 0.5))
  #+ theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.2914 (NR3C2 (Melanocortin Receptor))
genecounts <- plotCounts(dds3, gene = "MSTRG.2914", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("NR3C2counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("NR3C2 (Melanocortin Receptor)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.28953 (GRM1, mGluR1)
genecounts <- plotCounts(dds3, gene = "MSTRG.28953", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("GRM1counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("GRM1") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.16448 (GPER1)
genecounts <- plotCounts(dds3, gene = "MSTRG.16448", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("GPER1counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("GPER1") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.41018 (SRD5A1)
genecounts <- plotCounts(dds3, gene = "MSTRG.41018", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("SRD5A1counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("SRD5A1") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.41719 (SRD5A3)
genecounts <- plotCounts(dds3, gene = "MSTRG.41719", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("SRD5A3counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("SRD5A3") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.10190 (LCOR)
genecounts <- plotCounts(dds3, gene = "MSTRG.10190", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("LCORcounts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("LCOR") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.5356 (NLGN2)
genecounts <- plotCounts(dds3, gene = "MSTRG.5356", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("NLGN2counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("NLGN2") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.38619 (OGFR)
genecounts <- plotCounts(dds3, gene = "MSTRG.38619", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)
#tiff("OGFRcounts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("OGFR") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#MSTRG.49345 (PCDH1)
genecounts <- plotCounts(dds3, gene = "MSTRG.49345", 
                         intgroup = c("sampletype","test","species"),returnData = TRUE)

#tiff("PCDH1counts.tiff", units="in", width=2, height=2.5, res=300)
ggplot(genecounts, aes(x=sampletype, y=count,fill=test,shape=species)) +
  geom_point(aes(shape=species),position=position_jitter(w=0.025,h=0),size=3) + 
  xlab("") + # for the x axis label
  scale_fill_manual(values=c("midnightblue", "goldenrod"))+scale_shape_manual(values=c(21, 22))+scale_y_log10()+theme_bw() +
  ggtitle("PCDH1") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none", text = element_text(size=14, color='black'))
#dev.off()

#Fold enrichment for Model 3 (MTxIPxSpecies)
res3<-as.data.frame(res3)
res3$genename<-rownames(res3)
countData$genename<-rownames(countData)

total3 <- merge(countData,res3,by="genename")

total3$res3sig<-ifelse(total3$padj < 0.05,total3$res3sig<-1,total3$res3sig<-0)
total3$res3sig<-as.factor(total3$res3sig)
total3$res3sig<-factor(total3$res3sig,levels=c(0,1), labels=c("Not significant","p < 0.05"))

total3.ordered <- total3[order(total3$res3sig),] 

#tiff("MTxPetrovAllElseFEplot.tiff", units="in", width=6, height=6, res=300)
ggplot(total3.ordered, aes(x=log10(PetroMT_FE), y=log10(MCMT_FE), shape=res3sig, alpha=res3sig, color=log2FoldChange)) +
  geom_point()+
  xlim(-5, 5)+ylim(-5, 5)+
  scale_shape_manual(values=c(1, 16), name="Significance")+
  scale_alpha_manual(values=c(0.1,1), guide="none")+
  scale_color_viridis()+
  geom_vline(xintercept = 0, linetype="dashed", color = "black")+
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  geom_abline(intercept=0, linetype="dashed", color = "gray")+
  theme(legend.position = c(0.05,0.775), text = element_text(size=12, color='black'))
#dev.off()

#tiff("MTxPetrovAllElseFEplot.tiff", units="in", width=6, height=6, res=300)
ggplot(total3.ordered, aes(x=log10(PetroControl_FE), y=log10(MCControl_FE), shape=res3sig, alpha=res3sig, color=log2FoldChange)) +
  geom_point()+
  xlim(-5, 5)+ylim(-5, 5)+
  scale_shape_manual(values=c(1, 16), name="Significance")+
  scale_alpha_manual(values=c(0.1,1), guide="none")+
  scale_color_viridis()+
  geom_vline(xintercept = 0, linetype="dashed", color = "black")+
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  geom_abline(intercept=0, linetype="dashed", color = "gray")+
  theme(legend.position = c(0.05,0.775), text = element_text(size=12, color='black'))
#dev.off()

###Sample Distances Plot
rld3 <- rlog(dds3, blind = FALSE)
sampleDists <- dist(t(assay(rld3)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld3$dex, rld3$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

##################################################################################

#Making paper figures

setwd("~/Dropbox/My Documents/Manuscripts/Aggression & PhosphoTRAP cichlid paper/Images")

#Figure 3B
#png("InputvsIPVolcanoPlot.tiff", units="in", width=2.56, height=2.56, res=600)
ggplot(data=res1_table, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.1, size=1) +
  scale_colour_manual(values=c("slategray3","darkorchid2")) +
  theme(plot.background = element_blank()
        ,panel.background = element_blank()
        #,panel.border = element_rect(colour = "black", fill=NA)
        ,legend.position = "none",axis.text.x=element_text(size=11, color='black')
        ,axis.text.y=element_text(size=11, color='black')
        ,axis.line = element_line(colour = "black")
        ,text = element_text(size=12, color='black')) +
  geom_hline(yintercept = 1.301, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  xlim(c(-8, 8)) + ylim(c(0, 20)) +
  xlab(expression(log[2]("Fold Change"))) + ylab(expression(log[10]("p-value")))
#dev.off()

#Principle Components Analysis####################

#Figure 3A
#rlog PCS#########################################
rld3 <- rlog(dds3, blind = FALSE)
head(assay(rld3), 3)
pcaDatarlog<-plotPCA(rld3, intgroup = c("sampletype", "test", "species","batch"), returnData = TRUE)
pcaDatarlog$TestSample<-paste(pcaDatarlog$sampletype,pcaDatarlog$test,sep=":")
percentVar <- round(100 * attr(pcaDatarlog, "percentVar"))
pcaDatarlog
#write.table(pcaDatarlog, file="pcaDatarlog.csv",sep=",",row.names=T)

#png("pcaPlotrlog.tiff", units="in", width=4.13, height=2.56, res=600)
ggplot(pcaDatarlog, aes(x = PC1, y = PC2, color = species, shape = TestSample)) +
  geom_point(size = 2) +
  scale_shape_manual(name="Condition",values=c(1,2,16,17)) +
  scale_fill_manual(name="Species", values=c("midnightblue","goldenrod1"))+
  scale_colour_manual(name="Species",values=c("midnightblue","goldenrod1")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+
  theme_bw() +
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text = element_text(size=10)
        #,legend.position=c(.8,.65)
        ,legend.background = element_rect(fill="white", size=.5, linetype="solid")
        ,legend.margin=unit(0.5,"cm")
        #,legend.key.height=unit(0.75,"line")
        ,legend.key = element_blank())
#dev.off()

#Figure 4
#tiff("MTvControlFEplot.tiff", units="in", width=5.12, height=3.15, res=600)
ggplot(total2.ordered, aes(x=log10(Control_FE), y=log10(MT_FE), shape=res2sig, alpha=res2sig, color=log2FoldChange)) +
  geom_point()+
  xlim(-3.75, 3.75)+ylim(-3.75, 3.75)+
  scale_shape_manual(values=c(1, 16), name="Significance")+
  scale_alpha_manual(values=c(0.1,1), guide="none")+
  scale_color_viridis()+
  geom_vline(xintercept = 0, linetype="dashed", color = "black")+
  geom_hline(yintercept = 0, linetype="dashed", color = "black")+
  geom_abline(intercept=0, linetype="dashed", color = "gray")+
  theme_bw()+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,legend.margin=unit(0.0,"cm")
        ,text = element_text(size=10))
#dev.off()

###tSNE plot ###########################################
# From: https://www.analyticsvidhya.com/blog/2017/01/t-sne-implementation-r-python/

#Figure 5

threewaySigGenes <- read_csv("MTxIPxSpecies_SigCounts.csv")
Labels<-threewaySigGenes$gene_id
threewaySigGenes$gene_id<-as.factor(threewaySigGenes$gene_id)

## for plotting
#colors = viridis(length(unique(threewaySigGenes$gene_id)))
#names(colors) = unique(threewaySigGenes$gene_id)

## Executing the algorithm on curated data
tsne <- Rtsne(threewaySigGenes[,-1], dims = 2, perplexity=20, verbose=TRUE, max_iter = 1000)
#exeTimeTsne<- system.time(Rtsne(threewaySigGenes[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500))

## Plotting
tSNEplot<-data.frame(tsne$Y)

tSNEplot$gene_id <- threewaySigGenes$gene_id
tSNEplot<-tSNEplot[c(3,1,2)]

gene_id <- rownames(res3Sig)
rownames(res3Sig) <- NULL
res3Sig <- cbind(gene_id,res3Sig)

mergedtSNE <- merge(tSNEplot, res3Sig, by = c('gene_id'))

gene_id <- rownames(countData)
rownames(countData) <- NULL
countData <- cbind(gene_id,countData)
mergedtSNE2 <- merge(mergedtSNE, countData, by = c('gene_id'))

#write.table(mergedtSNE2, file="mergedtSNE.csv",sep=",",row.names=T)

#Note: tSNE is a stochastic process, so this final figure will differ from the the published figure.

#Figure 5a

#tiff("tSNELog2FoldChange.tiff", units="in", width=4.68, height=4,68, res=600)
ggplot(mergedtSNE2,aes(X1,X2,color = as.numeric(log2FoldChange))) + geom_point(size = 2, alpha = 0.6) +
  theme_bw()+ scale_color_viridis() +
#  xlim(-35, 35)+ ylim(-35, 35) +
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,legend.position=c(0.85, 0.1)
        ,legend.direction = "horizontal"
        ,legend.title=element_blank()
        ,text = element_text(size=10))
#dev.off()

#Figure 5b
#tiff("tSNELog10baseMean.tiff", units="in", width=2, height=1.56, res=600)
ggplot(mergedtSNE2,aes(X1,X2,colour = as.numeric(log10(baseMean)))) + geom_point(size = .5, alpha = 0.6) +
  theme_bw()+ scale_color_viridis()+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,axis.text.x = element_blank()
        ,axis.text.y = element_blank()
        ,axis.text = element_blank()
        ,axis.title = element_blank()
        ,legend.position="right"
        ,legend.direction = "vertical"
        ,legend.title=element_blank()
        ,legend.margin=margin(0,0,0,0)
        ,legend.box.margin=margin(0,0,0,0)
        ,text = element_text(size=8))
#dev.off()

#Figure 5c
#tiff("tSNELog10PetroMTvControl_FE.tiff", units="in", width=2, height=1.56, res=600)
ggplot(mergedtSNE2,aes(X1,X2,colour = as.numeric(log10(PetroMTvControl_FE)))) + geom_point(size = .5, alpha = 0.6) +
  theme_bw()+ scale_color_viridis()+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,axis.text.x = element_blank()
        ,axis.text.y = element_blank()
        ,axis.text = element_blank()
        ,axis.title = element_blank()
        ,legend.position="right"
        ,legend.direction = "vertical"
        ,legend.title=element_blank()
        ,legend.margin=margin(0,0,0,0)
        ,legend.box.margin=margin(0,0,0,0)
        ,text = element_text(size=8))
#dev.off()

#Figure 5d
#tiff("tSNELog10MCMTvControl_FE.tiff", units="in", width=2, height=1.56, res=600)
ggplot(mergedtSNE2,aes(X1,X2,colour = as.numeric(log10(MCMTvControl_FE)))) + geom_point(size = .5, alpha = 0.6) +
  theme_bw()+ scale_color_viridis()+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,axis.text.x = element_blank()
        ,axis.text.y = element_blank()
        ,axis.text = element_blank()
        ,axis.title = element_blank()
        ,legend.position="right"
        ,legend.direction = "vertical"
        ,legend.title=element_blank()
        ,legend.margin=margin(0,0,0,0)
        ,legend.box.margin=margin(0,0,0,0)
        ,text = element_text(size=8))
#dev.off()
