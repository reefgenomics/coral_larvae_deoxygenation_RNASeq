setwd()
library(DESeq2)
library(ggplot2)

map.pre=read.delim("./Input_files/Metadata_ASL_counts.txt", row.names = 1, header = TRUE)
map=map.pre[order(rownames(map.pre)), ]
rownames(map)=gsub("\\.", "-", rownames(map))
map.a=map[order(rownames(map)), ] # change from A to B when using lar
map2=subset(map, map$Time == "T1") 
map4=subset(map, map$Time == "T2") 
map6=subset(map, map$Time == "T3")
#comparing T3C vs T2C for gene upregulation
map20=subset(map.a, map.a$Time == "T2" & map.a$Condition == "Control"| map.a$Time == "T3" & map.a$Condition == "Control") 

cnt_S=read.table("./Input_files/FinalTable_counts_ASL.txt", header=TRUE, row.names = 1, sep="\t")
cnts=round(cnt_S, 0)
cnts.sub=subset(cnts, rowSums(cnts) > 50)
cnts.sor=cnts.sub[ , order(names(cnts.sub))]
colnames(cnts.sor)=gsub("\\.", "-", colnames(cnts.sor))
cnt2=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map2))])]
cnt4=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map4))])]
cnt6=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map6))])]
cnt20=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map20))])]

#comparison 2 LARVAE
dds1=DESeqDataSetFromMatrix(countData = as.matrix(cnt2), colData = map2, design = ~Condition)
dds1= DESeq(dds1)
res1= results(dds1, contrast=c("Condition","Deoxygenation","Control"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res1, alpha=0.05), "comparison 2: Larvae T1, treatment vs control")
results1=as.data.frame(subset(res1, res1$padj<0.05))
write.table(results1, "output1.txt", sep = "\t", quote = FALSE, row.names= TRUE)

#comparison 4 LARVAE
dds1=DESeqDataSetFromMatrix(countData = as.matrix(cnt4), colData = map4, design = ~Condition)
dds1= DESeq(dds1)
res1= results(dds1, contrast=c("Condition","Deoxygenation","Control"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res1, alpha=0.05), "comparison 4: Larvae T1, treatment vs control")
results1=as.data.frame(subset(res1, res1$padj<0.05))
write.table(results1, "output2.txt", sep = "\t", quote = FALSE, row.names= TRUE)

#comparison 6 LARVAE
dds1=DESeqDataSetFromMatrix(countData = as.matrix(cnt6), colData = map6, design = ~Condition)
dds1= DESeq(dds1)
res1= results(dds1, contrast=c("Condition","Deoxygenation","Control"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res1, alpha=0.05), "comparison 6: Larvae T1, treatment vs control")
results1=as.data.frame(subset(res1, res1$padj<0.05))
write.table(results1, "output3.txt", sep = "\t", quote = FALSE, row.names= TRUE)

#comparison 20
dds13=DESeqDataSetFromMatrix(countData = as.matrix(cnt20), colData = map20, design = ~Time)
dds13= DESeq(dds13)
res13= results(dds13, contrast=c("Time","T3","T2"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res13, alpha=0.05), "comparison 20: Larvae T3c vs T2c")
results13=as.data.frame(subset(res13, res13$padj<0.05))
write.table(results13, "output4.txt", sep = "\t", quote = FALSE, row.names= TRUE)
