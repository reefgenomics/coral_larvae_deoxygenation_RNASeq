setwd()
#selecting genes of interest from fpkm table and average/SE for groups

cnt=read.table("./Input_file/FinalTable_fpkms_ASL.txt", header=TRUE, row.names = 1, sep= "\t")
cnts.sor=cnt[ , order(names(cnt))]
colnames(cnts.sor)=gsub("\\.", "-", colnames(cnts.sor))
rownames(cnts.sor)=gsub("\\-RA", "", rownames(cnts.sor))

map.pre=read.delim("./Input_file/Metadata_ASL_fpkms.txt", row.names = 1, header = TRUE)
map=map.pre[order(rownames(map.pre)), ]
map.C=subset(map, map$Time == "T1" & map$Condition == "Control")
map.D=subset(map, map$Time == "T1" & map$Condition == "Deoxygenation")
map.E=subset(map, map$Time == "T2" & map$Condition == "Control")
map.F=subset(map, map$Time == "T2" & map$Condition == "Deoxygenation")
map.G=subset(map, map$Time == "T3" & map$Condition == "Control")
map.H=subset(map, map$Time == "T3" & map$Condition == "Deoxygenation")

gene=read.table("./Input_file/goi_list_ASL.txt", header=TRUE, row.names= 1, sep= "\t")
genes=subset(gene, gene$anno == "HIFa") #change gene accordingly

#T1 control
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.C))])]
cnts.a=subset(cnts, row.names(cnts) %in% row.names(genes))
cnts.a$Mean= rowMeans(cnts.a, na.rm = FALSE)
cnts.b=transform(cnts.a, SD=apply(cnts.a,1, sd, na.rm = TRUE))
x <- cnts.b$SD
cnts.a$SE= x/sqrt(ncol(cnts))
cnts.a$Timepoint = rep("T1",length(nrow(cnts.a)))
cnts.a$Condition = rep("Control",length(nrow(cnts.a)))
T1_C= cnts.a[,c(5:8)]

#T1 treatment
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.D))])]
cnts.a=subset(cnts, row.names(cnts) %in% row.names(genes))
cnts.a$Mean= rowMeans(cnts.a, na.rm = FALSE)
cnts.b=transform(cnts.a, SD=apply(cnts.a,1, sd, na.rm = TRUE))
x <- cnts.b$SD
cnts.a$SE= x/sqrt(ncol(cnts))
cnts.a$Timepoint = rep("T1",length(nrow(cnts.a)))
cnts.a$Condition = rep("Treatment",length(nrow(cnts.a)))
T1_T= cnts.a[,c(5:8)]

#T2 control
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.E))])]
cnts.a=subset(cnts, row.names(cnts) %in% row.names(genes))
cnts.a$Mean= rowMeans(cnts.a, na.rm = FALSE)
cnts.b=transform(cnts.a, SD=apply(cnts.a,1, sd, na.rm = TRUE))
x <- cnts.b$SD
cnts.a$SE= x/sqrt(ncol(cnts))
cnts.a$Timepoint = rep("T2",length(nrow(cnts.a)))
cnts.a$Condition = rep("Control",length(nrow(cnts.a)))
T2_C= cnts.a[,c(5:8)]

#T2 treatment
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.F))])]
cnts.a=subset(cnts, row.names(cnts) %in% row.names(genes))
cnts.a$Mean= rowMeans(cnts.a, na.rm = FALSE)
cnts.b=transform(cnts.a, SD=apply(cnts.a,1, sd, na.rm = TRUE))
x <- cnts.b$SD
cnts.a$SE= x/sqrt(ncol(cnts))
cnts.a$Timepoint = rep("T2",length(nrow(cnts.a)))
cnts.a$Condition = rep("Treatment",length(nrow(cnts.a)))
T2_T= cnts.a[,c(5:8)]

#T3 control
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.G))])]
cnts.a=subset(cnts, row.names(cnts) %in% row.names(genes))
cnts.a$Mean= rowMeans(cnts.a, na.rm = FALSE)
cnts.b=transform(cnts.a, SD=apply(cnts.a,1, sd, na.rm = TRUE))
x <- cnts.b$SD
cnts.a$SE= x/sqrt(ncol(cnts))
cnts.a$Timepoint = rep("T3",length(nrow(cnts.a)))
cnts.a$Condition = rep("Control",length(nrow(cnts.a)))
T3_C= cnts.a[,c(5:8)]

#T3 treatment
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.H))])]
cnts.a=subset(cnts, row.names(cnts) %in% row.names(genes))
cnts.a$Mean= rowMeans(cnts.a, na.rm = FALSE)
cnts.b=transform(cnts.a, SD=apply(cnts.a,1, sd, na.rm = TRUE))
x <- cnts.b$SD
cnts.a$SE= x/sqrt(ncol(cnts))
cnts.a$Timepoint = rep("T3",length(nrow(cnts.a)))
cnts.a$Condition = rep("Treatment",length(nrow(cnts.a)))
T3_T= cnts.a[,c(5:8)]

out=rbind(T1_T, T1_C, T2_T, T2_C, T3_T, T3_C)
out$Type = rep("hifa",length(nrow(out))) #change accordingly
write.table(out, "HIFa.txt", sep = "\t", quote = FALSE, row.names= FALSE, col.names = TRUE) #change accordingly

#controls only
out=rbind(T1_C, T2_C, T3_C)
out$Type = rep("HIFa",length(nrow(out))) #change accordingly
write.table(out, "HIFa_control.txt", sep = "\t", quote = FALSE, row.names= FALSE, col.names = TRUE) #change accordingly

#treatment only
out=rbind(T1_T, T2_T, T3_T)
out$Type = rep("HIFa",length(nrow(out))) #change accordingly
write.table(out, "HIFa_treatment.txt", sep = "\t", quote = FALSE, row.names= FALSE, col.names = TRUE) #change accordingly


library("ggplot2")
#scatterplots of fpkm expression patterns over timepoints

data <- read.delim("./Input_file/HIFa.txt") #change file accordingly for different genes
data$Timepoint <- as.character(data$Timepoint) #as.character() identifies as character then can be made into a factor
data$Timepoint <- factor(data$Timepoint, levels=unique(data$Timepoint))

#HIFa gene as example here, change scale/colour accordingly 
ggplot(data, aes(x = Timepoint, y = Mean, color = Condition, group = Condition)) + ylim(10, 50) + geom_point(aes(fill = Condition), shape=21, colour= "black", size=7) + geom_line(size=2) + ylab(label="FPKMs") + xlab("Time point") + theme_classic(base_size = 29) + theme(aspect.ratio = 1) + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2, color= "black") + scale_color_manual(values=c("#CCCCCC", "#FFCC33")) + scale_fill_manual(values=c("#CCCCCC", "#FFCC33")) 


