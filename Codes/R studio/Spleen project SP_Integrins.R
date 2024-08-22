library(pheatmap)
library(ggplot2)
library(VennDiagram)
library(limma)
library(RColorBrewer)
library(UniProt.ws)
library(heatmaply)

thr.fc=2
thr.pv =0.05
Shubham.color = c(brewer.pal(8,"Dark2"))

data.pr.all = read.csv("F:/Proteomics/Data/Spleen project_SP/Excel sheets/16. Spleen project _Integrins.csv")
annotation.all = data.pr.all[,4:5]
data.pr = data.pr.all[,6:23]
data.pr.log2 = log2(data.pr)
rownames(data.pr.log2) = annotation.all$Description
conditions = factor(c("LSK","LSK","LSK","LSK","LSK","LSK","ST","ST","ST","ST","ST","ST","TB","TB","TB","TB","TB","TB"))

#### 2. PCA analysis
data.pr.log2.2 = data.pr.log2[!apply(data.pr.log2,1,function(x) any(is.infinite(x))),]
data.pr.log2.2 = data.pr.log2.2[!apply(data.pr.log2.2,1,function(x) any(is.na(x))),]
pca <- prcomp(t(data.pr.log2.2), center = TRUE, scale = FALSE)
percent.var <- round(100*pca$sdev^2/sum(pca$sdev^2))
pca2 <- cbind(as.data.frame(pca$x),Samples=rownames(pca$x),conditions=conditions)

pca2$conditions

####pdf("PCA plot.pdf",width=6, height=6)
ggplot(pca2, aes(PC1, PC2, color=conditions)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
  ylab(paste0("PC2: ",percent.var[2],"% variance")) +
  coord_fixed()
   dev.off()

#### 3. design
design <- model.matrix(~0+conditions)
colnames(design) <- levels(conditions)
rownames(design) <- colnames(data.pr.log2)

#### 6. heatmap: based on all proteions
df <- as.data.frame(pca2$conditions)
rownames(df) = pca2$Samples
colnames(df) = "conditions"

####pdf(file = "pheatmap protein expression, all proteins.pdf",width=6, height=6);
pheatmap(data.pr.log2.2, cluster_rows=F, show_rownames=T,scale='row',
         cluster_cols=T,main="all protein heatmap",annotation_col = df,color=colorRampPalette(c("#ff5252","#303336","#66b2b2"))(100)) 
dev.off()  


#### 3. design
design <- model.matrix(~0+conditions)
colnames(design) <- levels(conditions)
rownames(design) <- colnames(data.pr.log2.2)

#### 4. lmfit
genotype = unique(conditions)

### YOU MUST SET YOUR REFERENCE CONDITION HERE
reference = "ST"

### MAKING CONTRAST
contrasts = paste0( genotype[genotype !=reference] ,"-", reference)
mc = makeContrasts(contrasts=contrasts, levels=design)
lm.fit <- lmFit(data.pr.log2.2, design)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)


#### 5. differentiually expressed proteins Adelta referene to WT
thr.pv=0.05
diffexp.pr.all = topTable(eb, adjust="BH",coef = 1,number = dim(data.pr.log2.2)[1])
diffexp.pr.all$P.Value[rownames(diffexp.pr.all)==rownames(data.pr.log2.2)[598]]
dim(diffexp.pr.all)
diffexp.pr.sig = diffexp.pr.all[(diffexp.pr.all$adj.P.Val < thr.pv),]
dim(diffexp.pr.sig)
diffexp.pr.sig = diffexp.pr.sig[!is.na(diffexp.pr.sig$adj.P.Val),]
dim(diffexp.pr.sig)

### Volcano plot
plot(diffexp.pr.all$logFC,-log10(diffexp.pr.all$adj.P.Val), pch=16,
     xlab="log2 fold change",ylab="-log10 adj.P.Val", 
     col=ifelse(diffexp.pr.all$adj.P.Val<0.05 & diffexp.pr.all$logFC>log2(2),"#66b2b2",
                ifelse(diffexp.pr.all$adj.P.Val<0.05 & diffexp.pr.all$logFC< -log2(2),"#ff5252","#303336")))
#main = "differentially expressed proteins")
abline(h= -log10(0.05),lty=2)
abline(v=1,lty=2)
abline(v=-1,lty=2)

#650:525

text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="H7BX38"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="H7BX38"]),
     labels="???",col="#004c4c")

text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="P35441"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="P35441"]),
     labels="???",col="#004c4c")

text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="Q62009"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="Q62009"]),
     labels="???",col="#004c4c")
#???
#INTEGRINS
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="Q3V3R4"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="Q3V3R4"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="A0A286YCM1"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="A0A286YCM1"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="P11688"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="P11688"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="A2ARA8"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="A2ARA8"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="B8JK39"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="B8JK39"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="Q3V0T4"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="Q3V0T4"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="Q9QUM0"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="Q9QUM0"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="P24063"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="P24063"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="G5E8F1"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="G5E8F1"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="Q9QXH4"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="Q9QXH4"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="Q542I8"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="Q542I8"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="P09055"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="P09055"]),
     labels="???",col="#004c4c")
text(diffexp.pr.all$logFC[rownames(diffexp.pr.all)=="O54890"],
     -log10(diffexp.pr.all$adj.P.Val[rownames(diffexp.pr.all)=="O54890"]),
     labels="???",col="#004c4c")


#### To export file
####Write the data file name which you want to export
write.table(diffexp.pr.all,file="Integrins for pathway names.csv",sep=",")


#Blue to grey- 
1. "#27408E","#304FAF","#536CB5","#7688BB","#98A5C0","#C0C6CB"
2. "#3D7C91","#303B52","#303336","#53595E","#6A6E76","#838A90"
3. "#27408E","#304FAF","#536CB5","#7688BB","#98A5C0","#C0C6CB","#838A90","#6A6E76","#53595E","#303336"


#234390 #00192b #004c4c

#Red Greens
#"#D91A1A","#B41C1C","#8F1D1E","#6B1F20","#010203","#014122","#018141","#00C060","#00FF7F"
