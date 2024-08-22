#######################################################
##########  PATHWAY Enrichment analysis
#######################################################
de.pr.tb.st = read.csv("F:/Proteomics/Data/Spleen project_SP/Data/TB vs ST/TB vs ST ALL PROTEINS.csv")
UniProt2Reactome = read.delim("F:/Proteomics/Data/Reactome pathway txt/UniProt2Reactome.txt", header=FALSE)
UniProt2Reactome = UniProt2Reactome[UniProt2Reactome$V6 == "Mus musculus",]
UniProt2Reactome = UniProt2Reactome[UniProt2Reactome$V1 %in% rownames(de.pr.tb.st),]

enrich.path = data.frame(Pathway = unique(UniProt2Reactome$V4),pval_up=1,pval_dn=1)
thr.pval = 0.05
thr.fc = 0

de.up = rownames(de.pr.tb.st)[de.pr.tb.st$P.Value<thr.pval & de.pr.tb.st$logFC>thr.fc]
de.dn = rownames(de.pr.tb.st)[de.pr.tb.st$P.Value<thr.pval & de.pr.tb.st$logFC<thr.fc]

m.up = length(de.up)
n.up = length(rownames(de.pr.tb.st))-m.up

m.dn = length(de.dn)
n.dn = length(rownames(de.pr.tb.st))-m.dn

for(i in 1:length(enrich.path$Pathway)) {
  my.path = enrich.path$Pathway[i]
  my.path.protein = UniProt2Reactome$V1[UniProt2Reactome$V4 == my.path]
  k = length(my.path.protein)
  
  x.up = sum(my.path.protein %in% de.up)
  x.dn = sum(my.path.protein %in% de.dn)
  
  my.path.up.protein = my.path.protein[my.path.protein %in% de.up]
  my.path.dn.protein = my.path.protein[my.path.protein %in% de.dn]
  
  enrich.path$pval_up[i] = phyper(x.up-1, m.up, n.up, k, lower.tail = FALSE)
  enrich.path$pval_dn[i] = phyper(x.dn-1, m.dn, n.dn, k, lower.tail = FALSE)
  
  print(paste(i,k,x.up,x.dn))
}

#enrich.path.sig=enrich.path[c(48,49,51,70,71,72,79,81,104,253,274,240),]
##If you remove the below # only significant will be counted orelse numbered
enrich.path.sig = enrich.path[apply(enrich.path[,2:3], 1, function(x) any(x<thr.pval)),1:3]
rownames(enrich.path.sig) = enrich.path.sig$Pathway
enrich.path.sig = enrich.path.sig[,-1]

##### plotting heatmap
enrich.path.sig.log10 = -log10(enrich.path.sig)

min(enrich.path.sig.log10)
max(enrich.path.sig.log10)

library(pheatmap)

#pdf(file = "/Users/mtalam/OneDrive - UAE University/UAEU/research/shubham/pathway_Enrichment_comp4.pdf",height = 10,width = 6)
pheatmap(enrich.path.sig[,1:2],
         color=colorRampPalette(c("#ff5252","#66b2b2","#303336"))(2000))
dev.off()


#### To export file
#Write the data file name which you want to export
write.table(enrich.path,file="TB vs ST SIG PATHWAYS (q-value).csv",sep=",")


## eg 1 is number of pathway
#for(i in 1:length(enrich.path$Pathway)) {
i=644
my.path = enrich.path$Pathway[i]
  my.path.protein = UniProt2Reactome$V1[UniProt2Reactome$V4 == my.path]
  my.path.protein 
  my.path
#}
  "#B41C1C","#B41C1C","#018141","#018141","#010203"
  