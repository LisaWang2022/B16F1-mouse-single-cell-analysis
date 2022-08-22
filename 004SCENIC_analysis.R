library(SCENIC)
library(Seurat)
library(GENIE3)
library(AUCell)
library(foreach)
library(RcisTarget)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)

counts<-data@assays$RNA@counts
counts<-as.matrix(counts)
Gene<-rownames(counts)
countsnew<-cbind(Gene,counts)
write.table(counts,"lym_counts.txt",quote=F,sep="\t",row.names=F)
lym<-data

data<-readRDS("../../mouse_all_remove_RBC.rds")
counts<-data@assays$RNA@counts
counts<-as.matrix(counts)
Gene<-rownames(counts)
countsnew<-cbind(Gene,counts)
write.table(countsnew,"all_counts.txt",quote=F,sep="\t",row.names=F)
cell_type<-data$Cluster_name
Cell<-rownames(data@meta.data)
meta<-data.frame(Cell,cell_type)
write.table(meta,"all_meta.txt",quote=F,sep="\t",row.names=F)

#SCENIC
Sample<-data$SAMPLE
Cluster<-data$seurat_clusters
Celltype<-data$Cluster_name
Group<-data$GROUP
cellinfo<-data.frame(Sample,Cluster,Celltype,Group)
scenicOptions <- initializeScenic(org="mgi", dbDir="./SCENIC_database", nCores=10)
### 鍏辫〃杈剧綉缁?
#genesKept <- geneFiltering(counts, scenicOptions)
#exprMat_filtered <- counts[genesKept, ]
#runCorrelation(exprMat_filtered, scenicOptions)
#exprMat_filtered_log <- log2(exprMat_filtered+1)
#runGenie3(exprMat_filtered_log, scenicOptions)#time consuming about 7 days.



#randomly select 5000 cells for GENIE3
names<-colnames(count)
number<-length(names)
snumber<-sample(number,5000)
exp_sub<-counts[,snumber]
genesKept <- geneFiltering(exp_sub, scenicOptions)
exp_sub_filtered <- exp_sub[genesKept, ]
runCorrelation(exp_sub_filtered, scenicOptions)
exp_sub_log <- log2(exp_sub_filtered+1)
runGenie3(exp_sub_log, scenicOptions)

### Build and score the GRN 
exprMat_log <- log2(counts+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
#optional
#runSCENIC_4_aucell_binarize(scenicOptions)

#results
coexp<-readRDS("int/1.6_tfModules_asDF.Rds")


#identify cell-specific regulon
#cellinfo<-data@meta.data
saveRDS(cellinfo,"cellinfo.rds")
cellinfo<-readRDS("cellinfo.rds")
regulonAUC<-loadInt(scenicOptions,"aucell_regulonAUC")
rss<-calcRSS(AUC=getAUC(regulonAUC),cellAnnotation=cellinfo[colnames(regulonAUC),"Celltype"])
rssPlot<-plotRSS(rss)
rssPlot$plot

rss1<-calcRSS(AUC=getAUC(regulonAUC),cellAnnotation=cellinfo[colnames(regulonAUC),"Group"])
rssPlot1<-plotRSS(rss1)
rssPlot1$plot

rss2<-calcRSS(AUC=getAUC(regulonAUC),cellAnnotation=cellinfo[colnames(regulonAUC),"SAMPLE"])
rssPlot2<-plotRSS(rss2)
rssPlot2$plot

pdf("rssPlot_by_sample.pdf",height=5,width=8)
rssPlot2$plot
dev.off()

#heatmap
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellinfo), cellinfo$Celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pdf("heatmap_by_cluster.pdf",height=7,width=9)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()

regulonActivity_byCellType1 <- sapply(split(rownames(cellinfo), cellinfo$Sample),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled1 <- t(scale(t(regulonActivity_byCellType1), center = T, scale=T))

regulonActivity_byCellType_Scaled1<-regulonActivity_byCellType_Scaled1[order(rownames(regulonActivity_byCellType_Scaled1)),]
#regulonActivity_byCellType_Scaled1<-data.frame(regulonActivity_byCellType_Scaled1)

colgroup<-c("Day0_Blood","Day0_Blood","Day15_Blood","Day15_Blood","Day15_Tumor","Day15_Tumor")
col1<-as.matrix(colgroup)
col2<-t(col1)
colnames(col2)<-colnames(regulonActivity_byCellType_Scaled1)
col2<-data.frame(t(col2))
colnames(col2)<-"Tissue"#
cname<-c("Normal","Sephin1","Normal","Sephin1","Normal","Sephin1")
pdf("heatmap_of_regulons_20210923.pdf",height=6,width=7)
pheatmap(regulonActivity_byCellType_Scaled1,cluster_rows=F,cluster_cols=F,
         gaps_col=c(2,4),annotation_col = col2,labels_col = cname)
dev.off()

#regulons写入文件
regnew<-c()
for(i in 1:length(regulons)){
  rname<-names(regulons)[i]
  rvec<-paste(regulons[[i]],collapse=",")
  deach<-data.frame(rname,rvec)
  regnew<-rbind(regnew,deach)
}
colnames(regnew)<-c("Name","Genes")
write.csv(regnew,"regulons.csv",quote=F)