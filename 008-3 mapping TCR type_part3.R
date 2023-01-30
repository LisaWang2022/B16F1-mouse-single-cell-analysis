library(scRepertoire)
library(stringr)
library(ggsci)
Clonotype<-read.csv("matched_tcr_result.csv",header=T)
mouse.data$TCR_type<-Clonotype$Var_type
pdf("cluster_by_tcr_type_graph.pdf",height=5,width=7)
DimPlot(mouse.data,group.by="TCR_type",cols=c("red","orange","blue","grey","green"))
dev.off()
pdf("cluster_by_tcr_type_graph_splitted.pdf",height=10,width=9)
DimPlot(mouse.data,split.by="SAMPLE",group.by="TCR_type",cols=c("red","orange","blue","grey","green"),ncol=2)
dev.off()


Idents(mouse.data)<-mouse.data$TCR_type
tcr_all_markers<-FindAllMarkers(mouse.data,min.pct=0.25)
write.csv(tcr_all_markers,"tcr_all_markers.csv",quote=F)

