library(tibble)
library(ggplot2)
library(cowplot)
library(Seurat)
library(fgsea)
library(msigdbr)
library(ClusterProfiler)
#differentiate genes were calculated by FindMarkers
filelist<-dir("diff_gene_by_group/tumor_diff/")
###################
#use fgsea
###################

#get gmt file
gobp_pathway<-gmtPathways("/home/wangrj/wangrj/wangrj/mouse_breasttumor/BD_second/gskb_gmt_files/unfiltered/MousePath_GO_BP_gmt.gmt")
kegg_pathway<-gmtPathways("/home/wangrj/wangrj/wangrj/mouse_breasttumor/BD_second/gskb_gmt_files/unfiltered/MousePath_Metabolic_KEGG_gmt.gmt")
#draw picture of each cell type in tumor sample
draw_gsva_barplot<-function(fl){
  #get ranks from foldchange data
  diffgene<-read.csv(paste0("diff_gene_by_group/tumor_diff/",fl),row.name=1)
  diffgene<-subset(diffgene,p_val_adj<0.05)
  diffgene<-subset(diffgene,p_val<0.05)
  diffgene$gene<-rownames(diffgene)
  diffgene$gene<-toupper(diffgene$gene)
  diffgene<-diffgene[!duplicated(diffgene$gene),]
  rownames(diffgene)<-diffgene$gene
  rownames(diffgene)<-toupper(rownames(diffgene))
  diff_frame<-data.frame(rownames(diffgene),diffgene$avg_logFC)
  ranks<-deframe(diff_frame)

  #gsea analysis
  fgseaRes <- fgsea(pathways = gobp_pathway, 
                  stats    = ranks,
                  nperm=1000,
                  minSize  = 3,
                  maxSize  = 500)

  fgdata<-fgseaRes[order(fgseaRes$NES),]
  upgene<-fgdata[(length(rownames(fgdata))-20):length(rownames(fgdata)),]
  downgene<-fgdata[1:20,]
  geneall<-rbind(downgene,upgene)
  #geneall<-subset(fgdata,pval<=0.05)
  geneall$Group<-"NA"
  for(i in 1:length(geneall$ES)){
    if(geneall$NES[i]>0 & geneall$pval[i]<=0.05){
      geneall$Group[i]="Up"
    }
    else if(geneall$NES[i]<0 & geneall$pval[i]<=0.05){
      geneall$Group[i]="Down"
    }
    else if(geneall$pval[i]>0.05){
      geneall$Group[i]="Unsig"
    }
  }
  colnames(geneall)[1]<-"Pathway"

  geneall$Group<-factor(geneall$Group,levels=c("Up","Down","Unsig"))
  geneall$Pathway<-gsub("GO_BP_MM_","",geneall$Pathway)
  geneall$Pathway<-factor(geneall$Pathway,levels=as.character(geneall$Pathway))
  #barplot
  clname<-gsub("_diff_genes.csv","",fl)
  clname<-gsub(" ","_",clname)
  pdfname<-paste0("gsva_plots/tumor_gsva_barplot_gobp_",clname,".pdf")
  titlename=paste0("GSEA_analysis_",clname)
  ggplot(geneall,aes(x=Pathway,y=NES,fill=Group))+
    geom_bar(stat="identity")+
    scale_fill_manual(values = c("#FF0099","#6699FF","lightgrey"))+
    #scale_fill_manual(values = c("lightgrey"))+
    coord_flip()+guides(fill = "none")+theme_classic()+
    labs(title=titlename)+
    geom_text(data = subset(geneall, NES>0 & pval<=0.05),aes(x=Pathway, y= 0, label= paste0(Pathway,"  ")),hjust = "outward")+
    geom_text(data = subset(geneall, NES>0 & pval>0.05),aes(x=Pathway, y= 0, label= paste0(Pathway,"  ")),color="darkgrey",hjust = "outward")+
    geom_text(data = subset(geneall, NES<0 & pval>0.05),aes(x=Pathway, y= 0, label= paste0("  ",Pathway)),color="darkgrey",hjust = "inward")+
    geom_text(data = subset(geneall, NES<0 & pval<=0.05),aes(x=Pathway, y= 0, label= paste0("  ",Pathway)),hjust = "inward")+
    theme(plot.title=element_text(hjust=0.5,size=18),axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title=element_text(size=15))
    
  ggsave(pdfname,height=4,width=9)
}

#TCR analysis
tcr<-subset(mouse.data,TCR_type!="None")
Idents(tcr)<-tcr$TCR_type
all_markers<-FindAllMarkers(tcr,min.pct=0.25)
diffgene<-subset(all_markers,p_val_adj<0.05 & cluster=="Hyperexpanded")
rownames(diffgene)<-toupper(rownames(diffgene))
diff_frame<-data.frame(rownames(diffgene),diffgene$avg_logFC)
ranks<-deframe(diff_frame)
fgseaRes <- fgsea(pathways = gobp_pathway, 
                  stats    = ranks,
                  nperm=1000,
                  minSize  = 3,
                  maxSize  = 500)

fgdata<-fgseaRes[order(fgseaRes$NES,decreasing=T),]
upgene<-fgdata[1:20,]
geneall<-upgene[order(upgene$NES),]
colnames(geneall)[1]<-"Pathway"
geneall$Pathway<-gsub("GO_BP_MM_","",geneall$Pathway)
geneall$Pathway<-factor(geneall$Pathway,levels=as.character(geneall$Pathway))
geneall$Significance<-log2(geneall$pval)*(-1)
ggplot(geneall,aes(x=Pathway,y=NES,fill=Significance))+
  geom_bar(stat="identity")+
  scale_fill_gradient(low="#FF99CC",high="#FF0099")+
  coord_flip()+theme_classic()+labs(title="GSEA analysis of hyperexpanded cluster")+
  theme(plot.title=element_text(hjust=1,size=18),axis.title=element_text(size=15))
ggsave("TCR_gsea_hyperexpanded.pdf",height=4,width=9)
