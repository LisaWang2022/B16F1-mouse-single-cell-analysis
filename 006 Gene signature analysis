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
