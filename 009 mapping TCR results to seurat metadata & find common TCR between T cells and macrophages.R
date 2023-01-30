library(stringr)
library(hash)
combined<-readRDS("/storage/work/wangrj/MOUSE_SEPHIN1_RESULTS/VDJ_analysis/scREpertoire_combined.rds")
exp_barcode<-colnames(mouse.data)
#merge combined to one single dataset
tcr_merged<-c()
for(i in 1:12){
  tdata<-combined[[i]]
  bar1<-tdata$barcode
  barlist<-str_split(bar1,"_",simplify=T)
  bar2<-barlist[,4]
  barcode<-paste(bar2,i,sep="_")
  tsub<-tdata[,4:14]
  tnew<-cbind(barcode,tsub)
  tcr_merged<-rbind(tcr_merged,tnew)
}
tcr_barcode<-tcr_merged$barcode
tcr_hash<-hash(keys=as.character(tcr_barcode),values=data.frame(t(as.matrix(tcr_merged))))
tcr_meta<-c()
for(eb in exp_barcode){
  if(eb %in% tcr_barcode){
    tcr_meta<-rbind(tcr_meta,as.character(tcr_hash[[as.character(eb)]]))
  }
  else{
    nalist<-rep("NA",12)
    tcr_meta<-rbind(tcr_meta,t(nalist))
  }
}
colnames(tcr_meta)<-colnames(tcr_merged)
tcr_meta_sub<-data.frame(tcr_meta[,2:12])
mouse.data@meta.data<-cbind(mouse.data@meta.data,tcr_meta_sub)
saveRDS(mouse.data,"mouse_all_Sephin1_refiltered.rds")
#find common TCR sequence between macrophage and Tcells
#compare each CTgene between Tcell and macrophages. If common, annotated as "Common_in_Cd4" or "Common_in_Cd8" 
msub<-subset(mouse.data,CTgene!="NA")
ctgenelist<-unique(msub$CTgene)
Cluster_name<-mouse.data$Cluster_name
Cluster_detailed<-mouse.data$Cluster_detailed
CTgene<-mouse.data$CTgene
TCR_type<-mouse.data$TCR_type
outdata<-data.frame(Cluster_name,Cluster_detailed,CTgene,TCR_type)
write.csv(outdata,"TCR_data.csv")

macdata<-subset(outdata,Cluster_name=="Macrophages" & CTgene!="NA")
cd4data<-subset(outdata,Cluster_name=="Cd4+ T cells" & CTgene!="NA")
cd8data<-subset(outdata,Cluster_name=="Cd8+ T cells" & CTgene!="NA")

ctgene_mac<-unique(macdata$CTgene)
ctgene_cd4<-unique(cd4data$CTgene)
ctgene_cd8<-unique(cd8data$CTgene)
Anno_mac<-c()
for(i in 1:length(outdata$Cluster_name)){
  if(Cluster_name[i]=="Macrophages"){
    if(CTgene[i] %in% ctgene_cd4){
      Anno_mac<-c(Anno_mac,"Common in Cd4")
    }
    else if(CTgene[i] %in% ctgene_cd8){
      Anno_mac<-c(Anno_mac,"Common in Cd8")
    }
    else{
      if(outdata$TCR_type[i]=="None"){
        Anno_mac<-c(Anno_mac,"Conventional mac")
      }
      else{
        Anno_mac<-c(Anno_mac,"Unique in mac")
      }
    }
  }
  else{
    Anno_mac<-c(Anno_mac,"Non mac")
  }
}
mouse.data$Anno_mac<-Anno_mac
mouse.data$Anno_mac<-factor(mouse.data$Anno_mac,levels=c("Common in Cd4","Common in Cd8","Unique in mac","Conventional mac","Non mac"))
DimPlot(mouse.data,group.by="Anno_mac",cols=c("red","blue","yellow","grey70","grey90"))

pdf("dimplot_of_common_tcr_mac.pdf",height=6,width=12)
DimPlot(mouse.data,group.by="Anno_mac",cols=c("blue","red","#996600","grey70","grey90"),
        split.by="SAMPLE",ncol=4,pt.size = 2)+labs(title="")
dev.off()

tumorsub<-subset(mouse.data,TISSUE=="Tumor")
pdf("dimplot_of_common_tcr_mac_tumor_only.pdf",height=5,width=9)
DimPlot(tumorsub,group.by="Anno_mac",cols=c("blue","red","#996600","grey70","grey90"),
        split.by="GROUP",ncol=4,pt.size = 0.2)+labs(title="")
dev.off()

Anno_mac<-gsub("Common in Cd4","Cd4-share",Anno_mac)
Anno_mac<-gsub("Common in Cd8","Cd8-share",Anno_mac)
Anno_mac<-gsub("Unique in mac","Unique",Anno_mac)
Anno_mac<-gsub("Conventional mac","Conventional-mac",Anno_mac)
Anno_mac<-gsub("Non mac","Non-mac",Anno_mac)
mouse.data$Anno_mac<-Anno_mac
mouse.data$Anno_mac<-factor(mouse.data$Anno_mac,levels=c("Cd4-share","Cd8-share","Unique","Conventional-mac","Non-mac"))
pdf("dimplot_of_common_tcr_mac.pdf",height=6,width=12)
DimPlot(mouse.data,group.by="Anno_mac",cols=c("blue","red","#996600","grey70","grey90"),
        split.by="SAMPLE",ncol=4,pt.size = 2)+labs(title="")
dev.off()

tumorsub<-subset(mouse.data,TISSUE=="Tumor")
pdf("dimplot_of_common_tcr_mac_tumor_only.pdf",height=5,width=9)
DimPlot(tumorsub,group.by="Anno_mac",cols=c("blue","red","#996600","grey70","grey90"),
        split.by="GROUP",ncol=4,pt.size = 0.2)+labs(title="")
dev.off()


