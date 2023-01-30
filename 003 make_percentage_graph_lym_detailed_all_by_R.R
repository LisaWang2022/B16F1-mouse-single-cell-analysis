#make lim percentage
library(ggplot2)
library(Seurat)
mouse.data<-readRDS()
cdetail<-mouse.data$Cluster_detailed
cdetailnew<-gsub("Cyt_NK","NK",cdetail)
cdetailnew<-gsub("Non_cyt_NK","NK",cdetailnew)
mouse.data$Cluster_detailed_new<-cdetailnew
saveRDS(mouse.data,"all_data_mouse_20210911.rds")

Cluster<-mouse.data$Cluster_detailed_new
Sample<-mouse.data$SAMPLE
Percentage<-c()
clu<-unique(Cluster)
sam<-unique(Sample)
for(sa in sam){
  sub1=subset(mouse.data,SAMPLE==sa)
  for(cl in clu){
    if(cl %in% sub1$Cluster_detailed_new){
      sub2=subset(sub1,Cluster_detailed_new==cl)
      per<-length(sub2$orig.ident)/length(sub1$orig.ident)
    }
    else{
      per=0
    }
    dsub<-data.frame(sa,cl,per)
    Percentage<-rbind(Percentage,dsub)
  }
}

colnames(Percentage)<-c("Sample","Cluster","Percentage")
Percentage$Group<-"NA"
Percentage$Tissue<-"NA"
for(i in 1:length(Percentage$Sample)){
  salist<-str_split(Percentage$Sample[i],pattern="_",simplify=T)
  sa1<-paste(salist[1],salist[2],sep="_")
  Percentage$Tissue[i]<-sa1
  Percentage$Group[i]<-salist[3]
}

msub<-subset(mouse.data,subset=Cluster_name=="NK cells" | Cluster_name=="NKT cells" | Cluster_name=="Cd4+ T cells" | Cluster_name=="Cd8+ T cells")
csub<-unique(msub$Cluster_detailed_new)
sub<-c()
for(i in 1:length(Percentage$Sample)){
  if(Percentage$Cluster[i] %in% csub){
    sub<-rbind(sub,Percentage[i,])
  }
}
sub$Cluster<-factor(sub$Cluster,levels=c("NK","NKT","Exhausted_Cd8","Cyt_Cd8","Naïve_Cd8","Effector_Cd4","Treg_Cd4","Naïve_Cd4"))
p1<-ggplot(sub,aes(x=Tissue,y=Percentage,fill=Group,group=Group))+
  geom_bar(stat="identity",position="dodge",color="black")+
  scale_fill_manual(values = c("#FF0099","#6699FF"))+theme_classic()+
  theme(text=element_text(size=15),axis.text.x=element_text(angle=90,size=12))+
  facet_wrap(~Cluster,scales="free_y",ncol=4)
celltypes<-unique(sub$Cluster)
tidays<-unique(mouse.data$Tissue_days)
pall<-list()
starlist_chisq<-list()
tilist<-unique(sub$Tissue)
for(ct in celltypes){
  plist<-c()
  slist<-c()
  for(ti in tidays){
    sub1<-subset(mouse.data,Tissue_days==ti)
    result<-make_fisher_chisq_test(sub1,ct,c("Normal","Sephin1"),"GROUP","Cluster_detailed_new")
    peach<-result[1]
    pvalue<-as.numeric(peach)
    stars<-result[2]
    plist<-c(plist,peach)
    slist<-c(slist,stars)
  }
  pall<-c(pall,list(plist))
  starlist_chisq<-c(starlist_chisq,list(slist))
}
names(starlist_chisq)<-celltypes
df<-list()
x1all<-list()
x2all<-list()
x3all<-list()
y1all<-list()
y2all<-list()
#ymax<-max(sub$Percentage)
for(cl in celltypes){
  sub1<-subset(sub,Cluster==cl)
  df<-c(df,list(sub1))
  x1list<-c()
  x2list<-c()
  x3list<-c()
  y1list<-c()
  y2list<-c()
  ymax<-max(sub1$Percentage)
  yy<-ymax/10
  per<-sub1$Percentage
  for(i in 1:3){
    x1<-i-0.3
    x2<-i+0.3
    x3<-i
    sub2<-subset(sub1,Tissue==tilist[i])
    y1<-max(sub2$Percentage)+yy
    y2<-y1+yy+yy+yy
    x1list<-c(x1list,x1)
    x2list<-c(x2list,x2)
    x3list<-c(x3list,x3)
    y1list<-c(y1list,y1)
    y2list<-c(y2list,y2)
  }
  x1all<-c(x1all,list(x1list))
  x2all<-c(x2all,list(x2list))
  x3all<-c(x3all,list(x3list))
  y1all<-c(y1all,list(y1list))
  y2all<-c(y2all,list(y2list))
}

geom_star<-function(df,xstartlist,xendlist,xlist,y1list,y2list,starlist){
  list(geom_text(data=df,aes(x=xlist[1],y=y2list[1],label=paste("",starlist[1],sep="\n"))),
       geom_segment(data=df,aes(x =xstartlist[1], y =y1list[1], xend = xendlist[1], yend = y1list[1])),
       geom_text(data=df,aes(x=xlist[2],y=y2list[2],label=paste("",starlist[2],sep="\n"))),
       geom_segment(data=df,aes(x =xstartlist[2], y =y1list[2], xend = xendlist[2], yend = y1list[2])),
       geom_text(data=df,aes(x=xlist[3],y=y2list[3],label=paste("",starlist[3],sep="\n"))),
       geom_segment(data=df,aes(x =xstartlist[3], y =y1list[3], xend = xendlist[3], yend = y1list[3])))
}
p2<-p1+geom_star(df[[1]],x1all[[1]],x2all[[1]],x3all[[1]],y1all[[1]],y2all[[1]],starlist_chisq[[1]])+
  geom_star(df[[2]],x1all[[2]],x2all[[2]],x3all[[2]],y1all[[2]],y2all[[2]],starlist_chisq[[2]])+
  geom_star(df[[3]],x1all[[3]],x2all[[3]],x3all[[3]],y1all[[3]],y2all[[3]],starlist_chisq[[3]])+
  geom_star(df[[4]],x1all[[4]],x2all[[4]],x3all[[4]],y1all[[4]],y2all[[4]],starlist_chisq[[4]])+
  geom_star(df[[5]],x1all[[5]],x2all[[5]],x3all[[5]],y1all[[5]],y2all[[5]],starlist_chisq[[5]])+
  geom_star(df[[6]],x1all[[6]],x2all[[6]],x3all[[6]],y1all[[6]],y2all[[6]],starlist_chisq[[6]])+
  geom_star(df[[7]],x1all[[7]],x2all[[7]],x3all[[7]],y1all[[7]],y2all[[7]],starlist_chisq[[7]])+
  geom_star(df[[8]],x1all[[8]],x2all[[8]],x3all[[8]],y1all[[8]],y2all[[8]],starlist_chisq[[8]])

ggsave("mouse_sephin1_cluster_percentage_lym_detailed_with_pvalue.pdf",height=4,width=11)

name="mouse_sephin1_cluster_percentage_lym_detailed.pdf"
p1
ggsave(name,height=4,width=11)






