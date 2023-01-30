#Fisher exact test
cname<-mouse.data$Cluster_name
group<-mouse.data$GROUP
mdata<-data.frame(cname,group)
tmdata<-table(mdata)
temp<-data.frame(tmdata)
numall<-length(cname)
num1<-temp$Freq[3]
num2<-temp$Freq[19]
#dsub<-rbind(c(num1,(22469-num1)),c(num2,(46062-num2)))
dsub<-rbind(c(num1,num2),c((22469-num1),(46062-num2)))
fisher.test(dsub)
#卡方检验，chisq检验细胞占比显著性
#matrix making and chisq pvalue return
mdata<-readRDS("/home/wangrj/wangrj/wangrj/MOUSE_SEPHIN1_RESULTS/data_reanalysis/mouse_all_data_20211017.rds")
make_fisher_chisq_test<-function(data,celltype,samplelist,samplename,clustername){
  Cluster<-data[[clustername]]
  Sample<-data[[samplename]]
  temp_data<-data.frame(Cluster,Sample)
  v_data<-length(rownames(subset(temp_data,Sample==samplelist[1] & Cluster==celltype)))
  adj_data<-length(rownames(subset(temp_data,Sample==samplelist[2] & Cluster==celltype)))
  v_nodata<-length(rownames(subset(temp_data,Sample==samplelist[1] & Cluster!=celltype)))
  adj_nodata<-length(rownames(subset(temp_data,Sample==samplelist[2] & Cluster!=celltype)))
  mat<-rbind(c(v_data,v_nodata),c(adj_data,adj_nodata))
  numN<-sum(mat)
  numT1<-(v_data+v_nodata)*min(v_data,v_nodata)/numN
  numT2<-(adj_data+adj_nodata)*min(adj_data,adj_nodata)/numN
  numT<-min(numT1,numT2)
  if(numT>=5 & numN>=40){
    result<-chisq.test(mat)
  }
  else{
    result<-fisher.test(mat)
  }
  pvalue<-result$p.value
  pvalue<-as.numeric(pvalue)
  if(pvalue>=0.05){
    stars<-"NS"
  }
  else if(pvalue<0.05 & pvalue>=0.01){
    stars<-"*"
  }
  else if(pvalue<0.01 & pvalue>=0.001){
    stars<-"**"
  }
  else if(pvalue<0.001 & pvalue>=0.0001){
    stars<-"***"
  }
  else{
    stars<-"****"
  }
  resultlist<-c(pvalue,stars,celltype,samplelist[2],celltype)
  return(resultlist)
}
