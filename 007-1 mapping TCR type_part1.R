library(scRepertoire)
library(stringr)
library(ggsci)
theCall <- function(x) {
  if (x == "gene") {
    x <- "CTgene"
  } else if(x == "nt") {
    x <- "CTnt"
  } else if (x == "aa") {
    x <- "CTaa"
  } else if (x == "gene+nt") {
    x <- "CTstrict"
  }
  return(x)
}
clonoTypes<-c(None = 0, Rare = 1e-04, Small = 0.001, Medium = 0.01, 
              Large = 0.1, Hyperexpanded = 1)
cloneCall<-theCall("gene")
df1<-lapply(combined,"[[",cloneCall)
df2<-lapply(combined,"[[","barcode")
data_out<-c()
for(i in 1:length(df1)){
  Barcode<-df2[[i]]
  CTgene<-df1[[i]]
  datasub<-data.frame(Barcode,CTgene)
  datasub$Num<-i
  data_out<-rbind(data_out,datasub)
}

write.csv(data_out,"VDJ_analysis/data_to_match.csv",quote=F,row.names=F)
