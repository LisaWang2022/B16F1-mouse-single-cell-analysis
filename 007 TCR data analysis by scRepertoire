library(scRepertoire)
library(stringr)
library(ggsci)
csvlist<-c("BN1d0.csv","BN2d0.csv","BS1d0.csv","BS2d0.csv","BN1d15.csv","BN2d15.csv","BS1d15.csv","BS2d15.csv","TN1d15.csv","TN2d15.csv","TS1d15.csv","TS2d15.csv")
contig_list<-list()
for(cl in csvlist){
  filename<-paste0("/storage/work/wangrj/MOUSE_SEPHIN1_RESULTS/VDJ_analysis/immunarch_files/TCR/",cl)
  csvfile<-read.csv(filename,stringsAsFactors=F)
  csvfilter<-csvfile[,c('barcode','is_cell','contig_id','high_confidence','length','chain',
                        'v_gene','d_gene','j_gene','c_gene','full_length','productive','cdr3',
                        'cdr3_nt','reads','umis','raw_clonotype_id','raw_consensus_id')]
  contig_list<-c(contig_list,list(csvfilter))
}

combined <- combineTCR(contig_list,samples=c("Day0_Blood","Day0_Blood","Day0_Blood","Day0_Blood","Day15_Blood","Day15_Blood","Day15_Blood","Day15_Blood","Day15_Tumor","Day15_Tumor","Day15_Tumor","Day15_Tumor"),ID=c("Normal","Normal","Sephin1","Sephin1","Normal","Normal","Sephin1","Sephin1","Normal","Normal","Sephin1","Sephin1"),cells="T-AB")
                      
combined<-addVariable(combined,name="TISSUE",variables=c("Blood","Blood","Blood","Blood","Blood","Blood","Blood","Blood","Tumor","Tumor","Tumor","Tumor"))
combined<-addVariable(combined,name="GROUP",variables=c("Normal","Normal","Sephin1","Sephin1","Normal","Normal","Sephin1","Sephin1","Normal","Normal","Sephin1","Sephin1"))
combined<-addVariable(combined,name="DAYS",variables=c("Day0","Day0","Day0","Day0","Day15","Day15","Day15","Day15","Day15","Day15","Day15","Day15"))
pdf("conbined_quantContig.pdf",height=4,width=6)
quantContig(combined, cloneCall="gene+nt", scale = T)+scale_fill_npg()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
pdf("conbined_abundanceContig.pdf",height=4,width=6)
abundanceContig(combined, cloneCall = "gene+nt", scale = F,group="GROUP")+
  scale_fill_manual(values = c("yellow","green","blue","red","purple","orange"))

dev.off()
pdf("combined_clonalproportion.pdf",height=4,width=5)
clonalProportion(combined, cloneCall = "gene")+scale_fill_npg()+theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
