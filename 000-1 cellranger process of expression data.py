#!/usr/bin/python
#-*-coding:UTF-8-*-
import os,re
filelist=os.listdir("/data2/project/sephin_mouse_data/")
def dataprocess(sample,data):
    datadir="/data2/project/sephin_mouse_data/"
    os.chdir(datadir)
    samlist1=[]
    samlist2=[]
    for da in data:
        if re.match(sample,da) and re.search('.fq.gz',da):
            if re.search('_1',da):
                samlist1.append(da)
            elif re.search('_2',da):
                samlist2.append(da)
    samline1=' '.join(samlist1)
    samline2=' '.join(samlist2)
    os.system('cat '+samline1+' > processed_fastqs/'+sample+'_S1_L001_R1_001.fastq.gz')
    os.system('cat '+samline2+' > processed_fastqs/'+sample+'_S1_L001_R2_001.fastq.gz')
    os.chdir('/home/wangrj/wangrj/wangrj/MOUSE_SEPHIN1_RESULTS')
    os.system('/home/wangrj/wangrj/wangrj/softwares/cellranger-6.0.0/cellranger count --localcores=12 --id='+sample+'_output --transcriptome=/data1/database/refdata-cellranger-mm10-1.2.0 --fastqs=/data2/project/sephin_mouse_data/processed_fastqs/ --sample='+sample+' --expect-cells=8000') 

#main
samplelist=['BN1d0V5','BN2d0V5','BS1d0V5','BS2d0V5','TN1d15V5','TN2d15V5','TS1d15V5','TS2d15V5','BN1d15V5','BN2d15V5','BS1d15V5','BS2d15V5']
for sa in samplelist:
    dataprocess(sa,filelist)
