#!/usr/bin/python
#-*-coding:UTF-8-*-
import os,re
filelist=os.listdir("/data2/project/sephin_mouse_data/TCR/")
def dataprocess(sample,data):
    datadir="/data2/project/sephin_mouse_data/TCR/"
    os.chdir(datadir)
    samlist1=[]
    samlist2=[]
    for da in data:
        if re.match(sample,da) and re.search('.fq.gz',da):
            if re.search('_1',da):
                samlist1.append(da)
            elif re.search('_2',da):
                samlist2.append(da)
    samlist1=sorted(samlist1)
    samlist2=sorted(samlist2)
    samline1=' '.join(samlist1)
    samline2=' '.join(samlist2)
    os.system('cat '+samline1+' > processed/'+sample+'_S1_L001_R1_001.fastq.gz')
    os.system('cat '+samline2+' > processed/'+sample+'_S1_L001_R2_001.fastq.gz')
    os.chdir('/home/wangrj/wangrj/wangrj/MOUSE_SEPHIN1_RESULTS')
    os.system('/home/wangrj/wangrj/wangrj/softwares/cellranger-6.0.0/cellranger vdj --localcores=12 --id='+sample+'_output --reference=/data1/database/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 --fastqs=/data2/project/sephin_mouse_data/TCR/processed/ --sample='+sample)

#main
samplelist=['BN1d0V5-T','BN2d0V5-T','BS1d0V5-T','BS2d0V5-T','TN1d15V5-T','TN2d15V5-T','TS1d15V5-T','TS2d15V5-T','BN1d15V5-T','BN2d15V5-T','BS1d15V5-T','BS2d15V5-T']
for sa in samplelist:
    dataprocess(sa,filelist)
