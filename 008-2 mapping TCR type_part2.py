#!/usr/bin/python
#-*-coding:UTF-8-*-
def cal_var_type(number):
    TYPE='None'
    if number==0:
        TYPE='None'
    elif number<=0.0001:
        TYPE='Rare'
    elif number>0.0001 and number<=0.001:
        TYPE='Small'
    elif number>0.001 and number<=0.01:
        TYPE='Medium'
    elif number>0.01 and number<=0.1:
        TYPE='Large'
    else:
        TYPE='Hyperexpanded'
    return TYPE
filename=open('data_to_match.csv').readlines()
filename.pop(0)
#barfile=open('../mouse_cluster_by_lym_detailed.csv').readlines()
barfile=open('../data_reanalysis/mouse_cluster.csv').readlines()
barfile.pop(0)
barcode=[]
for bf in barfile:
    bfline=bf.strip().split(',')
    barcode.append(bfline[0])
    
    
tcr_result=[]
for i in range(12):
    tsub=[]
    for fn in filename:
        fnline=fn.strip().split(',')
        tcr_bar=fnline[0]
        tcr_type=fnline[1]  
        num=int(fnline[2])-1
        if num==i:
            tsub.append([tcr_bar,tcr_type])
        #if re.match(salist[i],tcr_bar):
        #    tsub.append([tcr_bar,tcr_type])
    tcr_result.append(tsub)
    
    
tcr_final={}
for i in range(12):
    sub=tcr_result[i]
    barsub=[]
    typesub=[]
    for j in range(len(sub)):
        barsub.append(sub[j][0])
        typesub.append(sub[j][1])
    typeeach=list(set(typesub))
    type_per={}
    for te in typeeach:
        per=float(typesub.count(te))/len(typesub)
        type_per[te]=per
    for k in range(len(barsub)):
        barline=barsub[k].split('_')
        barfinal=barline[3]+'_'+barline[4]
        perfinal=type_per[typesub[k]]
        typefinal=cal_var_type(perfinal)
        per_type=str(perfinal)+','+typefinal
        tcr_final[barfinal]=per_type

outfile=open('../data_reanalysis/matched_tcr_result.csv','w+')
outfile.write('Barcode,Percentage,Var_type\n')
for bc in barcode:
    if bc in tcr_final:
        outline=bc+','+tcr_final[bc]+'\n'
    else:
        outline=bc+',0,None\n'
    outfile.write(outline)
    
outfile.close()

tcr_specific=open('tcr_specific_cells.csv','w+')
for tf in tcr_final:
    if tf in barcode:
        continue
    else:
        tcr_specific.write(tf+'\n')
tcr_specific.close()
    
    
    
