#!/usr/bin/env python
#coding: utf8
import column
import sys
import collections
import re

##
# Takes a FASTA file. Places the posistions where the columns there are more than 50% indels in a list.
# Call column. Caled by TakeAway
##

def Indels(fil):
    BadColum = []
    columns=column.column(fil)
    for k in columns:
        Indel = re.findall('-',columns[k])
        LenIndel = len(Indel)
        LenColumn = len(columns[k])
        if float(LenIndel)/LenColumn>0.5:
            Pos=int(k)
            BadColum.append(Pos)
        else:
            pass
    return BadColum

##
# Takes a FASTA file. Places the posistions where the columns have at least 50% of amino acids
# are unique in a list.
# Call column, Caled by TakeAway.
##

def Unique(fil):
    BadColum = []
    columns=column.column(fil)
    for k in columns:
        results=collections.Counter(columns[k])
        cunter=0
        for key in results:
            if results[key]==1:
                cunter=cunter+1
            else:
                pass
        if float(cunter)/len(columns[k])>=0.5:
            Pos=int(k)
            BadColum.append(Pos) 
        else:
            pass
    return BadColum

##
# Takes a FASTA file. Places the posistions where the columns that have no amino acid appers more
# than twice in a list.
# Call column. Caled by TakeAway.
##
def MoreThan2(fil):
    BadColum = [] 
    columns=column.column(fil)
    for k in columns:
        results=collections.Counter(columns[k])
        number=[]
        for key in results:
            number.append(results[key]) 
            greater2=[]
            for i in range(len(number)):
                if int(number[i])>2:
                    greater2.append(number[i])
                else:
                    pass
         
        if len(greater2)==0:
            Pos=int(k)
            BadColum.append(Pos)
        else:
            pass
    return BadColum
##
# Takes a FASTA file. Takes away the bad columns in the sequnce.
# Call Indels, Unique, MoreThan2 and SeqDic. 
##
def TakeAway(fil):
    seq_dic=column.SeqDic(fil)
    BadColum1 = Indels(fil)
    BadColum2 = Unique(fil)
    BadColum3 = MoreThan2(fil)
    BadColum = BadColum1+BadColum2+BadColum3
    Colums = list(set(BadColum))
    for k in range(len(Colums)):
        Pos=int(Colums[k])
        for key in seq_dic:
            seq_dic[key]=seq_dic[key][:Pos]+seq_dic[key][Pos+1:]
    if len(seq_dic)!=0:
        pass
    else:
        sys.stdout.write('ERROR: All the columns has been removed') 
    return seq_dic  
    
     
     

            
        
