#!/usr/bin/env python
#coding: utf8
import sys
from Bio import AlignIO

##
# Takes a file. Takes tke name of the sequnce and place it as a key in the directory. The value in the 
# directory is the sequence that responds to the key name. 
##
def SeqDic(lines_list):
    dic = {}	
    key = ''
    i=0    
    for line in lines_list:
    	if line[0]=='>' and len(line.split())>1:
            list_line=line.split()
            dic[list_line[0][1:]]=''
	    key = list_line[0][1:]
            i=i+1	
	elif line[0]=='>':
            dic[line[0][1:-1]]=''
	    key = line[0][1:-1]
            i=i+1		 
	else:
	    line=line.replace('\n','')
            dic[key]=dic[key]+line         	         
    return dic

##
# Takes a FASTA file. Places the order thata the spices is placed in the files as the key in a dictonary and the species
# names as the values.
##
def order(lines_list):
    dic = {}
    i = 0
    for line in lines_list:
    	if line[0]=='>' and len(line.split())>1:
            list_line=line.split()
            dic[str(i)]=list_line[0][1:]
            i=i+1	
	elif line[0]=='>':
            dic[str(i)]=line[0][1:-1]
            i=i+1		 
	else:
	    pass       	         
    return dic
        

##
# Takes a Dictonary and checks if the values are in the same lenght if not gives an error.
##
def SameLength(dic):
    key1 = dic.keys()[0] # takes out the first key in the dictonary
    for k in dic:
        if len(dic[key1])==len(dic[k]):
	    pass
	else:
	    sys.stdout.write('ERROR: The sequnces have not the same lenght')
	    break

##
# Take a file and takes the i:th argument in every value and puts it in a new directonary where the key is the a number from 0 to the # length of the sequence.
## 
def column(lines_list):
    dic=SeqDic(lines_list) 
    Col_Dic = {}
    SameLength(dic)
    Key1 = dic.keys()[0]
    for i in range(0,len(dic[Key1])):
        Col_Dic[str(i)]=dic[Key1][i]
    for k in dic.keys()[1:]:	
        for i in range(0,len(dic[k])):
            Col_Dic[str(i)]=Col_Dic[str(i)]+dic[k][i]
    return Col_Dic		 	

	
