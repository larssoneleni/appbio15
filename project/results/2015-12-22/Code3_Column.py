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
    Split_list = lines_list.split('>')		    
    for line in lines_list[1:]:
    	if line[0]=='>' and len(line.split())>1:
            list_line=line.split()
            dic[list_line[0][1:]]=''
	    key = list_line[0][1:]	
	    for i in range(0,len(Split_list)):
		if len(re.findall(key,Split_list[i]))>0:
		    seq_list = Split_list[i].split('\n')[1:-1]
		    for j in seq_list:
		        dic[key]=dic[key]+seq_list[i]
		else:
		    pass
	elif line[0]=='>':
            dic[line[0][1:-1]]=''
	    key = line[0][1:-1]	
	    for i in range(0,len(Split_list)):
		if len(re.findall(key,Split_list[i]))>0:
		    seq_list = Split_list[i].split('\n')[1:-1]
		    for j in seq_list:
		        dic[key]=dic[key]+seq_list[i]
		else:
		    pass	
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
    for k in dic:	
        for i in range(0,len(dic[k])):
	    Col_Dic[i]=dic[k][i]
    print Col_Dic		 	
		

file_path=sys.argv[1]
msl_file=open('/home/n/u1tw68bn/appbio15/Project/data/asymmetric_1.0/'+file_path,'r')
lines_list=msl_file.readlines()
msl_file.close()

##
# Checks if we have a FASTA file. If not gives an Error.
## 
try:
    input_handle = open('/home/n/u1tw68bn/appbio15/Project/data/asymmetric_1.0/'+file_path, "rU")
    alignments = AlignIO.parse(input_handle, "fasta")
    input_handle.close()
    column(lines_list)
except IOError:
    sys.stdout.write('This is not A FASTA file. This program can only take FASTA files.')
	
