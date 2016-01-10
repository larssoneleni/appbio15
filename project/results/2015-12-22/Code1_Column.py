#!/usr/bin/env python
#coding: utf8

##
# Takes a file. Checks if we have a FASTA file. Takes tke name of the sequnce and place it as a key in the directory. The value in the 
# directory is the sequence that responds to the key name. 
##
def SeqDic(lines_list):
    dic = {}
    Split_list = lines_list.split('>')		    
    for line in lines_list[1:]:
    	if line[0]=='>' and len(line.split())>1:
            list_line=line.split()
            dic[list_line[0][1:]=''
	    key = list_line[0][1:]	
	    for i in range(0,len(Split_lis)):
		if len(re.findall(key,Split_list[i])>0:
		    seq_list = Split_list[i].split('\n')[1:-1]
		    for j in seq_list:
		        dic[key]=dic[key]+seq_list[i]
		else:
		    pass
	elif line[0]=='>':
            dic[line[0][1:-1]=''
	    key = line[0][1:-1]	
	    for i in range(0,len(Split_lis)):
		if len(re.findall(key,Split_list[i])>0:
		    seq_list = Split_list[i].split('\n')[1:-1]
		    for j in seq_list:
		        dic[key]=dic[key]+seq_list[i]
		else:
		    pass	
	else:
            pass   	         
    return dic
##
# Take a dictonary and takes the i:th argument in every value and puts it in a new directonary where the key is the a number from 0 to the # length of the sequence.
## 
def column(dic):
    Col_Dic = {}
    length(dic)	
    key1 = dic.keys()[0]
    for i in range(0,len(dic[key1])):
	Col_Dic[i]=dic[i][i]
	key = i
        for k in dic.keys()[1:]:
	    for j in range(0,len(dic[k])):
                if j == key:
		    Col_Dic[key]=Col_Dic[key]+dic[k][j]
	        else:
		    pass
    return Col_Dic		 	
		
    		

##
# Takes a Dictonary and checks if the values are in the same lenght if not gives an error.
##
def length(dic):
    key1 = dic.keys()[0] # takes out the first key in the dictonary
    for k in dic:
        if len(dic[key1])==len(dic[k]):
	    pass
	else:
	    sys.stdout.write('ERROR: The sequnces have not the same lenght')
	    break	  	      		


