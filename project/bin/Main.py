#coding: utf8
import sys
import os
from os import listdir
import dendropy
from dendropy import Tree, TaxonNamespace
from ReducingNoise import TakeAway
from dendropy.calculate import treecompare
import re
from column import SeqDic, order
import tempfile
import math
##
# Reading in directories
##
file_path_list=sys.argv[1:]
##
# Counting the times the referense tree is recovered. For the none reducing file
##
NotFixedCount = 0
##
# Counting the times the reference tree is recovered for the reducing file.
##
FixedCount = 0
Total = 0
for i in range(len(file_path_list)):
    file_path = file_path_list[i]
    onlyfiles = []
    AllFiles = []
    files = listdir('/home/4/u1we1f44/Documents/Appbio15/projekt/data/'+file_path)
    for index in range(len(files)):
	if len(re.findall('py',files[index]))>0 or len(re.findall('txt',files[index]))>0:
            pass
        else:
            AllFiles.append(files[index])
    for j in range(len(AllFiles)):
        if len(re.findall('tree',AllFiles[j]))>0:
            RefTree = AllFiles[j] 
        else:
            onlyfiles.append(AllFiles[j])
    for j in range(len(onlyfiles)):
        try:
            path = file_path+'/'+onlyfiles[j] 
            fil = open('/home/4/u1we1f44/Documents/Appbio15/projekt/data/'+path, 'r')
            lines_list=fil.readlines()
            fil.close()
            test = SeqDic(lines_list) # If this dose not worke we do not have a FASTA file
            ##
            # Makes a newick tree and checks if the referense tree is recovered. The none reducing file.
            ##
            line = 'cat /home/4/u1we1f44/Documents/Appbio15/projekt/data/'+path+' | fastprot -I fasta -O phylip | fnj -I phylip -O "newick" -o "Treeout.txt"' 
            os.system(line)
            TreePath=file_path+'/'+RefTree   
            t1=Tree.get(file=open('/home/4/u1we1f44/Documents/Appbio15/projekt/data/'+TreePath,'r'),schema="newick",tree_offset=0)
            t2=Tree.get(file=open('/home/4/u1we1f44/Documents/Appbio15/projekt/src/Treeout.txt','r'),schema="newick",tree_offset=0,taxon_namespace=t1.taxon_namespace)
            t1.encode_bipartitions()
            t2.encode_bipartitions()
            if treecompare.symmetric_difference(t1, t2)==0:
                NotFixedCount += 1
                os.remove('/home/4/u1we1f44/Documents/Appbio15/projekt/src/Treeout.txt')
                Total += 1
            else:
                Total += 1
	        os.remove('/home/4/u1we1f44/Documents/Appbio15/projekt/src/Treeout.txt')
            ##  
            # Makes a temporary file. In the temporary file with data with the nosie columns remoeved. MAkes a newick tree and checks if the refernse tree is recovered.
            # The nosie columns removed.
            ##
            os.system("touch temp.fa")
            tempf = open('temp.fa','w')
            seq_dic = TakeAway(lines_list)
            Order = order(lines_list)
            for i in range(len(Order)):
                tempf.write('>'+Order[str(i)]+'\n'+seq_dic[Order[str(i)]]+'\n')
            tempf.close()
            line = 'cat temp.fa | fastprot -I fasta -O phylip | fnj -I phylip -O "newick" -o "ReducingTreeout.txt"'
            os.system(line)
            TreePath=file_path+'/'+RefTree   
            tree1=Tree.get(file=open('/home/4/u1we1f44/Documents/Appbio15/projekt/data/'+TreePath,'r'),schema="newick",tree_offset=0)
            tree2=Tree.get(file=open('/home/4/u1we1f44/Documents/Appbio15/projekt/src/ReducingTreeout.txt','r'),schema="newick",tree_offset=0,taxon_namespace=tree1.taxon_namespace)
            tree1.encode_bipartitions()
            tree2.encode_bipartitions() 
            if treecompare.symmetric_difference(tree1, tree2)==0:
                FixedCount += 1
		os.remove('/home/4/u1we1f44/Documents/Appbio15/projekt/src/ReducingTreeout.txt')
                os.remove('/home/4/u1we1f44/Documents/Appbio15/projekt/src/temp.fa') 	
	    else:
		os.remove('/home/4/u1we1f44/Documents/Appbio15/projekt/src/ReducingTreeout.txt')
                os.remove('/home/4/u1we1f44/Documents/Appbio15/projekt/src/temp.fa')  
	  
	except IOError:
            sys.stdout.write('"{0}/{1}" is not a FASTA file. This program can only take FASTA files.\n'.format(file_path,onlyfiles[j]))	
            quit()

print 'Number of times Reference tree is recovered with none reducing file:'
print NotFixedCount
print 'Number of times Reference tree is recovered with reducing file:'
print FixedCount
print 'Total number of files:'
print Total


