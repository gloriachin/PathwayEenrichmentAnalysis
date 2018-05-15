# coding: utf-8

from scipy import stats

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector


stats_r = importr('stats')

import sys
if sys.argv[3]=='cp':
    fin = open('MsigDB_c2.cp.v5.1.symbols.gmt','r')
if sys.argv[3]=='halmark':
    fin = open('h.all.v5.1.symbols.gmt','r')
if sys.argv[3]=='KEGG':
    fin = open('KEGG_pathway_genes.gmt','r')
if sys.argv[3]=='KEGG_signal':
    fin = open('KEGG_Gene_new_signal.txt.strip.txt')
if sys.argv[3]=='MSigDB_KEGG':
    fin = open('MSigDB_KEGG_updata.txt')
if sys.argv[3]=='CHR':
    fin = open('c1.all.v5.1.symbols.gmt.txt', 'r')
if sys.argv[3]=='mouse_KEGG':
    fin = open('mouse_kegg_geneset_full.txt','r')
dic_reac_uni = {}
all_unip = set()

for line in fin.readlines():
    single_unip = set()
    words = line.strip().split('\t')
    key = words[0].strip()
    length = len(words)
    for i in range(2,length):
        all_unip.add(words[i])
        single_unip.add(words[i])
    dic_reac_uni[key] = single_unip
fin.close()

all_u = len(all_unip)
print("all genes in background datasets:")
print(all_u)

import re
input_file = open(sys.argv[1],'r')
input_list = set()
for line in input_file.readlines():
     input_list.add((line.strip().split('\t')[int(sys.argv[2])]).strip('"'))
all_match = 0
match_list_all=[]
for uni in input_list:
    if uni in all_unip:
        match_list_all.append(uni)
        all_match = all_match + 1
print("all matched uniprot")
print(all_match)

dic_match = {}
for pathway in dic_reac_uni:
    count = 0
    match_item=set()
    for word in input_list:
        if word in dic_reac_uni[pathway]:
            count = count + 1
            match_item.add(word)
    dic_match[pathway]=[count, len(dic_reac_uni[pathway]),','.join(match_item)]

ave_match = float(all_match)/all_u

if sys.argv[3]=='cp':
    PRE='CP_'
if sys.argv[3]=='halmark':
    PRE='Halmark_'
if sys.argv[3]=='KEGG':
    PRE='KEGG_'
if sys.argv[3]=='KEGG_signal':
    PRE='KEGG_signal_'
if sys.argv[3]=='MSigDB_KEGG':
    PRE='MSigDB_KEGG_'
if sys.argv[3]=='CHR':
    PRE='CHR'
if sys.argv[3]=='mouse_KEGG':
    PRE='MOUSE_KEGG_'
fou1 = open(PRE+sys.argv[1]+'.fdr.txt','w')
fou1.write('Pathway\tMatched_gene\tAll_genes\tRatio\tP-value\tFDR\tLable\tGenes\n')
result = {}
Name_list = []
un_Name_list = []
FDR_gene_list = []
for i in dic_match:
    if dic_match[i][0] != 0:
        Name_list.append(i)
        table = [(dic_match[i][0]),(dic_match[i][1]),(all_match)-(dic_match[i][0]),(all_u-dic_match[i][1])]
        p = stats.fisher_exact([[table[0],table[1]],[table[2],table[3]]],alternative='greater')
        FDR_gene_list.append(p[1])
        result[i]=[dic_match[i][0],dic_match[i][1],float(dic_match[i][0])/float(dic_match[i][1]),p[1]]
    else:
        un_Name_list.append(i)
        result[i]=[dic_match[i][0],dic_match[i][1],float(dic_match[i][0])/float(dic_match[i][1]),'NA']

FDR_gene_list_adjust = stats_r.p_adjust(FloatVector(FDR_gene_list), method='BH')
n = 0
for i in Name_list:
        fou1.write(i+'\t')
        for m in result[i]:
            fou1.write(str(m)+'\t')
        fou1.write(str(FDR_gene_list_adjust[n])+'\t')
        if result[i][2] > ave_match:
            fou1.write('Enriched\t')
        else:
            fou1.write('BelowAverage\t')
        fou1.write(dic_match[i][2]+'\t')
        fou1.write("\n")
        n = n+1
fou1.close()
