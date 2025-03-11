# use as python make_ddSeq_RUC_table.py <rootname> <r1.fastq.gz> <read_feature.txt>

import sys
import pandas as pd
import gzip
import re
import subprocess
import random

rootname = sys.argv[1]
R1_fastq_file  = sys.argv[2]
read_feature_file = sys.argv[3]

# loads the bc whitelists  
bc_wl=list()
for seq in open("Whitelist_96_bc.txt"):
    seq=seq.strip()
    bc_wl.append(seq)

# calculates the one-off whitelists
def get_1_oneoffs(in_hexamer):
    oneoffs_list=list()
    for i in range(1,7):
        pref=in_hexamer[0:i-1]
        suff=in_hexamer[i:]
        for letter in ['T', 'C', 'G', 'A']:
          oneoffs_list.append(pref+letter+suff)
    return oneoffs_list 

def get_all_oneoffs(in_list):
    all_1offs=list()
    for item in in_list:
        all_1offs=all_1offs+get_1_oneoffs(item)
    return all_1offs

bc_oowl=get_all_oneoffs(bc_wl)


# makes a set of the whitelist
bc_wl_set=set(bc_wl)
bc_oowl_set=set(bc_oowl)


# a function to calculate the distance between two strings
def dist_bc(seq1,seq2):
    seq1 = list(seq1)
    seq2 = list(seq2)
    return sum(map(lambda pair: pair[0] != pair[1], zip(seq1, seq2)))


# function to correct barcode
def correct_bc(inseq):
    if inseq in bc_wl_set:
        return inseq
    
    elif inseq in bc_oowl_set:
        for item in bc_wl:
            d = int(dist_bc(inseq, item))
            if d<=1:
                corbac = item
                break
        return corbac

# R1 is read and barcodes extracted and corrected. Non-matching sequences are saved separately to a fasta file.
RL_bcumi_lol=list()
nomatch_sequences_lol=list()
i=1
pattern = re.compile(r'[ACGTN].....TAG.........TGC......TAC.........GAA......ACG........GAC')
for line in gzip.open(R1_fastq_file, 'rt'):
    fq_nr=(i+3)//4
    if i%4==1:
        headline='R'+str(fq_nr)
    elif i%4==2:
        seq=line.strip()
        match=pattern.search(seq)
        if match:
            bc1=match.group()[0:6]
            bc2=match.group()[21:27]
            bc3=match.group()[42:48]
            umi=match.group()[51:59]
            if bc1 in bc_oowl_set and bc2 in bc_oowl_set and bc3 in bc_oowl_set:
                # corrects barcodes
                bc1=correct_bc(bc1)
                bc2=correct_bc(bc2)
                bc3=correct_bc(bc3)
                bc=bc1+bc2+bc3
                RL_bcumi_lol.append([headline, bc, umi])
        else:
            nomatch_sequences_lol.append([headline, seq])
    i=i+1

# nomatch sequences are saved to fasta, including sequence in header
with open(rootname+"_nomatch.fas", 'w') as f:
    for item in nomatch_sequences_lol:
        headline=item[0]
        seq=item[1]
        seqtail=seq[40:]
        f.write(">"+headline+"@"+seqtail+"\n"+seq+"\n")
        
# a subset of barcodes is sampled to generate top 1000 models
subset_bc_list=random.sample([entry[1] for entry in RL_bcumi_lol], 100000)

# function to output a list of top n barcodes
def top_n_strings(input_list, n):
    df = pd.DataFrame(input_list, columns=['strings'])
    frq = df['strings'].value_counts()
    return frq.head(n).index.tolist()

top_1000_bc_list=top_n_strings(subset_bc_list, 1000)
top_1000_bc_set=set(top_1000_bc_list)

# barcodes belonging to the top 1000 most frequent are selected
RL_bcumi_lol=[entry for entry in RL_bcumi_lol if entry[1] in top_1000_bc_set]

# model sequences are saved to fasta
with open(rootname+'_models.fas', 'w') as g:
    i=int(1)
    for bc in top_1000_bc_list:
        g.write(">Cell_"+str(i)+"@"+bc+"\n"+bc[0:6]+"TAGCCATCGCATTGC"+bc[6:12]+"TACCTCTGAGCTGAA"+bc[12:18]+"ACG\n")
        i=i+1

# runs cd-hit
process = subprocess.run(f"cd-hit-est-2d -i {rootname}_models.fas -i2 {rootname}_nomatch.fas -o {rootname}_nomatch_cdhit.fas -M 0 -T 10 -n 10 -d 0 -p 1 -c 0.93 -s2 0.75 -S2 20 -g 1", shell=True, check=True)

# nomatch barcodes are salvaged from cd-hit output
salvaged_bcumi_lol=list()
for line in open(rootname+"_nomatch_cdhit.fas.clstr"):
    if '>Cell' in line:
        cell = line.split("@")[0].split(">")[1]
        bc = line.split("...")[0].split("@")[1]

    elif 'Cluster' not in line and '/+/' in line:
        readID = line.split("@")[0].split(">")[1]
        R1seqtail = line.split('...')[0].split('@')[1]
        posinfo = line.split('/+/')[0].split('at ')[1]
        readstart = int(posinfo.split(":")[0])
        readend = int(posinfo.split(":")[1])
        refstart = int(posinfo.split(":")[2])
        refend = int(posinfo.split(":")[3])
        pident = float(line.split('%')[0].split('/+/')[1])

        if (readend-readstart+1)>=50 and (refend-refstart+1)>=50 and pident>=96:
            bord_UMI = R1seqtail[readend-40-3:readend-40+11]
            bord = bord_UMI[0:3]+bord_UMI[-3:]
            UMI = bord_UMI[3:11]

            if dist_bc("ACGGAC", bord)<=1:
                salvaged_bcumi_lol.append([readID, bc, UMI])

# matched and salvaged barcodes are joined and transformed to a dataframe
bc_umi_df=pd.DataFrame(RL_bcumi_lol+salvaged_bcumi_lol, columns=['Read_ID', 'bc', 'umi'])

bc_umi_df.to_csv(path_or_buf=rootname+'_bc_umi_df', sep='\t', index=False)

# read info is loaded from the bedtools output file
read_feature_df=pd.read_table(read_feature_file, names=['Read_ID', 'feature', 'align_length'])
read_feature_df['Gene']=[item.split("_")[1] for item in read_feature_df['feature']]
read_gene_df=read_feature_df[['Read_ID', 'Gene', 'align_length']]
read_gene_df = read_gene_df.sort_values('align_length', ascending=False)
read_gene_df = read_gene_df.drop_duplicates(subset=['Read_ID', 'Gene'])

read_gene_df.to_csv(path_or_buf=rootname+'_read_gene_df', sep='\t', index=False)

# The barcode and gene info are merged into a single table
read_bc_umi_gene_df=pd.merge(bc_umi_df, read_gene_df, on='Read_ID', how='inner')
read_bc_umi_gene_df=read_bc_umi_gene_df.reset_index(drop=True)

read_bc_umi_gene_df.to_csv(path_or_buf=rootname+'_read_bc_umi_gene_df', sep='\t', index=False)

# The UMIs are flattened 
read_bc_gene_df=read_bc_umi_gene_df.drop_duplicates(['bc', 'umi', 'Gene'], keep="first").reset_index(drop=True).drop('umi', axis=1)

read_bc_gene_df.to_csv(path_or_buf=rootname+'_read_bc_gene_df', sep='\t', index=False)

# The number of UMI per bracode-gene pair is caculated
bc_Gene_UC=read_bc_gene_df.groupby(['bc', 'Gene']).size().rename('RUC').reset_index()
bc_Gene_UC.to_csv(path_or_buf=rootname+'_RUC_long.txt', sep='\t', index=False)
