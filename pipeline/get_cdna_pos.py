# use as python get_cdna_pos.py <in.sam> <in_cDNA_pos_table> <in.bed>  <max_read_length> <n_reads>

import sys
import re
import pandas as pd
import random

infile = sys.argv[1]
cDNA_pos_file = sys.argv[2]
bedfile=sys.argv[3]
max_read_length=int(sys.argv[4])
n_reads=int(sys.argv[5])

# a function to calculate the end of a read from start and cigar
def get_end(startofalig, cigar):
    cignum=re.split('[A-Z]', cigar)[0:-1]
    ciglet=[i for i in cigar if not i.isdigit()]
    readspan=sum([int(j) for (j,i) in zip(cignum, ciglet) if i in {"M", "N", "D"}])
    endofalig=int(startofalig)+readspan-1
    return endofalig

# takes the sam entries and extracts coordinates
samlol=list()
for line in open(infile):
    if line[0]!="@":
        readID=line.split("\t")[0]
        chrom=line.split("\t")[2]
        startpos=int(line.split("\t")[3])
        cigar=line.split("\t")[5]
        endpos=get_end(startpos, cigar)
        saminfo=[readID, chrom, startpos, endpos]
        samlol.append(saminfo)

sam_df=pd.DataFrame(samlol, columns=["readID", "Chrom", "Start", "End"])

# samples n reads
sam_df=sam_df.sample(n=n_reads)

# makes a set of pairs with sam start and end coordinates
sam_coord_pairs=set()
for read in samlol:
    sam_coord_pairs.add((read[1], read[2]))
    sam_coord_pairs.add((read[1], read[3]))

# extracts coordinates from cDNA_pos table
cdna_coord_lol=list()
for line in open(cDNA_pos_file):
    if line.split("\t")[0]!="Chrom":
        chrom=line.split("\t")[0]
        genomepos=int(line.split("\t")[1])
        if (chrom, genomepos) in sam_coord_pairs:
            gene=line.split("\t")[3]
            cdnapos=int(line.split("\t")[5])
            cdna_coord_lol.append([chrom, genomepos, gene, cdnapos])

cdna_coord_df=pd.DataFrame(cdna_coord_lol, columns=["Chrom", "genomepos", "Gene", "cdnapos"])


# extracts cDNA length from bedfile
bedlol=list()
for line in open(bedfile):
    start=int(line.split("\t")[1])
    end=int(line.split("\t")[2])
    exonlen=end-start+1
    gene=line.split("\t")[3].split("_")[1]
    bedlol.append([gene, exonlen])

bed_df=pd.DataFrame(bedlol, columns=["Gene", "size"])
genesize_df=bed_df.groupby('Gene')["size"].sum().reset_index()

# merges read coordinatesa and cDNA_pos info
sam_df=sam_df.merge(cdna_coord_df, left_on=['Chrom', 'Start'], right_on=['Chrom', 'genomepos'], how='left').dropna().drop('genomepos', axis=1)
sam_df.columns=sam_df.columns.str.replace('Gene', 'start_gene').str.replace('cdnapos', 'cdna_start')
sam_df=sam_df.merge(cdna_coord_df, left_on=['Chrom', 'End'], right_on=['Chrom', 'genomepos'], how='left').dropna().drop('genomepos', axis=1)
sam_df.columns=sam_df.columns.str.replace('Gene', 'end_gene').str.replace('cdnapos', 'cdna_end')
sam_df=sam_df[["start_gene", "end_gene", "cdna_start", "cdna_end"]]
sam_df=sam_df[sam_df["start_gene"]==sam_df["end_gene"]]
sam_df=sam_df[["start_gene", "cdna_start", "cdna_end"]]
sam_df.columns=sam_df.columns.str.replace('start_gene', 'Gene')

# reads with long distance between start and end are discarded
sam_df=sam_df.loc[(sam_df["cdna_end"]-sam_df["cdna_start"])<max_read_length]

# calculates mean position
def get_mean(x, y):
    z=(x+y)/2
    z=round(z)
    return z
sam_df['mean_pos']=sam_df.apply(lambda x: get_mean(x["cdna_start"], x["cdna_end"]), axis=1)

# calculates relative position and distance to end
sam_df=sam_df.merge(genesize_df, on='Gene', how='left').dropna()
def get_d_to_end(meanpos, totsize):
    dte=totsize-meanpos
    return dte

sam_df['d_to_end']=sam_df.apply(lambda x: get_d_to_end(x["mean_pos"], x["size"]), axis=1)
sam_df['relpos']=round(sam_df['mean_pos']/sam_df['size'], 3)
sam_df=sam_df[["Gene", "relpos", "d_to_end"]]

sam_df.to_csv("sam_cDNA_pos_df", sep='\t', index=False, header=True)
