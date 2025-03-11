# use as python symplicity.py <in.sam> <pval_cutoff> <FC_cutoff>

import sys
import re
import pandas as pd
from scipy.stats import poisson
import random

infile_name = sys.argv[1]
pvalCutoff = float(sys.argv[2])
FCCutoff = float(sys.argv[3])

rootname=infile_name[:-4]
outfile_name=rootname+"_noDupl.sam"

# functions
# get_end
def get_end(startofalig, cigar):
    cignum=re.split('[A-Z]', cigar)[0:-1]
    ciglet=[i for i in cigar if not i.isdigit()]
    readspan=sum([int(j) for (j,i) in zip(cignum, ciglet) if i in {"M", "N", "D"}])
    endofalig=int(startofalig)+readspan-1
    return endofalig

# get_readexons_pos
def get_readexons_pos(startofalig, cigar):
    pos=int(startofalig)
    cignum=re.split('[A-Z]', cigar)[0:-1]
    ciglet=[i for i in cigar if not i.isdigit()]
    readexons=list()
    exon_length=0
    for (j,i) in zip(cignum, ciglet):
        if i in {"M"}:
            exon_length=exon_length+int(j)
        elif i in {"N", "D"}:
            readexons.append([pos, "start"])
            readexons.append([pos+exon_length-1, "end"])
            pos=pos+exon_length+int(j)
            exon_length=0
    readexons.append([pos, "start"])
    readexons.append([pos+exon_length-1, "end"])
    return readexons

# get_splice_sites
def get_splice_sites(startofalig, cigar):
    donors=list()
    acceptors=list()
    pos=int(startofalig)
    cignum=re.split('[A-Z]', cigar)[0:-1]
    ciglet=[i for i in cigar if not i.isdigit()]
    readexons=list()
    exon_length=0
    for (j,i) in zip(cignum, ciglet):
        if i in {"M", "D"}:
            exon_length=exon_length+int(j)
        elif i == "N":
            acceptors.append(pos)
            donors.append(pos+exon_length-1)
            pos=pos+exon_length+int(j)
            exon_length=0
    acceptors.append(pos)
    donors.append(pos+exon_length-1)
    donors=set(donors[:-1])
    acceptors=set(acceptors[1:])
    return [donors, acceptors]

# get_adjcov
def get_adjcov(in_lol, in_extremity):
    # makes a dataframe with upstream and downstream adjacent coverage
    readexon_coords=list()
    for read in in_lol:
        for pos in get_readexons_pos(read[3], read[4]):
            readexon_coords.append(pos)
    posdf = pd.DataFrame(readexon_coords, columns=['Position', 'Type'])

    start_counts=posdf.loc()[posdf["Type"]=="start"].groupby("Position").size().reset_index(name="scount")
    end_counts=posdf.loc()[posdf["Type"]=="end"].groupby("Position").size().reset_index(name="ecount")
    term_counts=pd.DataFrame({"Position": end_counts["Position"]+1, "tcount": end_counts["ecount"]})
    adjuptostarts=pd.DataFrame({"Position": start_counts["Position"]-1})
    counts_df=start_counts.merge(end_counts, how="outer").merge(term_counts, how="outer").merge(adjuptostarts, how="outer").fillna(0).sort_values(by=["Position"]).reset_index(drop=True)
    counts_df["cov"]=(counts_df["scount"].cumsum())-(counts_df["tcount"].cumsum())
    counts_df["prevcov"]=counts_df["cov"].shift()
    counts_df["nextcov"]=counts_df["cov"].shift(-1)

    # get adj cov of start doublets and end doublets
    sdoublet_df=counts_df[counts_df["scount"]>1]
    scov_df=pd.DataFrame({"Start": sdoublet_df["Position"],  "StartCov": sdoublet_df["cov"], "PrevCov": sdoublet_df["prevcov"]})
    edoublet_df=counts_df[counts_df["ecount"]>1]
    ecov_df=pd.DataFrame({"End": edoublet_df["Position"], "EndCov": edoublet_df["cov"], "NextCov": edoublet_df["nextcov"]})

    if in_extremity=="start":
        return scov_df
    elif in_extremity=="end":
        return ecov_df

# get_Cov
def get_Cov(Start, End, PrevCov, NextCov, spl_acc, spl_don):
    if Start in spl_acc and End in spl_don:
        Cov=pd.NA
    elif Start in spl_acc:
        Cov=NextCov
    elif End in spl_don:
        Cov=PrevCov
    else:
        Cov=max([PrevCov, NextCov])
    return Cov

# get_Epl
def get_Epl(cov, readseq):
    ASD=cov/MAL
    f=dict(AL_frq)[len(readseq)]
    EPL=ASD*f
    return EPL

# get_pcount
def get_pcount(obscount, explicity):
    if pd.isna(explicity):
        pcount=1
    else:
        pcount = 1-poisson.cdf(k=(obscount-1), mu=explicity)
    return pcount

# get_DeltaPlic
def get_DeltaPlic(Obs_plic, Exp_plic):
    if pd.isna(Exp_plic):
        DeltaPlic=0
    else:
        DeltaPlic=Obs_plic - round(Exp_plic)
    return DeltaPlic

# get_minDeltaCov
def get_minDeltaCov(start, end, prevcov, startcov, endcov, nextcov, spl_acc, spl_don):
    startNotch=startcov-prevcov
    endNotch=endcov-nextcov
    if start in spl_acc and end in spl_don:
        minDeltaCov=pd.NA
    elif start in spl_acc:
        minDeltaCov=endNotch
    elif end in spl_don:
        minDeltaCov=startNotch
    else:
        minDeltaCov=min(startNotch, endNotch)
    return minDeltaCov

# get_minCovFC
def get_minCovFC(start, end, prevcov, startcov, endcov, nextcov, spl_acc, spl_don):
    def safe_div(x, y, alternate):
        if y == 0:
            return alternate
        return x / y

    startFC=safe_div(startcov, prevcov, float('inf'))
    endFC=safe_div(endcov, nextcov, float('inf'))
    if start in spl_acc and end in spl_don:
        minCovFC=pd.NA
    elif start in spl_acc:
        minCovFC=endFC
    elif end in spl_don:
        minCovFC=startFC
    else:
        minCovFC=min(startFC, endFC)
    return minCovFC

# get_n_to_remove
def get_n_to_remove(oplicity, mindeltacov, deltaplic):
    # define n
    n=0
    if pd.isna(mindeltacov):
        mindeltacov=0
    if pd.isna(deltaplic):
        deltaplic=0
    if min(deltaplic, mindeltacov)>0:
        if oplicity==2:
            n=1
        elif oplicity>2 and mindeltacov>deltaplic:
            n=int(deltaplic)
        elif oplicity>2 and mindeltacov<deltaplic:
            n=int(mindeltacov)
    return n

# get_samlol
def get_samlol(in_losam, in_strand):
    pos_lol=list()
    neg_lol=list()
    for line in in_losam:
        readID=line.split("\t")[0]
        samflag=line.split("\t")[1]
        chrom=line.split("\t")[2]
        startpos=line.split("\t")[3]
        cigar=line.split("\t")[5]
        SEQ=line.split("\t")[9]
        saminfo=[readID, samflag, chrom, startpos, cigar, SEQ]
        binsa=str(bin(4096+int(samflag)))
        if binsa[-5]=="0":
            pos_lol.append(saminfo)
        elif binsa[-5]=="1":
            neg_lol.append(saminfo)
    if in_strand=="pos":
        return pos_lol
    elif in_strand=="neg":
        return neg_lol

# get_q
def get_q(in_vect):
    sorted_ps=sorted(in_vect)
    l=len(in_vect)
    prev_bh_value = 0
    bh_list=list()
    for i in range(l):
        bh_value = sorted_ps[i] * l / (i + 1)
        bh_value = min(bh_value, 1)
        bh_value = max(bh_value, prev_bh_value)
        bh_list.append(bh_value)
        prev_bh_value = bh_value

    ps_df=pd.DataFrame({"ps": in_vect})
    ps_df=ps_df.sort_values("ps")
    ps_df["qs"]=bh_list
    ps_df=ps_df.sort_index()
    qs_orig_order=ps_df["qs"].values.tolist()

    return qs_orig_order

# ban_n_reads
def ban_n_reads(Reads, n):
    banned_reads=random.sample(Reads, n)
    return banned_reads

# get_multiplet_info
def get_multiplet_info(in_lol):
    #The splice donors and acceptors are extracted
    spl_donors=set()
    spl_acceptors=set()
    for alig in in_lol:
        startpos=int(alig[3])
        cigar=alig[4]
        donacc=get_splice_sites(startpos, cigar)
        spl_donors=spl_donors.union(donacc[0])
        spl_acceptors=spl_acceptors.union(donacc[1])

    # multiplet df is created and multiplicates are identified
    possam_df = pd.DataFrame(in_lol, columns=['readID', 'samflag', 'Chrom', 'Start', 'cigar', 'readSeq'])
    multiplet_df  =possam_df.groupby(["readSeq", "Start", "cigar"]).agg(ReadIDs=("readID",lambda x: list(x)), plicity=('readID', 'size')).reset_index()
    multiplet_df = multiplet_df[multiplet_df["plicity"]>1]
    multiplet_df = multiplet_df.astype({"Start": int})

    # the function is interrupted and returns empty df if no multiplicates are found
    empty_df=pd.DataFrame(columns=["ReadIDs", "pcount", "minCovFC", "n_to_remove"])

    if len(multiplet_df.index)==0:
        return empty_df

    # end of alignment is added
    multiplet_df["End"] = multiplet_df.apply(lambda x: get_end(x["Start"], x["cigar"]), axis=1)

    # Start, End and Adjacent coverage are added
    startcov_df=get_adjcov(in_lol=in_lol, in_extremity="start")
    endcov_df=get_adjcov(in_lol=in_lol, in_extremity="end")
    multiplet_df=multiplet_df.merge(startcov_df, how="left", on="Start")
    multiplet_df=multiplet_df.merge(endcov_df, how="left", on="End")

    # coverage and plicity parameters are calculated
    multiplet_df["Cov"]=multiplet_df.apply(lambda x: get_Cov(Start=x["Start"], End=x["End"], PrevCov=x["PrevCov"], NextCov=x["NextCov"], spl_don=spl_donors, spl_acc=spl_acceptors), axis=1)
    multiplet_df["EPL"]=multiplet_df.apply(lambda x: get_Epl(cov=x["Cov"], readseq=x["readSeq"]), axis=1)
    multiplet_df["pcount"]=multiplet_df.apply(lambda x: get_pcount(obscount=x["plicity"], explicity=x["EPL"]), axis=1)

    # coverage and plicity parameters are calculated
    multiplet_df["DeltaPlic"]=multiplet_df.apply(lambda x: get_DeltaPlic(Obs_plic=x["plicity"], Exp_plic=x["EPL"]), axis=1)
    multiplet_df["minDeltaCov"]=multiplet_df.apply(lambda x: get_minDeltaCov(start=x["plicity"], end=x["End"], prevcov=x["PrevCov"], startcov=x["StartCov"], endcov=x["EndCov"], nextcov=x["NextCov"], spl_don=spl_donors, spl_acc=spl_acceptors), axis=1)
    multiplet_df["minCovFC"]=multiplet_df.apply(lambda x: get_minCovFC(start=x["Start"], end=x["End"], prevcov=x["PrevCov"], startcov=x["StartCov"], endcov=x["EndCov"], nextcov=x["NextCov"], spl_don=spl_donors, spl_acc=spl_acceptors), axis=1)

    # reads to remove are counted and essential info is returned
    multiplet_df["n_to_remove"]=multiplet_df.apply(lambda x: get_n_to_remove(oplicity=x["plicity"], mindeltacov=x["minDeltaCov"], deltaplic=x["DeltaPlic"]), axis=1)
    multiplet_df=multiplet_df[multiplet_df["n_to_remove"]>0]

    multiplet_df=multiplet_df[["ReadIDs", "pcount", "minCovFC", "n_to_remove"]]
    multiplet_df=multiplet_df.astype({"minCovFC": float})
    multiplet_df=multiplet_df.astype({"n_to_remove": int})

    return multiplet_df


# get_unwanted
def get_unwanted(in_df, in_FC_cutoff, in_q_cutoff):
    in_df["qval"]=get_q(in_df["pcount"].values.tolist())
    sig_df=in_df[in_df["minCovFC"]>in_FC_cutoff]
    sig_df=sig_df[sig_df["qval"]<in_q_cutoff]
    sig_df["banned_reads"]=sig_df.apply(lambda x: ban_n_reads(Reads=x["ReadIDs"], n=x["n_to_remove"]), axis=1)
    banned_lol=sig_df["banned_reads"].values.tolist()
    banned_set=set([x for y in banned_lol for x in y])
    return banned_set



# files pre-processing
with  open(infile_name, 'r') as infile:
    alsizes=list()
    for line in infile:
          if line[0]!="@":
              alsize=len(line.split("\t")[9])
              alsizes.append(alsize)
AL_frq=[(l, alsizes.count(l) / len(alsizes)) for l in set(alsizes)]
MAL=sum(alsizes)/len(alsizes)

infile.close()

# sam processing
infile = open(infile_name, 'r')
outfile = open(outfile_name, 'w')
logfile = open('log_anal.txt', 'w')
ends=set()
losam=[]
plicity_issue=False
prevstart=int()
prevchrom=str("")
prev_endmax=int(0)
pos_mp_info_df=pd.DataFrame(columns=["ReadIDs", "pcount", "minCovFC", "n_to_remove"])
neg_mp_info_df=pd.DataFrame(columns=["ReadIDs", "pcount", "minCovFC", "n_to_remove"])
for line in infile:
    if line[0]!="@":
        start=int(line.split("\t")[3])
        if start==prevstart:
            plicity_issue=True
        cigar=line.split("\t")[5]
        end=get_end(start, cigar)
        ends.add(end)
        endmax=max(ends)
        chrom=line.split("\t")[2]

        if (start<prev_endmax and chrom==prevchrom) or prevchrom=="":
            losam.append(line)

        else:
            if plicity_issue:
                # Gets multiplicate info df
                pos_mp_info_chunk_df=get_multiplet_info(in_lol=get_samlol(losam, "pos"))
                neg_mp_info_chunk_df=get_multiplet_info(in_lol=get_samlol(losam, "neg"))

                # Appends chunk info_df to main info_df
                pos_mp_info_df=pd.concat([pos_mp_info_df, pos_mp_info_chunk_df], axis=0, ignore_index=True).reset_index(drop=True)
                neg_mp_info_df=pd.concat([neg_mp_info_df, neg_mp_info_chunk_df], axis=0, ignore_index=True).reset_index(drop=True)

            # Resets plicity_issue
            plicity_issue=False

            # Empties ends and adds back end of first read of upcoming block
            ends=set([end])
            endmax=max(ends)

            # Empties losam and adds first read of upcoming block
            losam=[line]

        prevstart=start
        prevchrom=chrom
        prev_endmax=endmax


# after last read
if plicity_issue:
  # Gets multiplicate info df
  pos_mp_info_chunk_df=get_multiplet_info(in_lol=get_samlol(losam, "pos"))
  neg_mp_info_chunk_df=get_multiplet_info(in_lol=get_samlol(losam, "neg"))

  # Appends chunk info_df to main info_df
  pos_mp_info_df=pd.concat([pos_mp_info_df, pos_mp_info_chunk_df], axis=0, ignore_index=True).reset_index(drop=True)
  neg_mp_info_df=pd.concat([neg_mp_info_df, neg_mp_info_chunk_df], axis=0, ignore_index=True).reset_index(drop=True)

# makes false discovery rate correction and extracts reads to remove
pos_unwanted=get_unwanted(pos_mp_info_df, in_FC_cutoff=FCCutoff, in_q_cutoff=pvalCutoff)
neg_unwanted=get_unwanted(neg_mp_info_df, in_FC_cutoff=FCCutoff, in_q_cutoff=pvalCutoff)
unwanted=pos_unwanted.union(neg_unwanted)
infile.close()

infile = open(infile_name, 'r')
for line in infile:
    if line[0]=="@":
        outfile.write(line)

    else:
        readid=line.split("\t")[0]
        if readid not in unwanted:
            outfile.write(line)

outfile.close()
logfile.close()
