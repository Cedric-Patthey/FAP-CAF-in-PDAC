import sys
from Bio import SeqIO
import gzip

fastq_file = sys.argv[1]  # Input fastq file
ids_file = sys.argv[2] # Input file, one per line
result_file = sys.argv[3] # Output fastq file

unwanted = set()
with open(ids_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            unwanted.add(line)

fastq_sequences = SeqIO.parse(gzip.open(fastq_file, 'rt'),'fastq')
end = False
with open(result_file, "w") as f:
    for seq in fastq_sequences:
        if seq.id not in unwanted:
            SeqIO.write([seq], f, "fastq")
