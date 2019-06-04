from Bio import SeqIO

# count = 0
# for rec in SeqIO.parse("SRR5832185.fastq", "fastq"):
#     count += 1
# print("%i reads" % count)

good_reads = (rec for rec in \
              SeqIO.parse("SRR5832186.fastq", "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 30)
count = SeqIO.write(good_reads, "good_quality_SRR5832186.fastq", "fastq")
print("Saved %i reads" % count)

trimmed_SRR5832185 = (rec[13:85] for rec in SeqIO.parse("good_quality_SRR5832186.fastq", "fastq"))
count = SeqIO.write(trimmed_SRR5832185, "trimmed_SRR5832186.fasta", "fasta")
print("Saved %i reads" % count)