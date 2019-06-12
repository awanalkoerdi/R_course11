from Bio import SeqIO

# count = 0
# for rec in SeqIO.parse("SRR5832185.fastq", "fastq"):
#     count += 1
# print("%i reads" % count)

good_reads = (rec for rec in \
              SeqIO.parse("SRR5832185.fastq", "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 30)
count = SeqIO.write(good_reads, "good_quality_SRR5832185.fastq", "fastq")
print("Saved %i reads" % count)

trimmed_SRR5832185 = (rec[13:96] for rec in SeqIO.parse("good_quality_SRR5832185.fastq", "fastq"))
count = SeqIO.write(trimmed_SRR5832185, "trimmed_SRR5832185.fastq", "fastq")
print("Saved %i reads" % count)