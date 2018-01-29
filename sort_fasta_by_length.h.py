import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'rU')
Seqdict=SeqIO.to_dict(SeqIO.parse(FastaFile, 'fasta'))
seqlengths = [(x,len(Seqdict[x].seq)) for x in Seqdict]
sorted_seqlengths = sorted(seqlengths,key=lambda x: x[1],reverse=True)
sys.stderr.write("\n".join(["{}: {}".format(x[0],x[1]) for x in sorted_seqlengths[:5]]))
for contig in sorted_seqlengths:
	SeqIO.write(Seqdict[contig[0]],sys.stdout,'fasta')

#for rec in SeqIO.parse(FastaFile, 'fasta'):
#    name = rec.id
#    seq = rec.seq
#    seqLen = len(rec)
#    print name, seqLen

#FastaFile.close()
