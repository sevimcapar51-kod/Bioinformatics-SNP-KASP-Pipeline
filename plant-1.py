from Bio import SeqIO

genomic = SeqIO.read("gene_region.fasta", "fasta").seq

best_protein = ""
best_seq = ""

for strand, seq in [("+", genomic), ("-", genomic.reverse_complement())]:
    for frame in range(3):
        translated = seq[frame:].translate()
        
        # stoplara böl
        proteins = translated.split("*")
        
        for p in proteins:
            if len(p) > len(best_protein):
                best_protein = p
                best_seq = seq

print("Best protein length:", len(best_protein))

from Bio import SeqIO

genomic = SeqIO.read("gene_region.fasta", "fasta").seq

best_protein = ""
best_cds = ""

for strand, seq in [("+", genomic), ("-", genomic.reverse_complement())]:
    for frame in range(3):
        translated = seq[frame:].translate()
        proteins = translated.split("*")
        
        start = 0
        for p in proteins:
            if len(p) > len(best_protein):
                best_protein = p
                best_cds = seq[frame+start:frame+start+len(p)*3]
            start += (len(p)+1)*3

print("CDS length:", len(best_cds))
print("Protein length:", len(best_protein))

from Bio import SeqIO
from Bio.Align import PairwiseAligner

ref = SeqIO.read("Solyc02g084610_CDS.fasta", "fasta").seq
wild = best_cds  # az önce bulduğun CDS

aligner = PairwiseAligner()
aligner.mode = 'global'

alignment = aligner.align(ref, wild)[0]

ref_aln = alignment.aligned[0]
wild_aln = alignment.aligned[1]

# basit karşılaştırma
snps = []
for i in range(min(len(ref), len(wild))):
    if ref[i] != wild[i]:
        snps.append((i, ref[i], wild[i]))

print("REAL SNP count:", len(snps))
print(snps[:20])

clean_snps = [(pos, ref, alt) for pos, ref, alt in snps if ref.upper() != alt.upper()]

filtered = []

for i, (pos, ref, alt) in enumerate(clean_snps):
    isolated = True
    for j, (pos2, _, _) in enumerate(clean_snps):
        if i != j and abs(pos - pos2) < 10:
            isolated = False
            break
    if isolated:
        filtered.append((pos, ref, alt))

print("Isolated SNPs:", filtered)

from Bio import SeqIO

ref = SeqIO.read("Solyc02g084610_CDS.fasta", "fasta").seq

pos = 204

left = ref[pos-25:pos]
snp = ref[pos]
right = ref[pos+1:pos+26]

print("Context:")
print(left + "[" + snp + "/C]" + right)

from Bio import SeqIO
from Bio.Seq import Seq

record = SeqIO.read("Solyc02g084610_CDS.fasta", "fasta")
ref = record.seq

# SNP pozisyonu (0-based!)
pos = 204

# referans protein
ref_prot = ref.translate()

# mutant oluştur
alt_list = list(ref)
alt_list[pos] = "C"

# tekrar Seq objesine çevir
alt_seq = Seq("".join(alt_list))

alt_prot = alt_seq.translate()

# karşılaştır
for i in range(len(ref_prot)):
    if ref_prot[i] != alt_prot[i]:
        print("AA change at:", i, ref_prot[i], "→", alt_prot[i])