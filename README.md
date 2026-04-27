# Bioinformatics-Pipeline-for-SNP-Detection-and-KASP-Marker-Preparation
## Note
This project presents an ongoing bioinformatics pipeline for SNP detection and KASP marker development in tomato. It was developed to explore and implement key steps including sequence alignment, SNP identification, filtering, and marker candidate generation. The workflow combines computational analysis with biological interpretation and is continuously being improved.

# Description
This project was developed to detect Single Nucleotide Polymorphisms (SNPs) in the tomato genome (Solanum lycopersicum and Solanum arcanum) and to generate KASP (Competitive Allele-Specific PCR) marker candidates based on these SNPs. The workflow includes translating CDS sequences to proteins, genomic alignment, SNP extraction, and genomic context analysis. Core bioinformatics operations were performed using Biopython, supported by terminal-based tools such as BLAST and samtools. The project presents a prototype pipeline merging manual analysis with automation.

# Objectives
This project was developed with the following goals:
To perform SNP detection for a specific gene in the tomato genome.
To evaluate the impact of SNPs at the protein level.
To create KASP marker candidates from the SNP flanking regions.
To establish a basic marker development pipeline using bioinformatics tools.

# Technologies Used
Biopython
BLAST
samtools

# Installation

### Homebrew installation (macOS)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

### Setting the PATH
echo 'eval "$(/opt/homebrew/bin/brew shellenv zsh)"' >> /Users/macbook/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv zsh)"

### BLAST installation
brew install blast

### samtools installation
brew install samtools

# Execution

## Creating BLAST database
makeblastdb -in S_arcanum_genome.fna -dbtype nucl -out arcanum_db

## Running BLAST (protein vs genome)

## Extracting genomic region
samtools faidx S_arcanum_genome.fna CM050934.1:52329256-52330500 > gene_region.fasta

## Running the Python script
python3 plant-1.py

# Sample Usage

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

## Reading the CDS
record = SeqIO.read("Solyc02g084610_CDS.fasta", "fasta")

## Translating to protein
protein = record.seq.translate(to_stop=True)

print("Protein length:", len(protein))
print("First 100 aa:")
print(protein[:100])

## Saving as FASTA
protein_record = SeqRecord(
    protein,
    id="Solyc02g084610_protein",
    description="(RLP) protein translated from CDS"
)

SeqIO.write(protein_record, "protein.fasta", "fasta")

## Project Structure
plaintext
work_1/
│
├── plant-1.py
├── protein.fasta
├── Solyc02g084610_CDS.fasta
├── S_arcanum_genome.fna
├── S_arcanum_genome.fna.fai
├── gene_region.fasta
├── results.out
├── dizi.fasta
├── primer_blast.out
│
├── arcanum_db.nhr
├── arcanum_db.nin
├── arcanum_db.nsq
├── arcanum_db.ndb
├── arcanum_db.not
├── arcanum_db.njs
├── arcanum_db.ntf
└── arcanum_db.nto

# Troubleshooting (Errors and Solutions)

### Error: zsh: command not found: brew

Cause: Homebrew is not added to the PATH.

Solution: eval "$(/opt/homebrew/bin/brew shellenv zsh)"

### Error: command not found: makeblastdb

Cause: BLAST is not installed.

Solution: brew install blast

### Error: zsh: command not found: samtools

Cause: samtools is not installed.

Solution: brew install samtools

### Error: Error: Argument "-outfmt". Value is missing

Cause: The -outfmt parameter is missing in the BLAST command.

Solution: blastn -query dizi.fasta -db arcanum_db -out primer_blast.out -outfmt 6

### Error: BLAST Database error: No alias or index file found

Cause: Running in the wrong directory or incorrect database path.

Solution:cd /Users/macbook/work_1
blastn -query dizi.fasta -db arcanum_db -out primer_blast.out -outfmt 6

### Error: [E::fai_build3_core] Failed to open the file S_arcanum_genome.fna : No such file or directory

Cause: Running in the wrong folder.

Solution: cd /Users/macbook/work_1

### Error: TypeError: str.translate() takes exactly one argument (0 given)

Cause: Incorrect usage of the Python string method (instead of Biopython).

Solution: alt_prot = alt_seq.translate()

### Error: IndexError: list index out of range

Cause: Alignment parsing process performed incorrectly.

Solution: (Missing Information – the correct parsing method is not yet defined in the project).

## Learnings / Notes

CDS → protein translation can be reliably performed with Biopython.
Genomic locus detection via BLAST is possible, but strand information is of critical importance.
SNP → genomic coordinate conversion requires precision (especially on the reverse strand).
Lowercase sequences usually represent masked/repeat regions.
In KASP marker design, the genomic context is as decisive as the SNP itself.
Short primer sequences may yield multiple hits in BLAST (specificity issue).

## Next Steps

Integrate exon-intron aware mapping using GFF/GTF annotations.
Implement a more robust variant analysis for SNP validation.
Add thermodynamic analysis (Tm, GC%) to primer design.
Develop an automated SNP → KASP pipeline.

## Future Development & Roadmap

Full parameters for BLAST commands (protein alignment part).
Full code for the SNP extraction algorithm.
Final correct version of the alignment parsing.
Detailed implementation of the primer design algorithm.
KASP marker validation steps (in silico PCR or experimental validation).

### Developer: Sevim ÇAPAR
### Molecular Biologist | Biotechnology MSc. Researcher
