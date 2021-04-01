import os
import sys
import Bio
import pprint
import collections
from Bio import Entrez
from Bio import GenBank
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
cwd = os.getcwd()

codon_count = {}
CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

SynonymousCodons = { 
 'CYS': ['TGT', 'TGC'], 
 'ASP': ['GAT', 'GAC'], 
 'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
 'GLN': ['CAA', 'CAG'], 
 'MET': ['ATG'], 
 'ASN': ['AAC', 'AAT'], 
 'PRO': ['CCT', 'CCG', 'CCA', 'CCC'], 
 'LYS': ['AAG', 'AAA'], 
 'STOP': ['TAG', 'TGA', 'TAA'], 
 'THR': ['ACC', 'ACA', 'ACG', 'ACT'], 
 'PHE': ['TTT', 'TTC'], 
 'ALA': ['GCA', 'GCC', 'GCG', 'GCT'], 
 'GLY': ['GGT', 'GGG', 'GGA', 'GGC'], 
 'ILE': ['ATC', 'ATA', 'ATT'], 
 'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 
 'HIS': ['CAT', 'CAC'], 
 'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'], 
 'TRP': ['TGG'], 
 'VAL': ['GTA', 'GTC', 'GTG', 'GTT'], 
 'GLU': ['GAG', 'GAA'], 
 'TYR': ['TAT', 'TAC'] 
}

def read_file(filename):
    sequences = []
    descr = None
    with open(filename) as file:
        line = file.readline()[:-1]                     # always trim newline
        while line:
            if line[0] == '>':
                if descr:                               # any sequence found yet?
                    sequences.append((descr, seq))
                descr = line[1:].split('|')
                seq = ''                                # start a new sequence
            else:
                seq += line
            line = file.readline()[:-1]
        sequences.append((descr, seq))                  # easy to forget!; this add the third sequence
    return sequences


AGF = read_file('Piromy1_GeneModels_FilteredModels2_cds.fasta')
#Files for P. indianae genome can be found on Mycocosm (https://genome.jgi.doe.gov/portal/Piromy1/download/Piromy1_GeneCatalog_CDS_20200129.fasta.gz)

seqs = [AGF]
codon_count = CodonsDict
CDS_count = 0
#print (seqs)
#for genomes
for CDS in seqs:
    #print (CDS)
    for i in CDS:
        CDS_count += 1
        #print (i)
        dna_sequence = i[1]
        #print (dna_sequence)
#print (CDS_count)
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            #print (codon)
            if codon in codon_count:
                codon_count[codon] += 1
od = collections.OrderedDict(sorted(codon_count.items()))
pprint.pprint (od)
'''
#for 1 sequence
seqs = [iLOV]
codon_count = CodonsDict
for i in seqs:
    dna_sequence = i[0][1]
    #print (dna_sequence)
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        print (codon)
        if codon in codon_count:
            codon_count[codon] += 1
od = collections.OrderedDict(sorted(codon_count.items()))
pprint.pprint (od)'''
