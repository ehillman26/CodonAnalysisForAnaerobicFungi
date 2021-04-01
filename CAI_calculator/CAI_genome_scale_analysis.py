# -*- coding: utf-8 -*-
"""
Created on May 8 2020

@author: Ethan Hillman
"""
from __future__ import print_function

import os
import re
import sys
import pprint
import collections
from collections import OrderedDict
import csv
import math
import copy
from pprint import pprint as pp
from pathlib import Path

from Bio import SeqIO
from Bio.File import as_handle
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter
from Bio.SeqIO.Interfaces import _clean, _get_seq_string

from CAI import CAI, relative_adaptiveness


### -- PATH --- ###
cwd = os.getcwd()

codon_count = {}
CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
         'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0,
        #'TGA': 0,'TAA': 0, 'TAG': 0, 
              }

SynonymousCodons = { 
 'CYS': ['TGT', 'TGC'], 
 'ASP': ['GAT', 'GAC'], 
 'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
 'GLN': ['CAA', 'CAG'], 
 'MET': ['ATG'], 
 'ASN': ['AAC', 'AAT'], 
 'PRO': ['CCT', 'CCG', 'CCA', 'CCC'], 
 'LYS': ['AAG', 'AAA'], 
# 'STOP': ['TAG', 'TGA', 'TAA'], 
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

#from .CodonUsageIndices import SharpEcoliIndex
SharpEcoliIndex = {
    "GCA": 0.586, "GCC": 0.122, "GCG": 0.424, "GCT": 1, "AGA": 0.004,
    "AGG": 0.002, "CGA": 0.004, "CGC": 0.356, "CGG": 0.004, "CGT": 1, "AAC": 1,
    "AAT": 0.051, "GAC": 1, "GAT": 0.434, "TGC": 1, "TGT": 0.5, "CAA": 0.124,
    "CAG": 1, "GAA": 1, "GAG": 0.259, "GGA": 0.01, "GGC": 0.724, "GGG": 0.019,
    "GGT": 1, "CAC": 1, "CAT": 0.291, "ATA": 0.003, "ATC": 1, "ATT": 0.185,
    "CTA": 0.007, "CTC": 0.037, "CTG": 1, "CTT": 0.042, "TTA": 0.02,
    "TTG": 0.02, "AAA": 1, "AAG": 0.253, "ATG": 1, "TTC": 1, "TTT": 0.296,
    "CCA": 0.135, "CCC": 0.012, "CCG": 1, "CCT": 0.07, "AGC": 0.41,
    "AGT": 0.085, "TCA": 0.077, "TCC": 0.744, "TCG": 0.017, "TCT": 1,
    "ACA": 0.076, "ACC": 1, "ACG": 0.099, "ACT": 0.965, "TGG": 1, "TAC": 1,
    "TAT": 0.239, "GTA": 0.495, "GTC": 0.066, "GTG": 0.221, "GTT": 1}

RARE_CODONS = {'CTA': 0,'ATA': 0,'CCC': 0,'TGT': 0,'CGA': 0,'CGG': 0,'AGA': 0,'AGG': 0}
SemiRARE_CODONS = {'TTA': 0,'CTT': 0,'CTC': 0,'TCC': 0,'TCA': 0,'TCG': 0,'AGT': 0,'CCT': 0, 'CCA': 0}

def cai_for_gene(dna_sequence):
    """Calculate the CAI (float) for the provided DNA sequence (string).
    This method uses the Index (either the one you set or the one you 
    generated) and returns the CAI for the DNA sequence. 
    """
    cai_value, cai_length = 0, 0 
    
    # if no index is set or generated, the default SharpEcoliIndex will 
    # be used.
    index = SharpEcoliIndex
    codon_count = copy.deepcopy(CodonsDict)
    gene_index = copy.deepcopy(CodonsDict)
    
    for i in range(0, len(dna_sequence), 3):
        if len(dna_sequence) > 299:
            #print (codon_count)
            codon = str(dna_sequence[i:i + 3])
            if codon in index:
                if codon in codon_count:
                    codon_count[codon] += 1
            od = collections.OrderedDict(sorted(codon_count.items()))
        else:
            None

    #print ("done")
    if len(dna_sequence) > 299:
        pprint.pprint (od)
    else:
        None

    for aa in SynonymousCodons:
        if len(dna_sequence) > 299:
            total = 0.0 
            # RSCU - relative synonymous codon usage
            rscu = []
            codons = SynonymousCodons[aa]
            
            for item in codons:
                total += od[item]
                rscu.append(od[item])
                #print (item)
                #print (od[item])
                #print (total)

            # calculate the RSCU value for each of the codons
            rscu_max = max(rscu) 
            for item in codons:
                if rscu_max > 0:
                    w = od[item] / rscu_max
                    gene_index[item] = w
                    #denominator = float(total) / len(codons)
                    #rcsu.append(od[item] / denominator)
                
            # now generate the index W=RCSUi/RCSUmax: 
            #rcsu_max = max(rcsu) 
            #for codon_index, codon in enumerate(codons): 
            #    gene_index[codon] = rcsu[codon_index] / rcsu_max
        else:
            None

    #pprint.pprint (gene_index)
    #pprint.pprint (od)

    #for codon in gene_index:
        #if len(dna_sequence) > 299:

            #print (codon)
            #print (index[codon])
            #print (gene_index[codon])
        
'''     if codon not in ["ATG", "TGG"]:
cai_value += math.log(index[codon])
cai_length += 1
            # some indices may not include stop codons:
            elif codon not in ["TGA", "TAA", "TAG"]:
                raise TypeError("illegal codon in sequence: %s.\n%s"
                                % (codon, index))
                
    return math.exp(cai_value / (cai_length - 1.0))
'''



fasta_file = 'C:/Users/ethil/OneDrive/Documents/Mevalonate/E.coli_CodonUsage/E_coli_CDS.txt'
ref_seqs = 'C:/Users/ethil/OneDrive/Documents/Mevalonate/E.coli_CodonUsage/e_coli-HEG.fasta'
records = list(SeqIO.parse(fasta_file, "fasta"))

sequences = [seq.seq for seq in SeqIO.parse(ref_seqs, "fasta")]
weights = relative_adaptiveness(sequences=sequences)

'''ref_seqs = []
for item in refs:
    if((len(item.seq))%3 != 0):
        print (item.id, "not divisible by 3")
    else:
        ref_seqs.append(str(item.seq))
        #print (len(item.seq))
        #print (item.id)

#print (len(ref_seqs))
#print (ref_seqs)
#mystring =  ','.join(ref_seqs)
#print (mystring)
'''
CAI_GENE_Info = []
for item in records:
    if((len(item.seq))%3 != 0):
        print (item.id, "not divisible by 3")
    else:
        #print (str(item.seq))
        cai = CAI(str(item.seq), weights=weights)#ref_seqs[0],ref_seqs[1],ref_seqs[2],ref_seqs[3],ref_seqs[4],ref_seqs[5],ref_seqs[6],ref_seqs[7],ref_seqs[8],ref_seqs[9],ref_seqs[10],ref_seqs[11],ref_seqs[12],ref_seqs[13],ref_seqs[14],ref_seqs[15],ref_seqs[16],ref_seqs[17],ref_seqs[18],ref_seqs[19],ref_seqs[20],ref_seqs[21],ref_seqs[22],ref_seqs[23],ref_seqs[24],ref_seqs[25],ref_seqs[26],ref_seqs[27],ref_seqs[28],ref_seqs[29],ref_seqs[30],ref_seqs[31],ref_seqs[32],ref_seqs[33],ref_seqs[34],ref_seqs[35],ref_seqs[36],ref_seqs[37],ref_seqs[38],ref_seqs[39]))
        print (item.id, cai)
        #Gene_info = [item.id, cai_for_gene(item.seq)]
        #print (str(item.seq[0:3]))
        #CAI_GENE_Info.append(Gene_info)
        #gene_cai = cai_for_gene(item.seq)
        #print (gene_cai)
        #print (CAI_GENE_Info[0])
