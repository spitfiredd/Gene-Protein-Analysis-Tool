"""provides analysis on DNA sequences"""
from itertools import zip_longest


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


class DnaString:
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    contents = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    rna = {'A': 'A', 'T': 'U', 'G': 'G', 'C': 'C'}
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    @classmethod
    def update(cls, dnastring):
        for nucleotide in dnastring:
            cls.contents[nucleotide] += 1

    def __init__(self, dnastring):
        self.dnastring = dnastring
        self.update(dnastring)

    def complement(self):
        st = ''
        for nucleotide in self.dnastring:
            st += self.complement_dict[nucleotide]
        return st

    def reverse_complement(self):
        return self.complement()[::-1]

    def rna_translation(self):
        st = ''
        for nucleotide in self.dnastring:
            st += self.rna[nucleotide]
        return st

    def cg_pct(self):
        total = sum(self.contents.values())
        gc = 0
        for key, val in self.contents.items():
            if key in ('GC'):
                gc += val
        return gc / total


def hamming_differences(DNA1, DNA2):
    """Returns info on positional/nucleotide differences in 2 DNA strands in
    the form of a list"""
    hDiff = []

    if len(DNA1) > len(DNA2):
        for i in range(len(DNA1)):
            DNA2 += 'X'
    else:
        for i in range(len(DNA2)):
            DNA1 += 'X'
    for n1, n2 in list(zip(DNA1, DNA2)):
        if n1 != n2:
            locDiff = list(zip(DNA1, DNA2)).index((n1, n2))
            nucDiff = n1 + n2
            print(str(locDiff) + nucDiff)
            hDiff.append(str(locDiff) + nucDiff)
    return hDiff
