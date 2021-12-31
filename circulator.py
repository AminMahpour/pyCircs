#!/usr/bin/env python3
import re
from os.path import exists
import logging
import multiprocessing
from multiprocessing import Pool
from itertools import groupby
import argparse
version = "0.1a"
logo = """


       .__                     .__                       __   
  ____ |__|______   ____  __ __|  | _____ _______  _____/  |_ 
_/ ___\|  \_  __ \_/ ___\|  |  \  | \__  \\_  __ \/  _ \   __\\
\  \___|  ||  | \/\  \___|  |  /  |__/ __ \|  | \(  <_> )  |  
 \___  >__||__|    \___  >____/|____(____  /__|   \____/|__|  
     \/                \/                \/                   
Translate potential CDSs in circRNAs.\nBy Amin Mahpour.\n"""
examples = """
Examples:
Circulator.py input.fa output.fa #Simple usage
Circulator.py -t 4 input.fa output.fa #Add 4 threads
Circulator.py -t 4 -l input.fa output.fa #Add 4 threads and keep the longest CDS.
Circulator.py -t 4 -l -m 25 input.fa output.fa #Add 4 threads and keep the longest CDS and set a cutoff of 25 amino acids.
if you used this application in your project please cite:
XXX.XXXX 
\n
"""

def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """
    "first open the file outside "
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield headerStr, seq.upper().replace("U", "T")


parser = argparse.ArgumentParser(prog="Circulator.py", formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=logo, epilog=examples)
parser.add_argument("fasta_file", metavar="Input", help="Input FASTA file.")
parser.add_argument("output_file", metavar="Output", help="Output protein FASTA file.")
parser.add_argument("-m", "--min", type=int, default=20,
                    help='Cutoff CDS size. default: 20.')

parser.add_argument("-t", "--threads", type=int, default=0,
                    help='Requested processor numbers. Default(0) uses all available cores. example: -t 4')

parser.add_argument("-l", "--longest", action="store_true", default=False, help="Only keep longest ORFs.")
parser.add_argument("-r", "--rna", action="store_true", default=False, help="Input fasta file is RNA sequence.")

parser.add_argument("-v", "--version", action="version", version=version)

args = parser.parse_args()
logging.basicConfig(level=logging.INFO)
# print(args.accumulate(args.integers))
logging.info(rf"input file: {args.fasta_file}")
logging.info(rf"output file: {args.output_file}")



input_file = args.fasta_file
output_file = args.output_file

if not (exists(input_file)):
    logging.critical("Input file not found.")
    exit(1)
keep_longest = args.longest
min_aa = args.min
fastaRNA = args.rna


if args.threads == 0:
    NUM_PROCS = multiprocessing.cpu_count()
else:
    NUM_PROCS = args.threads

len_limit = 0

# NUM_PROCS = multiprocessing.cpu_count()


def translate(seq):
    table = {
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
    protein: str = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


def cdsSeeker(seq: str) -> [int, int]:
    # print("Starting...")
    startList: [re.Match] = []
    endList: [re.Match] = []
    orfList = set()

    for starts in re.finditer("ATG", string=seq):
        # print(starts.start())
        startList.append(starts)

    for ends in re.finditer("TAA|TAG|TGA", string=seq):
        # print(ends.start())
        endList.append(ends)

    for i in startList:
        for j in endList:
            if abs(i.start() - j.end()) % 3 == 0 and \
                    abs(i.start() - j.end()) != 0 and \
                    i.start() < j.start():
                orfList.add((i.start(), j.end()))
                # print("---")
                # print(i.start(), j.end())
                # print(seq[i.start(): j.end()])
                # print("---")

    return orfList


def findCircleCDS(seq: [str]) -> [str]:
    name: str = seq[0]
    duplicateSequence: str = seq[1] * 2
    CDSsequence = set()

    for i in cdsSeeker(duplicateSequence):
        if abs(i[0] - i[1]) >= min_aa * 3:
            CDSsequence.add(duplicateSequence[i[0]:i[1]])

    translated = set()
    for i in CDSsequence:
        thisCDS = translate(i).strip("_")
        if thisCDS.find("_") == -1:
            if len(thisCDS) >= len_limit:
                translated.add((name, i, translate(i).strip("_")))
                # print(name,i,translate(i).strip("_"))

    final_list = list(dict.fromkeys(translated))
    # print (final_list)
    if not keep_longest:
        return final_list
    else:

        res = max(final_list, key=lambda item: len(item[2]), default=0)
        if res == 0:
            return []

        return [res]


if __name__ == '__main__':
    # example = "TAGTTTTTTTTTTATGGGAGGA"
    # cdses = findCircleCDS(example)
    # print(cdses)
    print(rf"Running with {NUM_PROCS} CPUs...")
    fastaIterator = fasta_iter(input_file)
    with Pool(processes=NUM_PROCS) as pool:
        x = pool.map(findCircleCDS, fastaIterator)

    cleaned = [i for i in x if i != []]

    f = open(file=output_file, mode="w")

    for i in cleaned:
        for j in i:
            f.write(rf">{j[0]}")
            f.write("\n")
            f.write(j[2])
            f.write("\n")

    f.close()
#   for x in fasta_iter(input_file):
#   print(x[0])
