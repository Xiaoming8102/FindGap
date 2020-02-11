#!/usr/bin/env python
import os
import gzip
import sys
import re

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)

def read_fasta (fasta_file) :
    """
        Iterates over all sequences in a fasta file. One at a time,
        without reading the whole file into the main memory.
    """

    try :
        INPUT = (gzip.open if fasta_file.endswith('.gz') else open)(fasta_file)
    except IOError:
        print "[Error] Cannot find Fasta file : %s !" % fasta_file
        exit(-1)
    sanitize = re.compile(r'[^ACTGN]')
    sanitize_seq_id = re.compile(r'[^A-Za-z0-9]')

    chrome_seq = []
    chrome_id = None
    seen_ids = set() 

    for line in INPUT :
        if line[0] == '>':
            if chrome_id is not None:
                yield chrome_id, ''.join(chrome_seq)

            chrome_id = sanitize_seq_id.sub('_', line.split()[0][1:])
            if chrome_id in seen_ids:
                error('Found identical sequence ids (id: %s) in the fasta file: %s.'
                      ' Please, make sure that all sequence ids are unique and '
                      ' contain only alphanumeric characters: A-Za-z0-9_' % (chrome_id, fasta_file))
            seen_ids.add(chrome_id)
            chrome_seq = []
        else:
            chrome_seq.append(sanitize.sub('N', line.strip().upper()))
    yield chrome_id, ''.join(chrome_seq)
    INPUT.close()


def Findgap ( FastaFn, outputFn ):
    # the input fasta filename appended

    cut_format = "NNNN"

    total_chr = 0
    len_chr = dict()
    OUT = open(outputFn, 'w') if outputFn else sys.stdout

    for chr, seq in read_fasta(FastaFn):
        total_chr += 1
        L = len(seq)
        len_chr[chr] = L

        start_pos = 0
        counter = 0
        gap=False
        gap_length = 0
        for i in range(0,L):
            if seq[i] == "N":
                if gap_length == 0:
                    start_pos=counter
                    gap_length = 1
                    gap = True
                else:
                    gap_length += 1
            else:
                if gap:
                    OUT.write("%s\n" % "\t".join( [chr, str(start_pos), str(start_pos+gap_length), str(gap_length)] ))
                    gap_length = 0
                    gap = False
            counter += 1

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog -i <genome.fa> [-o <output>]\n" \
            "Author : Guo, Weilong Xie Xiaoming; guoweilong@gmail.com 2016301010312@cau.edu.cn; 2020-02-11\n" \
            "Last Update: 2020-02-11\n" \
            "Description: Get the positions of all the NNNNNN(gap) fragments\n"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Genome sequence file in Fasta format", metavar="FILE")
    parser.add_option("-o", dest="outfile",
                  help="Name of the output file (standard output if not specified). Format: chr Nstart_pos Nend_pos (0-base)\n", metavar="FILE")
    (options, args) = parser.parse_args()

    if (options.infile is None) :
        print parser.print_help()
        exit(-1)
    Findgap (options.infile, options.outfile)


# ===========================================
if __name__ == "__main__":
    main()
