#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import re
import my

#http://biostar.stackexchange.com/questions/2145/bulk-download-of-ncbi-gene-summary-field

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'Convert RefGene Genbank protein files to tab-delimited database files',
        epilog = 'pypeline.refgene2db version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--ref', '-r', required=True,
        help='path to gene_RefSeqGene file')
    parser.add_argument('--input', '-i', nargs='+',
        help='downloaded refseqgeneN.genomic.gbff.gz files')
    parser.add_argument('--output', '-o', required=True,
        help='output file')
    args = parser.parse_args()
    
    locus2gene = {}
    with open(args.ref) as ref:
        for r in ref:
            if r.startswith('#'):
                continue
            r = r.rstrip('\n').split('\t')
            locus2gene[r[3].split('.')[0]] = r[2]
    
    inputs = my.unglob(args.input)
    locus2comment = {}
    locus_count = 0
    refseq_re = re.compile(r"\s*\[provided by RefSeq.*\]",re.I)
    for input in inputs:
        f = my.open_gz_or_text(input)
        in_comment=False
        for line in f:
            if line[0:5] == "LOCUS":
                locus = line.split()[1]
                comment = ""
                locus_count += 1
            elif line[0:7] == "COMMENT":
                in_comment=True
                comment += line.split("    ")[1].replace("\n", " ")
            elif line[0:7] == "PRIMARY":
                in_comment = False
                try:
                    summary = comment.split("Summary:")[1]#.strip().split('[provided by RefSeq].')[0].rstrip()
                except:
                    summary = comment#.strip().split('[provided by RefSeq].')[0].rstrip()
                locus2comment[locus] = refseq_re.split(summary)[0]
            elif in_comment:
                comment += line.split("            ")[1].replace("\n", " ")
    with open(args.output,'w') as output:
        output.write('#NG_ID\tGeneSymbol\tSUMMARY\n')
        for locus in sorted(locus2comment):
            output.write('{}\t{}\t{}\n'.format(locus, locus2gene.get(locus,''), locus2comment[locus]))

if __name__ == "__main__": sys.exit(main())
