#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import warnings
import copy
import my

class MGIMousePhenotypeFilesError(Exception): pass

#-H /scratch1/tmp/myourshaw/mouse_phenoypes/HMD_HumanPhenotype.rpt -V /scratch1/tmp/myourshaw/mouse_phenoypes/VOC_MammalianPhenotype.rmwhitespace.rpt

name = 'pypeline.MGI_mouse_phenotype_files'
version = 1
copyright = '©2013 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

def run(HMD, VOC):
    
    hmd_header = ('#Human_Gene', 'Human_Entrez_Gene_ID', 'HomoloGene_ID', 'Mouse_Marker_Symbol', 'MGI_Marker_ID', 'Mammalian_Phenotype_ID')
    line_count = 0
    data_input_line_count = 0
    data_output_line_count = 0
    print 'processing ' + HMD
    fixed_rpt = HMD+'.flat'
    with open(HMD) as rpt, open(fixed_rpt, 'w') as fixed:
        fixed.write('\t'.join(hmd_header)+'\n')
        for line in rpt:
            line_count += 1
            line = line.rstrip('\n')
            if not bool(line.strip()) or line.startswith('#'):
                continue
            else: #data line
                data_input_line_count += 1
                fields = [f.strip() for f in line.split('\t')]
                phenotypes = [p.strip() for p in fields[5].split()]
                if not phenotypes:
                    phenotypes = ['']
                for phenotype in phenotypes:
                    output_line = '\t'.join(fields[:5]) + '\t' + phenotype + '\n'
                    fixed.write(output_line)
                    data_output_line_count += 1
    print 'read {} input lines, {} data input lines\nwrote {} data output lines to {}'.format(line_count, data_input_line_count, data_output_line_count, fixed_rpt)

    voc_header = ('#Mammalian_Phenotype_ID', 'Mammalian_Phenotype_Name', 'Description')
    line_count = 0
    data_input_line_count = 0
    data_output_line_count = 0
    print 'processing ' + VOC
    fixed_rpt = VOC+'.clean'
    with open(VOC) as rpt, open(fixed_rpt, 'w') as fixed:
        fixed.write('\t'.join(voc_header)+'\n')
        for line in rpt:
            line_count += 1
            line = line.rstrip('\n')
            if not bool(line.strip()) or line.startswith('#'):
                continue
            else: #data line
                data_input_line_count += 1
                fields = [f.strip() for f in line.split('\t')]
                output_line = '\t'.join(fields) + '\n'
                fixed.write(output_line)
                data_output_line_count += 1
    print 'read {} input lines, {} data input lines\nwrote {} data output lines to {}'.format(line_count, data_input_line_count, data_output_line_count, fixed_rpt)

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'flatten HMD_HumanPhenotype.rpt file and cleanup VOC_MammalianPhenotype.rpt file from MGI-Mouse Genome Informatics',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--HMD', '-H', required = True,
        help='input VOC_MammalianPhenotype.rpt file')
    parser.add_argument('--VOC', '-V', required = True,
        help='input VOC_MammalianPhenotype.rpt file')
    args = parser.parse_args()

    run(HMD=args.HMD, VOC=args.VOC)

if __name__ == "__main__": sys.exit(main())
