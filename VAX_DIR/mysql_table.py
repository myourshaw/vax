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

class SqlColumnsError(Exception): pass

#dbsnp135
#-i /scratch1/tmp/myourshaw/resources/dbsnp135/00-All.vcf

#-i /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulationNoGeno/YRI-1412-Y-nogeno.vcf.gz
#-i /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulationNoGeno/ByPopulationNoGeno.variants.txt
#-i /scratch1/tmp/myourshaw/resources/dbsnp135/00-All.vcf.flat.txt
#--skip ## --header_line #CHROM -i /scratch0/tmp/myourshaw/1000genomes/phase1_integrated_release_version3_20120430/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#-i /scratch0/tmp/myourshaw/hgnc/hgnc_xref_20120902.txt --header_line 1
#-i /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.snps_indels.vcf_allele_counts.txt.gz
#-i /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.snps_indels.vcf_genotype_counts.txt.gz
#-i /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.chrY.snps_indels.vcf
#-i /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.snps_indels.vcf_allele_counts.txt.gz /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.snps_indels.vcf_genotype_counts.txt.gz
#-i /scratch0/tmp/myourshaw/gene_intervals/get_genomic_regions_ALLGENES.txt.genes.txt /scratch0/tmp/myourshaw/gene_intervals/get_genomic_regions_ALLGENES.txt.genomic_region_details.txt /scratch0/tmp/myourshaw/gene_intervals/get_genomic_regions_ALLGENES.txt.genomic_target_intervals.txt
#--column_names chrom chromStart chromEnd name score strand -i /scratch0/tmp/myourshaw/gene_intervals_test2/ensembl69_known_putative_novel_pseudo_non-coding_CDS_ess_5sr4_3sr13_5utr_3utr.bed
#--no_column_names -i /scratch0/tmp/myourshaw/gene_intervals_test2/ensembl69_known_putative_novel_pseudo_non-coding_CDS_ess_5sr4_3sr13_5utr_3utr.bed
#-i /data/storage-1-01/archive/myourshaw/ensembl/gene_intervals/known_putative_novel_pseudo_non-coding_CDS_ess_5sr4_3sr13_5utr_3utr.genes.txt /data/storage-1-01/archive/myourshaw/ensembl/known_putative_novel_pseudo_non-coding_CDS_ess_5sr4_3sr13_5utr_3utr.genomic_region_details.txt
#-i /scratch1/tmp/myourshaw/mmjj_20130514/vcfs/mmjj_20130514.hc.raw.indels.vcf


def run(input, table=None, skip=None, header_line='#', column_names=None, no_column_names=False):
    
    header_line_number = int(header_line) if my.is_int(header_line) and int(header_line) > 0 else None
    header_chars = header_line if not my.is_int(header_line) else None
    
    output = os.path.join(os.path.dirname(input),table+'.mysql') if table else input+'.table.mysql'
    with my.open_gz_or_text(input) as input_file, open(output, 'w') as mysql:
        column_widths = {}
        column_names = {}
        column_names_list = None
        line_count = 0
        first_data_line = None
        header_found = False if header_line_number or header_chars else True
        if column_names:
            header_found = True
            column_names_list = column_names
            column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
            sql_column_spec = my.get_sql_column_spec(column_names_list)
        for line in input_file:
            line_count += 1
            line = line.rstrip('\n')
            if not bool(line.strip()):
                continue
            if not header_found:
                if no_column_names:
                    header_found = True
                    column_names_list = ['col'+str(i+1) for i in range(len(line.split('\t')))]
                    column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
                    sql_column_spec = my.get_sql_column_spec(column_names_list)
                elif (header_line_number and header_line_number == line_count) or (header_chars and header_chars != '#' and line.startswith(header_chars)) or (header_chars and header_chars == '#' and line.startswith(header_chars) and line[1] != '#'):
                    header_found = True
                    column_names_list = line.split('\t')
                    column_names = {i:column_names_list[i] for i in range(len(column_names_list))}
                    sql_column_spec = my.get_sql_column_spec(column_names_list)
                continue
            elif skip and line.startswith(skip):
                continue
            else: #data line
                if first_data_line == None:
                    first_data_line = line_count
                fields = line.split('\t')
                indices = [i for i in range(len(fields))]
                column_widths = {i: max(len(fields[i]),column_widths.setdefault(i,0)) for i in indices}
                sql_data_dict = my.get_sql_data_dict(column_names, fields)
                my.update_sql_column_spec(sql_column_spec, sql_data_dict)
                pass
        #output.write('col_number\tcol_name\tcol_width\n')
        #output.write(''.join(['{}\t{}\t{}\n'.format(col+1, column_names[col], column_widths[col]) for col in sorted(column_widths)]))

        if not (column_names_list and sql_column_spec):
            raise SqlColumnsError('Could not identify columns. Try using --column_names or --no_column_names.')

        mysql_scripts = my.get_mysql_scripts(
            table_name=table if table else os.path.basename(input).replace('.','_')[:64],
            index_base_name=os.path.basename(input).replace('.','_')[:16],
            indices=[],
            columns_out=column_names_list,
            columns_out_spec=sql_column_spec,
            rows_to_delete=first_data_line-1
            )
        mysql.write(mysql_scripts[0])

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'calculate maximum column widths for tab-delimited file',
        epilog = 'pypeline.column_widths version 1.0β1 ©2011-2013 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required = True,
        help='input  file(s)[.gz]')
    parser.add_argument('--table', '-t',
        help='input  file[.gz]')
    parser.add_argument('--skip',
        help='skip lines that start with character(s)')
    parser.add_argument('--header_line', default = '#',
        help='1-based row number of column names or # if the first row that starts with a single # character')
    parser.add_argument('--column_names', nargs='*',
        help='list of column names (when file does not contain column name header)')
    parser.add_argument('--no_column_names', action='store_true', default=False,
        help='file does not contain column name header, so use "col1, col2, ..."')
    args = parser.parse_args()

    if not my.file_exists(args.input):
        raise SqlColumnsError('The file '+args.input+' does not exist.')
    if args.table and len(args.table) > 64:
        raise SqlColumnsError('Table name '+args.table+' is too long (max 64).')
    run(input=args.input,table=args.table,skip=args.skip,header_line=args.header_line,column_names=args.column_names,no_column_names=args.no_column_names,)

if __name__ == "__main__": sys.exit(main())
