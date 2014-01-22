#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import csv
import datetime
from time import strftime
import re

class Omim2DbError(Exception): pass
#mysql import
#WORKDIR="/Volumes/scratch/OMIM";
#table=omim_genemap;
#COLLIST=`head -1 $WORKDIR/$table.txt`
#mysql -u $USER -p$PASSWORD --execute "TRUNCATE TABLE vw.$table";
#mysql -u $USER  -p$PASSWORD --execute "LOAD DATA LOCAL INFILE '$WORKDIR/$table.txt' INTO TABLE vw.$table IGNORE 1 LINES; SHOW WARNINGS;" > $table.output;

#-d /scratch1/tmp/myourshaw/resources/omim/20111230
def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'Convert ftp://anonymous:myourshaw%40ucla.edu@grcf.jhmi.edu/OMIM/ to tab-delimited database files',
        epilog = 'pypeline.omim2db version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--dir', '-d', required=True,
        help='directory of downloaded OMIM files')
    args = parser.parse_args()
    genemap_in = os.path.join(args.dir, 'genemap')
    morbidmap_in = os.path.join(args.dir, 'morbidmap')
    mim2gene = os.path.join(args.dir, 'mim2gene.txt')
    if not (os.path.exists(args.dir) and os.path.exists(genemap_in)  and os.path.exists(morbidmap_in) and os.path.exists(mim2gene)):
        raise  Omim2DbError('{} must contain genemap, morbidmap, and mim2gene.txt')
    genemap_in_fields = [
        'Chromosome_Map_Entry_Number',
        'Month_entered',
        'Day_entered',
        'Year_entered',
        'Cytogenetic_location',
        'Gene_Symbols',
        'Gene_Status',
        'Title',
        'Title2',
        'MIM',
        'Method',
        'Comments',
        'Comments2',
        'Disorders',
        'Disorders2',
        'Disorders3',
        'Mouse_correlate',
        'Reference'
        ]
    genemap_out_fields = [
        '#chromosome',
        'map_entry',
        'date_entered',
        'cytogenetic_location',
        'gene',
        'gene_status',
        'title',
        'mim_locus',
        'method',
        'comments',
        'disorder',
        'mim_phenotype',
        'phene_mapping_key',
        'mouse_correlate',
        'reference'
        ]
    morbidmap_in_fields = (
        'Disorder',
        'symbols',
        'MIM',
        'cytogenetic_location'
        )
    morbidmap_out_fields = (
        '#disorder',
        'mim_phenotype',
        'phene_mapping_key',
        'symbols',
        'mim_locus',
        'cytogenetic_location'
        )
    disorder_re = re.compile(r'(?P<Disorder>.+?), (?P<MIM_phenotype>\d+) \((?P<phene_mapping_key>\d+)\)', re.I)
    disorder2_re = re.compile(r'(?P<Disorder>.+?) \((?P<phene_mapping_key>\d+)\)', re.I)
    genemap_out = os.path.join(args.dir, 'omim_genemap.txt')
    gm_out = csv.DictWriter(open(genemap_out, "w"), genemap_out_fields, delimiter='\t')
    gm_out.writeheader()
    with open(genemap_in, 'rb') as f:
        gm_in = csv.DictReader(f, fieldnames=genemap_in_fields, delimiter='|')
        for row in gm_in:
            if row['Gene_Symbols'] == 'HTC2, HCG, CGH, CXINSq27.1,':
                foo = row['Gene_Symbols']
                bar = re.split(r', *', row['Gene_Symbols'])
                baz = [g for g in re.split(r', *', row['Gene_Symbols']) if g]
            genes = [g for g in re.split(r', *', row['Gene_Symbols']) if g]
            disorders = ' '.join([row['Disorders'].strip(), row['Disorders2'].strip(), row['Disorders3'].strip()]).strip()
            disorders = re.split(r'; *', disorders)
            for gene in genes:
                for disorder in disorders:
                    gm = {}
                    gm['#chromosome'] = row['Chromosome_Map_Entry_Number'].split('.')[0]
                    gm['map_entry'] = row['Chromosome_Map_Entry_Number'].split('.')[1]
                    #gm['Date_entered'] = strftime("%Y-%m-%d", datetime.date.timetuple(datetime.date(int(('19' if int(row['Year_entered']) >= 60 else '20')+row['Year_entered']), int(row['Month_entered']), int(row['Day_entered']))))
                    gm['date_entered'] = '-'.join([row['Year_entered'],row['Month_entered'],row['Day_entered']])
                    gm['cytogenetic_location'] = row['Cytogenetic_location']
                    gm['gene'] = gene
                    gm['gene_status'] = row['Gene_Status']
                    gm['title'] = ' '.join([row['Title'], row['Title2']]).strip()
                    gm['mim_locus'] = row['MIM']
                    gm['method'] = ','.join(re.split(r', *', row['Method']))
                    gm['comments'] = ' '.join([row['Comments'], row['Comments2']]).strip()
                    d_re = disorder_re.match(disorder)
                    if d_re != None:
                        gm['disorder'] = d_re.group('Disorder')
                        gm['mim_phenotype'] = d_re.group('MIM_phenotype')
                        gm['phene_mapping_key'] = d_re.group('phene_mapping_key')
                    else:
                        d_re = disorder2_re.match(disorder)
                        if d_re != None:
                            gm['disorder'] = d_re.group('Disorder')
                            gm['phene_mapping_key'] = d_re.group('phene_mapping_key')
                        else:
                            gm['disorder'] = disorder
                    gm['disorder'] = gm['disorder'].lstrip('{').rstrip('}').lstrip('[').rstrip(']')
                    gm['mouse_correlate'] = row['Mouse_correlate']
                    gm['reference'] = row['Reference']
                    gm_out.writerow(gm)
    morbidmap_out = os.path.join(args.dir, 'omim_morbidmap.txt')
    mm_out = csv.DictWriter(open(morbidmap_out, "w"), morbidmap_out_fields, delimiter='\t')
    mm_out.writeheader()
    with open(morbidmap_in, 'rb') as f:
        mm_in = csv.DictReader(f, fieldnames=morbidmap_in_fields, delimiter='|')
        for row in mm_in:
            mm ={}
            disorder = disorder_re.match(row['Disorder'])
            if disorder != None:
                mm['#disorder'] = disorder.group('Disorder')
                mm['mim_phenotype'] = disorder.group('MIM_phenotype')
                mm['phene_mapping_key'] = disorder.group('phene_mapping_key')
            else:
                disorder = disorder2_re.match(row['Disorder'])
                if disorder != None:
                    mm['#disorder'] = disorder.group('Disorder')
                    mm['phene_mapping_key'] = disorder.group('phene_mapping_key')
                else:
                    mm['#disorder'] = row['Disorder']
            mm['#disorder'] = mm['#disorder'].lstrip('{').rstrip('}').lstrip('[').rstrip(']')
            mm['symbols'] = ','.join(re.split(r', *', row['symbols']))
            mm['mim_locus'] = row['MIM']
            mm['cytogenetic_location'] = row['cytogenetic_location']
            mm_out.writerow(mm)

if __name__ == "__main__": sys.exit(main())
