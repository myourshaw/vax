#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import csv
import argparse

#-i /Users/myourshaw/apps/vax/downloads/hpa/normal_tissue.csv

#command line arguments
parser = argparse.ArgumentParser(
    description = 'convert csv file to tab-delimited',
    epilog = 'pypeline.csv2tab version 1.0β1 ©2011-2013 Michael Yourshaw all rights reserved')
parser.add_argument('--input', '-i', required = True,
    help='input  file')
parser.add_argument('--output', '-o',
    help='output  file (default: <input>.tab)')
parser.add_argument('--delimiter', '-d', default=',',
    help='input field delimiter (default: ,')
parser.add_argument('--quotechar', '-q', default='"',
    help='input quote character (default: ")')
parser.add_argument('--output_delimiter', default = '\t',
    help='output field delimiter (default: ,')
parser.add_argument('--output_quotechar', default=None,
    help='input quote character (default: csv.QUOTE_NONE)')
args = parser.parse_args()

if not args.output:
    if args.output_delimiter == '\t':
        args.output = args.input + '.tab'
    elif args.output_delimiter == ',':
        args.output = args.input + '.csv'
    else:
        args.output = args.input + '.converted'
with open(args.input, 'rb') as csvfile, open(args.output, 'wb') as tabfile:
    reader = csv.reader(csvfile, delimiter=args.delimiter, quotechar=args.quotechar)
    writer = csv.writer(tabfile, delimiter=args.output_delimiter, quotechar=args.output_quotechar) if args.output_quotechar \
    else csv.writer(tabfile, delimiter=args.output_delimiter, quoting=csv.QUOTE_NONE)
    for row in reader:
        writer.writerow(row)
    
    
    