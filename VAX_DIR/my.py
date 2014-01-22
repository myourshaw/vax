#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import tempfile
from glob import glob
import gzip
from time import localtime, strftime
import math
import re
import collections
import argparse
from ConfigParser import SafeConfigParser
import csv
import pysam
import subprocess

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

class MyError(Exception): pass
class CheckFilesError(Exception): pass

def makedir(dir):
	try:
		os.makedirs(dir)
	except OSError:
		if os.path.isdir(dir):
			pass
		else:
			raise

def split_file(file, lines, suffix_length=8, prefix=None, repeat_header=False, compress=False):
    #returns list of output files
    if not file_exists(file):
        raise MyError("split_file can't find {}".format(file))
    if not is_number(lines) or lines < 1:
        raise MyError("split_file lines {} invalid".format(lines))
    if not is_number(suffix_length) or suffix_length < 1:
        raise MyError("split_file suffix_length {} invalid".format(suffix_length))

    #return value
    parts = []
    
    if not prefix:
        prefix = file
    file_suffix = prefix+'{:0='+str(suffix_length)+'}'+('.gz' if compress else '')
    lines_out = 0
    suffix = 0
    file_out = None
    header = ''
    if repeat_header:
        with open_gz_or_text(file) as f:
            for line in f:
                if line.startswith('#'):
                    header += line
                else:
                    break
    with open_gz_or_text(file) as f:
        for line in f:
            if not file_out or lines_out >= lines:
                lines_out = 0
                if file_out:
                    file_out.close()
                filename_out = file_suffix.format(suffix)
                suffix += 1
                parts += [filename_out]
                if compress:
                    file_out = gzip.open(filename_out, 'w')
                else:
                    file_out = open(filename_out, 'w')
                if repeat_header:
                    file_out.write(header)
            if repeat_header and line.startswith('#'):
                continue
            file_out.write(line)
            lines_out += 1
        
    if file_out:
        file_out.close()

    return parts

def get_immediate_subdirectory_paths(dir):
	return [os.path.join(dir, name) for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

def get_immediate_subdirectory_names(dir):
	return [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

def is_gz(filename):
	GZIP_MAGIC = b"\x1F\x8B"
	with open(filename, 'rb') as fh:
		magic = fh.read(len(GZIP_MAGIC))
	if magic == GZIP_MAGIC:
		return True
	else:
		return False

def open_gz_or_text(filename):
	if is_gz(filename):
		return gzip.open(filename, 'r')
	else:
		return open(filename, 'r')

def open_bam(filename):
	if is_gz(filename):
		return pysam.Samfile(filename, "rb")
	else:
		return pysam.Samfile(filename, "r")

def unique_dir(prefix='', dir=os.path.join(os.path.expanduser('~'),'tmp')):
	makedir(dir)
	return tempfile.mkdtemp(prefix=prefix+"_", dir=dir)

#http://stackoverflow.com/questions/183480/is-this-the-best-way-to-get-unique-version-of-filename-w-python
def unique_file(file_name):
	dirname, filename = os.path.split(file_name)
	prefix, suffix = os.path.splitext(filename)
	fd, filename = tempfile.mkstemp(suffix, prefix+"_", dirname)
	return os.fdopen(fd,'w'), filename

#deprecated
def list_files(pathnames):
    return sorted(set(flatten([glob(os.path.expanduser(os.path.expandvars(p))) for p in pathnames])))

def unglob(globs,sort_unique=True):
    if isinstance(globs, str):
        globs = [globs]
    if sort_unique:
        return sorted(set(flatten([glob(os.path.expanduser(os.path.expandvars(p))) for p in globs]))) if globs and isinstance(globs, (list, tuple, set)) else []
    else:
        return flatten([glob(os.path.expanduser(os.path.expandvars(p))) for p in globs]) if globs and isinstance(globs, (list, tuple, set)) else []

def get_files(names):
	return (file for file in names if os.path.isfile(file))

def is_int(s):
	try:
		n = int(s)
		if (math.isinf(n) or math.isnan(n)): return False
		return True
	except ValueError:
		return False

def is_number(s):
	try:
		n = float(s)
		if (math.isinf(n) or math.isnan(n)): return False
		return True
	except ValueError:
		return False

#http://code.activestate.com/recipes/363051/
#http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def file_exists(file_path):
	if isinstance(file_path, str) and os.path.isfile(file_path):
		try:
			open(file_path)
			return True
		except IOError as e: return False
	else: return False

def is_exe(file_path):
	return file_exists(file_path) and os.access(file_path, os.X_OK)

def list_existing_files(files):
	"""list_existing_files(sequence of paths) -> list of paths

	Returns a single, flat list of all paths that openable files
	or None if all files are missing"""

	result = []
	for f in flatten(files):
		if file_exists(f): result.append(f)
	return None if len(result) == 0 else result

def list_missing_files(files):
	"""list_missing_files(sequence of paths) -> list of paths

	Returns a single, flat list of all paths that are not files or cannot be opened
	or None if all files are present and accessible"""

	result = []
	for f in flatten(files):
		if not file_exists(f): result.append(f)
	return None if len(result) == 0 else result

def print_log(s):
	print '[{}] {}\n'.format(strftime("%Y-%m-%d %H:%M:%S", localtime()),s.rstrip('\n'))

def localtime_stamp():
    return '{}'.format(strftime("%Y-%m-%d %H:%M:%S", localtime()))
    
def localtime_squish():
	return '{}'.format(strftime("%Y%m%d%H%M%S", localtime()))

def print_err(s):
	sys.stderr.write('[{}] {}\n'.format(strftime("%Y-%m-%d %H:%M:%S", localtime()),s.rstrip('\n')))

def rcut(string, cut):
    return re.sub()

def bam_strip(string):
    s = string[:-3] if string.endswith('.gz') else string
    return s[:-4] if s.endswith('.bam') or s.endswith('.sam') else s

def qseq_strip(string):
    s = string[:-3] if string.endswith('.gz') else string
    return s[:-9] if s.endswith('_qseq.txt') else s

def vcf_strip(string):
    s = string[:-3] if string.endswith('.gz') else string
    return s[:-4] if s.endswith('.vcf') else s

def r_strip(string, strip):
    return string[:-len(strip)] if string.endswith(strip) else string

def get_config(args=None):
    config = SafeConfigParser()
    config.readfp(open(os.path.join(os.path.dirname(__file__),'pypeline.cfg')))
    cfg_files = [os.path.expanduser('~/.pypeline.cfg')]
    if args and 'config' in dir(args) and args.config:
        cfg_files.append(args.config)
    config.read(cfg_files)
    return config

#sql columns
def get_sql_column_spec(columns):
    return {c: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for c in columns}

def get_sql_data_dict(columns,values):
    return {columns[i]: values[i] for i in range(len(values))}

def update_sql_column_spec(sql_column_spec, sql_data_dict):
    for c in sql_data_dict.keys():
        this_value = sql_data_dict[c]
        spec = sql_column_spec.get(c)
        if spec and this_value != '':
            spec['size_min'] = len(this_value) if spec['size_min'] == None else min(len(this_value),spec['size_min'])
            spec['size_max'] = len(this_value) if spec['size_max'] == None else max(len(this_value),spec['size_max'])
            if is_number(this_value):
                spec['min'] = float(this_value) if spec['min'] == None else min(float(this_value),spec['min'])
                spec['max'] = float(this_value) if spec['max'] == None else max(float(this_value),spec['max'])
            if not spec['type'] == 'varchar':
                if not is_number(this_value):
                    spec['type'] = 'varchar'
                elif not spec['type'] == 'float':
                    if not is_int(this_value):
                        spec['type'] = 'float'
                    else:
                        spec['type'] = 'int'
#microsoft sql server
def get_microsoft_sql_server_scripts(
    table_name='[dbo].[table]',
    index_base_name='table',
    clustered_index=[],
    indices=[],
    columns_out=[],
    columns_out_spec={},
    rows_to_delete=0
    ):
    #returns [create_table_script, import_data_script, index_scripts]
    #create table script
    create_table_script = """SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO
--substitute your [schema].[table_name] for {0}
IF OBJECT_ID(N'{0}', N'U') IS NOT NULL 
DROP TABLE {0};
GO
CREATE TABLE {0} (
""".format(table_name)

    for c in columns_out:
        spec = columns_out_spec.get(c)
        if spec:
            if c == columns_out[0] and c.startswith('#'):
                c = c[1:]
            if spec['type'] == 'varchar':
                if spec['size_min'] == spec['size_max'] and spec['size_max']<=1000:
                    ms_type = '[CHAR] ({})'.format(spec['size_max'])
                else:
                    ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
            elif spec['type'] == 'float':
                if spec['min'] >= -3.4e+38 and spec['max'] <= 3.4e+38:
                    ms_type = '[REAL]'
                elif spec['min'] >= -1.79e+308 and spec['max'] <= 1.79e+308:
                    ms_type = '[FLOAT]'
                else:
                    ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
            elif spec['type'] == 'int':
                if spec['min'] >= 0 and spec['max'] <= 1:
                    ms_type = '[BIT]'
                elif spec['min'] >= 0 and spec['max'] <= 255:
                    ms_type = '[TINYINT]'
                elif spec['min'] >= -32768 and spec['max'] <= 32767:
                    ms_type = '[SMALLINT]'
                elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
                    ms_type = '[INT]'
                elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
                    ms_type = '[BIGINT]'
                else:
                    ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
            else:
                if c in ('compound_het'):
                    ms_type = '[BIT]'
                else:
                    ms_type = '[CHAR] (1)'
            ms_spec = '\t[{}] {} {}NULL{}\n'.format(c, ms_type, 'NOT ' if c in clustered_index else '', ',' if c != columns_out[-1] else '')
            create_table_script += ms_spec
    create_table_script += ') ON [PRIMARY]\nGO\n\nSET ANSI_PADDING OFF\nGO\n'
        
        #import data script
    import_data_script = """
--substitute your [schema].[table_name] for {0}
--substitute your path to unzipped file for {1}
--delete first {2} rows of comments before importing

DECLARE @bulk_cmd varchar(1000)
SET @bulk_cmd = 'BULK INSERT {0}
FROM ''{1}'' 
WITH (ROWTERMINATOR = '''+CHAR(10)+''',
FIRSTROW=2)'
EXEC(@bulk_cmd)
GO
""".format(table_name, r'C:\path\to\file\to\be\imported', rows_to_delete)

    #create index scripts
    index_scripts = ''
    #clustered_index is a list of columns
    if clustered_index:
        ix_name = 'IX_'+index_base_name+'_'+'_'.join(clustered_index)
        ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in clustered_index])
        index_scripts += """CREATE CLUSTERED INDEX [{0}] ON {1}
(
    {2}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(ix_name, table_name, ix_cols)
    #indices is a list of indices containing lists of columns
    for i in indices:
        ix_name = 'IX_'+index_base_name+'_'+'_'.join(i)
        ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in i])
        index_scripts += """CREATE NONCLUSTERED INDEX [{0}] ON {1}
(
    {2}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(ix_name, table_name, ix_cols)
    return [create_table_script, import_data_script, index_scripts]

def get_mysql_scripts(
    table_name='table',
    index_base_name='table',
    indices=[],
    columns_out=[],
    columns_out_spec={},
    rows_to_delete=0):
    #returns [create_table_script, import_data_script_unix, import_data_script_windows]
    
    create_table_script = \
"""delimiter $$
DROP TABLE IF EXISTS `{0}`$$
CREATE TABLE `{0}` (
""".format(table_name)  #table name length limit of 64 characters
    table_specs = []
    for c in columns_out:
        spec = columns_out_spec.get(c)
        if spec:
            if c == columns_out[0] and c.startswith('#'):
                c = c[1:]
            if spec['type'] == 'varchar':
                if spec['size_max']<=500:  #row size limit of 65535 bytes, therefore text preferred over varchar, choosing arbitrary limit
                    my_type = 'varchar({})'.format(spec['size_max'])
                else:
                    my_type = 'text'
            elif spec['type'] == 'float':
                my_type = 'float'
            elif spec['type'] == 'int':
                if spec['min'] >= 0 and spec['max'] <= 1:
                    my_type = 'tinyint(1)'
                elif spec['min'] >= 0 and spec['max'] <= 255:
                    my_type = 'tinyint unsigned'
                elif spec['min'] >= -32768 and spec['max'] <= 32767:
                    my_type = 'smallint'
                elif spec['min'] >= -8388608 and spec['max'] <= 8388607:
                    my_type = 'mediumint'
                elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
                    my_type = 'int'
                elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
                    my_type = 'bigint'
                elif spec['size_max']<=500:
                    my_type = 'varchar({})'.format(spec['size_max'])
                else:
                    my_type = 'text'
            else:
                if c in ('compound_het'):
                    my_type = 'tinyint(1)'
                else:
                    my_type = 'char(1)'
            my_spec = '\t`{}` {} NULL'.format(c.strip(), my_type)
            table_specs.append(my_spec)
    create_table_script += ',\n'.join(table_specs)
    if indices:
        keys = []       
        for i in indices:
            ix_name = 'IX_'+index_base_name+'_'+'_'.join(i)
            ix_cols = '`,`'.join(i)
        keys.append("""
KEY `{0}` (`{1}`)""".format(ix_name, ix_cols))
        create_table_script += ','.join(keys)
    create_table_script += """
) ENGINE=MyISAM DEFAULT CHARSET=latin1$$"""

    #mysqlimport
    import_data_script_unix = 'mysqlimport --password=password \\\n --delete --local --verbose --ignore-lines={} \\\nschema_name \\\n/path/to/{}/file.{}\n'.format(rows_to_delete, index_base_name, index_base_name)
    import_data_script_windows = '"C:\Program Files\MySQL\MySQL Server 5.5\bin\mysqlimport" --ignore-lines={} schema_name "C:\/path/to/{}/file.{}" -L -u root -ppassword'.format(rows_to_delete, index_base_name, index_base_name)

    return [create_table_script, import_data_script_unix, import_data_script_windows]
