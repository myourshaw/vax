#!/bin/bash

#download_refseq_data.sh
USAGE="usage: $0 [/top/level/output/directory]"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};

#output directory
if (( $# > 0 )); then d=$1; else d=${wd}; fi

REFSEQ_DIR=${d};
mkdir -p $REFSEQ_DIR;

cmd="ncftpget ftp.ncbi.nih.gov ${REFSEQ_DIR} "'/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene';
${cmd};
cmd="ncftpget -R ftp.ncbi.nih.gov ${REFSEQ_DIR} "'/refseq/H_sapiens/RefSeqGene/refseqgene*.genomic.gbff.gz';
${cmd};

