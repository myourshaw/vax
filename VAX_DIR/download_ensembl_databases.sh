#!/bin/bash
#download_ensembl_databases.sh
USAGE="usage: $0 Ensembl_version_number [/top/level/output/directory]"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -f '%Y' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};

wd=`pwd`;

#Ensembl version, e.g., 72
if (( $# < 1 )); then echo ${USAGE}; exit; fi
v=$1;
#output directory
if (( $# > 1 )); then d=$2; else d=${wd}; fi
local_directory=${d}/release-${v}/mysql;
mkdir -p ${local_directory};

#download via rsync
cmd='rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_mysql/homo_sapiens_*'" ${local_directory}";
${cmd};

cmd='rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_mysql/ensembl_*'" ${local_directory}";
${cmd};

#unzip
for d in ${local_directory}/*; do
cmd="gunzip -rf ${d}/";
${cmd}; 
done

cd ${wd};
