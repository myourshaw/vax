#!/bin/bash

#download_hpa_data.sh
USAGE="usage: $0 [/top/level/output/directory]"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};
wd=`pwd`;

#output directory
if (( $# > 0 )); then d=$1; else d=${wd}; fi

HPA_DIR=$d;
mkdir -p $HPA_DIR;

wget -N http://www.proteinatlas.org/download/normal_tissue.csv.zip  -P ${HPA_DIR};
unzip -o ${HPA_DIR}/normal_tissue.csv.zip -d ${HPA_DIR};
wget -N http://www.proteinatlas.org/download/subcellular_location.csv.zip  -P ${HPA_DIR};
unzip -o ${HPA_DIR}/subcellular_location.csv.zip -d ${HPA_DIR};
