#!/bin/bash

#download_dbnsfp_data.sh
USAGE="usage: $0 </top/level/output/directory> <version> </path/to/install/directory"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};

#output directory
if (( $# < 3 )); then echo USAGE; exit 1; fi

DOWNLOAD_DIR=$1;
VERSION=$2;
DBNSFP_DIR=$3;

mkdir -p $DBSNFP_DIR;

wget -N http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv${VERSION}.zip -P ${DOWNLOAD_DIR}
unzip ${DOWNLOAD_DIR}/dbNSFPv${VERSION}.zip
cat ${DOWNLOAD_DIR}/dbNSFP${VERSION}_variant.chr* | bgzip -c > ${DBNSFP_DIR}/dbNSFP.gz
tabix -f -s 1 -b 2 -e 2 ${DBNSFP_DIR}/dbNSFP.gz
