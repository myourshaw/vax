#!/bin/bash
 
#ensembl_database_install_server.sh
#runs on mysql server
#assumes mysql files downloaded and unzipped
USAGE="usage: $0 <hgmd directory> <hgmd version> <host> <port> <user> [<password>]"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi
 
#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");
 
HERE=${script_dir};
wd=`pwd`;
 
if (( $# < 5 )) || [[ ! -d ${1} ]]; then echo ${USAGE}; exit; fi

HGMD_DIR=${1}; shift;
VERSION=${1}; shift;
HOST=${1}; shift;
PORT=${1}; shift;
USER=${1}; shift;
if [ -z "${1}" ]; then
	read -p "Enter MySQL password for database installation (or run using the -p option): " -s -r
	if [ -z ${REPLY} ]
	then
		echo;
		log "Installer aborted."
		exit 1;
	else
		PASSWORD=${REPLY};
	fi
else
  PASSWORD=${1}; shift;
fi

MYSQL="$(which mysql) -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";
MYSQLIMPORT="mysqlimport -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";

#drop/recreate databases and import data
#hgmd_views should be last
for DB in {hgmd_pro,hgmd_snp,hgmd_phenbase,hgmd_views}; do
  echo "dropping and recreating ${DB}";
  ${MYSQL} -e "DROP DATABASE IF EXISTS ${DB};";
  ${MYSQL} -e "CREATE DATABASE ${DB};";
  echo "importing ${DB}";
  cat ${HGMD_DIR}/${DB}-${VERSION}.dump | ${MYSQL} ${DB}
done

#add indices
echo "indexing hgmd_pro.hg19_coords";
${MYSQL} -e "USE hgmd_pro;
ALTER TABLE hgmd_pro.hg19_coords
ADD INDEX IX_chromosome (chromosome ASC)
, ADD INDEX IX_coordSTART (coordSTART ASC)
, ADD INDEX IX_coordEND (coordEND ASC);";

#stored procedures
echo "creating hgmd_pro.coord2hgmd stored procedure";
${MYSQL} -e "USE hgmd_pro;
DROP procedure IF EXISTS hgmd_pro.coord2hgmd;
DELIMITER \$\$
CREATE DEFINER=CURRENT_USER PROCEDURE coord2hgmd(chromosome VARCHAR(2), coordSTART INT(11), coordEND INT(11))
BEGIN
select distinct
allmut.acc_num,
allmut.disease,
allmut.tag,
allmut.base,
allmut.hgvs,
allmut.codon,
allmut.amino,
allmut.deletion,
allmut.insertion,
allmut.descr
FROM hgmd_pro.hg19_coords
JOIN hgmd_pro.allmut
ON hg19_coords.acc_num = allmut.acc_num
where hg19_coords.chromosome = chromosome
and (hg19_coords.coordSTART BETWEEN  coordSTART and coordEND
or hg19_coords.coordEND BETWEEN  coordSTART and coordEND);
END\$\$
DELIMITER ;";

echo "creating hgmd_pro.gene2hgmd_disease stored procedure";
${MYSQL} -e "USE hgmd_pro;
DROP procedure IF EXISTS gene2hgmd_disease;
DELIMITER \$\$
CREATE DEFINER=CURRENT_USER PROCEDURE gene2hgmd_disease(gene VARCHAR(10))
BEGIN
SELECT distinct disease
from hgmd_pro.allgenes
where allgenes.gene = gene;
END\$\$
DELIMITER ;";

#permissions
echo "setting permissions";
for DB in {hgmd_pro,hgmd_views,hgmd_snp,hgmd_phenbase}; do
  ${MYSQL} -e "GRANT SELECT, INSERT, UPDATE, CREATE TEMPORARY TABLES ON ${DB}."'*'" TO 'hgmduser'@'%' IDENTIFIED BY 'mootation';";
  ${MYSQL} -e "GRANT SELECT, EXECUTE ON ${DB}."'*'" TO 'vax'@'%' IDENTIFIED BY 'vax';";
done

#test
echo "testing hgmd_pro.gene2hgmd_disease";
$(which mysql) -h ${HOST} --user=vax -P ${PORT} -pvax -e "CALL hgmd_pro.gene2hgmd_disease('EXOSC3');";
echo "testing hgmd_pro.coord2hgmd";
$(which mysql) -h ${HOST} --user=vax -P ${PORT} -pvax -e "CALL hgmd_pro.coord2hgmd('9', 37780792, 37785040);";


