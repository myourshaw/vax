-- MySQL
-- Table structure for table `HumanMitoCartaAll`
-- Note that this table contains records for all human genes; use MITOCARTA_LIST=1 for the 1023 MitoCarta genes.
--
DELIMITER $$
DROP TABLE IF EXISTS `HumanMitoCartaAll`$$
CREATE TABLE HumanMitoCartaAll (
	TRAINING_DATASET varchar(50),
	HUMAN_ENTREZ int(10),
	MOUSE_ENTREZ int(10),
	MAPPING varchar(50),
	SYM varchar(50),
	SYNONYMS varchar(255),
	DESCRIPTION varchar(255),
	TARGETP int(5),
	YEASTMITOHOMOLOG int(5),
	RICKETTSIA_BLASTP_EXPECT float,
	MITODOMAIN int(5),
	PGC_INDUCTION float,
	COEXPRESSION_GNF_N50 int(10),
	PROTEOMICS varchar(50),
	MAESTRO_SCORE float,
	MITOCARTA_cFDR float,
	MITOCARTA_LIST int(1),
	MCARTA_TYPE varchar(50),
	NUM_MSMS_TISSUES int(10),
	CHROMOSOME varchar(50),
	START int(15),
	STOP int(15),
	PFAM_DOMAINS varchar(255),	
  PRIMARY KEY (`SYM`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1$$