# pancanmafmerge
A repo for code specifically tailored to the merger of the TCGA PanCan VCF (per-caller) files into a MAF file.  This repository borrows from code available in other repositories and has some custom code to handle the merger.  The project's main purpose is simply to do the merge for this one project and therefore has no features to make it more generic or maintainable beyond this goal.

Dependencies

This project depends on variant effect predictor (VEP) to be installed and functioning.  We are using VEP-82 with the grch37 version of the human reference.  This should be downloaded and cached as per the instructions in VEP.
This project depends on COSMIC annotations (http://cancer.sanger.ac.uk/cosmic) specifically the COSMIC mutations VCF file.  We are unable to distribute this with the project and users will have to get this file on their own.  We use a filtered version of this file with only the COSMIC variants with 3 or greater entries in COSMIC.
Thsi project depends on a dbSNP validation status database built in sqlite3 and the associated sqlite3 jdbc driver (sqlite-jdbc-3.8.11.2/sqlite-jdbc-3.8.11.2.jar).  You may use dbsnp_valstatus.py to generate the initial validation status maping .txt file and then load this into an sqlite3 database.  You should then index the rsid column.  The schema for our dbsnp rsid database is:

sqlite> .schema                              
CREATE TABLE dbsnpvalstat (rsid, valstatus); 
CREATE INDEX rsid on dbsnpvalstat (rsid);    

This project depends on a .fa file of the human reference (grch37) and an associated .dict file for this reference.  Annotation of the context for each variant also depends on the htsjdk .jar file and picard (picard-tools-1.129/picard.jar, picard-tools-1.129/picard-lib.jar, picard-tools-1.129/htsjdk-1.129.jar), be sure that you have these.

You should also have the required interpreters and virtual machines installed:
perl, python, java, jython
