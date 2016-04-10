# pancanmafmerge

A repo of code used to merge [TCGA PancanAtlas VCFs](https://wiki.nci.nih.gov/x/2gcYAw) from various somatic callers, into a deduplicated merged multisample [MAF format](https://wiki.nci.nih.gov/x/eJaPAQ). This repository borrows from code available in other repositories and has some custom code to handle the merger. The project's main purpose is simply to do the merge for this one project and therefore has no features to make it more generic or maintainable beyond this goal.

Members of the TCGA Network can get more information about the PancanAtlas effort [at this link](https://wiki.nci.nih.gov/display/TCGAM/PanCancerAtlas).

### Dependencies

- Ensembl's variant effect predictor (VEP) v82 installed for offline use with the grch37 cache [Instructions](https://gist.github.com/ckandoth/9d6ad6a7fd3b058e5bc98a1ce884641a)
- COSMIC annotations (http://cancer.sanger.ac.uk/cosmic) specifically the COSMIC mutations VCF file.We are unable to distribute this with the project and users will have to get this file on their own.  We use a filtered version of this file with only the COSMIC variants with 3 or greater entries in COSMIC.
- dbSNP validation status database built in sqlite3 and the associated sqlite3 jdbc driver (sqlite-jdbc-3.8.11.2/sqlite-jdbc-3.8.11.2.jar). You may use dbsnp_valstatus.py to generate the initial validation status maping .txt file and then load this into an sqlite3 database. You should then index the rsid column. The schema for our dbsnp rsid database is:

```
    sqlite> .schema
    CREATE TABLE dbsnpvalstat (rsid, valstatus);
    CREATE INDEX rsid on dbsnpvalstat (rsid);
```

- A .fa fasta file of the human reference (grch37) and an associated .dict file for this reference.
- Annotation of the context for each variant also depends on the htsjdk .jar file and picard (picard-tools-1.129/picard.jar, picard-tools-1.129/picard-lib.jar, picard-tools-1.129/htsjdk-1.129.jar), be sure that you have these.
- You should also have the required interpreters and virtual machines installed: perl, python, java, jython

### Reporting bugs

What?! You found bugs?

<img src="http://i.giphy.com/ph6ewybUlGbW8.gif" width="240">

That's impossible! Please double check and triple check, and if you're super certain that we're not infallible, then [click here](https://github.com/covingto/pancanmafmerge/issues) to report the bug, along with a clear explanation of how we can find or reproduce it.
