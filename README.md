# pancanmafmerge

## Dependencies

A repo of code used to merge [TCGA PancanAtlas VCFs](https://wiki.nci.nih.gov/x/2gcYAw) from various somatic callers, into a deduplicated merged multisample [MAF format](https://wiki.nci.nih.gov/x/eJaPAQ). This repository borrows from code available in other repositories and has some custom code to handle the merger. The project's main purpose is simply to do the merge for this one project and therefore has no features to make it more generic or maintainable beyond this goal.

Members of the TCGA Network can get more information about the PancanAtlas effort [at this link](https://wiki.nci.nih.gov/display/TCGAM/PanCancerAtlas).


```sql
.schema                              
CREATE TABLE dbsnpvalstat (rsid, valstatus); 
CREATE INDEX rsid on dbsnpvalstat (rsid);    
```

- Ensembl's variant effect predictor (VEP) v82 installed for offline use with the grch37 cache [Instructions](https://gist.github.com/ckandoth/9d6ad6a7fd3b058e5bc98a1ce884641a)
- COSMIC annotations (http://cancer.sanger.ac.uk/cosmic) specifically the COSMIC mutations VCF file.We are unable to distribute this with the project and users will have to get this file on their own.  We use a filtered version of this file with only the COSMIC variants with 3 or greater entries in COSMIC.
- dbSNP validation status database built in sqlite3 and the associated sqlite3 jdbc driver (sqlite-jdbc-3.8.11.2/sqlite-jdbc-3.8.11.2.jar). You may use dbsnp_valstatus.py to generate the initial validation status maping .txt file and then load this into an sqlite3 database. You should then index the rsid column. The schema for our dbsnp rsid database is:

You should also have the required interpreters and virtual machines installed:
perl, python, java, jython

## Running

### merge.py

merge.py is the main application of this project and is used to merge VCF files from a variety of formats into a single VCF file and then converts that into a MAF file.  The following is an example of the command line for running merge.py

```bash
python scripts/merge.py --vcfs radia.vcf muse.vcf --callers RADIA MUSE --tmpdir tmp merged.vcf
```

merge.py performs the following opperations (in this order):

1. filter: registered filter subroutines are called depending on the indicated caller
2. sort: the VCF files are sorted based on the order of contigs in the associated .fa.dict file for the reference
3. vcf2vcf: all vcfs are converted to a "lowest common denominator" vcf using vcf2vcf
4. merge: all vcfs are merged using `vcf-merge.py` to create a single vcf.  Note that the vcf created by `vcf-merge.py` has a very strange header and should only be used within this workflow unless the user excercises extreme caution.
5. annotation: the merged file is first annotated with VEP and then with a post-VEP canannotation utility that adds COSMIC, context, and dbSNP information to the merged VCF
6. vcf2maf: the merged and annotated vcf is converted using vcf2maf (with some modifications, see the fork covingto/vcf2maf)

### dispatch_server.py

Understandably, formatting, monitorying, and error revovery with such a large project is essential.  It is envisioned that this system may be run to merge together very large collections of VCF files across different infrastructures.  To that end, `dispatch_server.py` helps to organize and run the large job sets for this project.

Three main functions exist within `dispatch_server.py`; dispatcher, worker, and queue.

#### dispatcher

To begin you must set up manifest files containing a primary key (same across all manifests) and file paths to locate the vcf files for your project (relative to some base directory).  Then, you must set up a configuration file which is used by dispatcher to load a job list.  An outline of the configuration file is:

```json
{
    "fmaps": [
        {
            "basepath": "/path/to/base/directory",
            "fmap": "/path/to/manifest/file",
            "pkey": "primary_key_column_name",
            "mapping": {"pindel": "PINDEL", "varscan.somatic": "VARSCANS"}
        }
    ]
}
```

Once you have loaded the dispatcher, you are ready to start workers.

#### worker

The worker should be given the host and port of the dispatcher and be started in an environment with about 2 CPUs and 16gb of available RAM and have access to the file system that contains the vcf data and output directories.

Once started, the workers will communicate with the dispatcher to take jobs.  If an error is detected in the execution of a merge, the worker will communicate back to the dispatcher that there is a non-zero exit code and that job will be moved to "problems", which can be queried by sending {"request": "problems"} to the dispatcher (host:port) and parsing the returned json.  An email is also sent to kylecovington1 at gmail dot com, this is currently hard coded into the `merge.py` program.  Please change this so that kylecovington1 at gmail dot com is not saturated with emails should you use this application.

Workers stop working once there are no more jobs.

## Reporting bugs

What?! You found bugs?

<img src="http://i.giphy.com/ph6ewybUlGbW8.gif" width="240">

That's impossible! Please double check and triple check, and if you're super certain that we're not infallible, then [click here](https://github.com/covingto/pancanmafmerge/issues) to report the bug, along with a clear explanation of how we can find or reproduce it.
