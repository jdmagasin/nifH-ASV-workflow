#!/usr/bin/bash

OUTDIR=SraMetadata
##
## From the SRA, download metadata for all runs from each of the projects
## specified, one project accession per line in the input file.
##
## Usage:
##    downloadSraMetadata.sh <projects.txt>
##
## Input file: Comment lines in the input file begin with '#'.  Blank lines are
## ignored.
##
## Output: in directory $OUTDIR, one file per project.  Because the project
##         accession is one of the fields, you can simply concatenate all the
##         outfiles if you like.
##
## Requires NCBI Entrez Programming Utilities (E-utilities). [Jonathan, you have
## these in your NCBI_tools conda environment.]
##

## Exit immediately upon any error
set -e

PROJFILE=$1
if [ ! -f "$PROJFILE" ] ; then
    echo "Need a file listing the project accessions.  See docs at top of script."
    exit -1
fi
if [ -z `which esearch` ] ; then
    echo "I cannot find esearch (from NCBI E-Utilities)."
    exit -1
fi
if [ -d "$OUTDIR" ] ; then
    echo "$OUTDIR already exists. Abort!"
    exit -1
fi
mkdir "$OUTDIR"

## - Drop comments, blank lines. Also remove dups even though this will change
##   the order in which projects are downloaded.
## - Prevent esearch from hijacking stdin with echo trick described here:
##     https://stackoverflow.com/questions/346445/bash-while-read-loop-breaking-early
## - To be safe from NCBI ignoring too many/frequent requests, sleep.
while read -r ACC; do
    OFILE="${OUTDIR}/${ACC}.csv"
    echo "Getting on project $ACC"
    echo "" | esearch -db sra -query "$ACC" | efetch -format runinfo > "$OFILE"
    sleep 1s
done < <(cat "$PROJFILE" | egrep -v '(^$|^#)' | sort | uniq)
echo "Done!"
exit 0


## Simple test file with these two accessions should get 71 total records:
## PRJNA476143
## PRJNA111397
