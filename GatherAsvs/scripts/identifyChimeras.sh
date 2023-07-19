#!/usr/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

## Identify chimeras using uchime3_denovo (via my post DADA2 pipeline script).
## Creates replacement ASV fasta files with just the non-chimeras.  Each fasta
## includes as part of its name the tag from the INPUTTABLE.

WDIR="NonChimericAsvs"
LOG="chimera.log"

## Exit when any command fails, ande give an error message on EXIT.  Also track
## the last command executed (DEBUG).
## 
set -e
trap 'echo -e "The following command failed with exit code $?:\n${last_command}\n\nCheck $WDIR/$LOG for clues.\n"' EXIT
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG


INPUTTABLE="asvs.noChimera.fasta_table.tsv"   # Input
if [ ! -f "$INPUTTABLE" ] ; then
    echo "Missing $INPUTTABLE"
    exit -1
fi

## Requires appropriate conda environment.
MYCHIMERASCRIPT=`which check_chimera_denovo.sh`
if [ ! -x "$MYCHIMERASCRIPT" ] ; then
    echo "Cannot find $MYCHIMERASCRIPT"
    exit -1
fi

echo "Identifying chimera using $MYCHIMERASCRIPT"
date
echo
mkdir -p "$WDIR"
cat "$INPUTTABLE" | egrep -v '(^#|^$)' | tr -s "\t" > "$WDIR/fastatab.tmp"
cd "$WDIR"
echo "" > "$LOG"
while read line ; do
    dsnam=`echo "$line" | cut -f1`
    dstag=`echo "$line" | cut -f3 | tr -d "\n"`
    dsnam="${dsnam}_tag_${dstag}"
    asvpath=`echo "$line" | cut -f2`
    if [ ! -f "$asvpath" ] ; then
        ## Maybe failed because path is relative and we moved to WDIR.  Try again.
        asvpath="../$asvpath"
        if [ ! -f "$asvpath" ] ; then
            echo "ERROR: $asvpath does not exist.  Check your $INPUTTABLE"
            exit -1
        fi
    fi
    if [ `cat "$asvpath" | grep -c '^>'` -eq 0 ] ; then
        echo "No sequences in $asvpath ?"
        exit -1
    fi
    ## Assume abundance table is alongside and similarly named as
    ## the FASTA
    tsvpath=`echo "$asvpath" | sed 's/\.fasta$/.tsv/'`

    if [ -f "$dsnam.chimera.out" ] ; then
        echo "Already have chimera results for $dsnam" >> "$LOG"
        continue
    fi
    
    echo -e "#######################################\n" >> "$LOG"
    echo "## Checking for chimeras in $dsnam" >> "$LOG"
    echo "## Checking for chimeras in $dsnam"
    command="$MYCHIMERASCRIPT $tsvpath $asvpath"
    eval "$command" &>> "$LOG"
    mv nonchimera.fasta "${dsnam}.noChimera.fasta"
    rm chimera.fasta
    mv Work_on_chimera "$dsnam.Work_on_chimera"
    echo >> "$LOG"
done < fastatab.tmp
rm fastatab.tmp

echo
echo "Finished!"
date
exit 0
