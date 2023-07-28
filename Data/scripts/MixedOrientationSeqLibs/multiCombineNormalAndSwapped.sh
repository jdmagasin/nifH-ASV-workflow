#!/usr/bin/bash

set -e

echo
echo "Will combine Normal and Swapped pipeline outputs for each of the processing groups."
echo

## The directory of this script also must have combineAsvsFromNormalAndSwappedRuns.R
SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
if [ ! -x "$SDIR/combineAsvsFromNormalAndSwappedRuns.R" ] ; then
    echo "Cannot find combineAsvsFromNormalAndSwappedRuns.R"
    exit -1
fi

## Loop over processing groups
## The "sort" is required and corresponds to sort(list.files()) in the R script.
## It ensure that $i here and indexOfAsvTables in the R script index into the
## same list.
i=0
for pg in `find Normal -name asvs.noChimera.tsv | sort`; do
    i=$((i+1))
    echo
    if [ `echo "$pg" | grep -c 'DNA'` -gt 0 ] ; then desc=DNA ; else desc=RNA ; fi
    echo "#${i} run ~~~~~~~[ $desc ASVs ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo
    $SDIR/combineAsvsFromNormalAndSwappedRuns.R $i
    mv "asvs.noChimera.combined.$i.fasta" "asvs.noChimera.combined.$desc.fasta"
    mv "asvs.noChimera.combined.$i.tsv"   "asvs.noChimera.combined.$desc.tsv"
    mv "cluster_normal_vs_swapped.$i.png" "cluster_normal_vs_swapped.$desc.png"
    echo
    echo "Renamed outputs to use \"$desc\" rather than index $i"
    echo
done

echo "Done!"
exit 0


