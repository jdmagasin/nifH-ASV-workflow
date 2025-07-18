#!/usr/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

## Check whether each sample in the AUID abundance table (column name) can be
## mapped to a SAMPLEID in the metadata table.  Recall that AUID table columns
## have tag info, then "___", then a sample name as understood/used by the DADA2
## pipeline.  This script greps for the "sample name" (everything after the
## "___") in the metadata SAMPLEIDs.  AUID table "sample names" will usually
## fail to have found for them a SAMPLEID due to:
##  -- Inconsistency in how samples were named as seen by the pipeline
##     (e.g. based on the FASTQ names) versus in metadata tables.
##  -- X was prepended by R to numeric sample identifiers somewhere during the
##     pipeline but are missing in the metadata SAMPLEIDs.
##  -- Metadata SAMPLEIDs were assumed (and required )to be in the first column
##     of the metadata files but they were something else (
##  -- Metatranscriptomic samples have _transcriptomic appended to their sample
##     IDs.  Drop the suffix before trying to find the sample.  (Note that there
##     could be a metagenomic sample with the same ID but without the suffix.)
##

cat("Checking whether sample names (columns) in the AUID abundance table can be found in the\n",
    "metadata SAMPLEIDs.  Will match everything after the \"___\" in the abundance column sample\n",
    "names (i.e., will exclude the tags that precede \"___\") to the metadata SAMPLEIDs.\n")

mSampIds <- read.table('metadata.csv', sep=',', header=T)$SAMPLEID
stopifnot(!(table(mSampIds) > 1))
mSampIds <- unique(sub('_transcriptomic$', '', mSampIds))  # Use original sample ID for metaT's

aSampIds <- colnames(read.table('../FilterAuids/auid.abundances.filtered.tsv.gz', nrows=3, check.names=F))
stopifnot(table(aSampIds) == 1)  # uniqueness check
x <- table( sub('^.+___(.+)$','\\1', aSampIds) )
if (any(x > 1)) {
    cat("Some column names in the AUID abundance table have identical",
        "strings at their ends, after the '___':\n")
    print(names(which(x > 1)))
}
aSampIds <- names(which(x == 1))  # makes aSampIds be just the part after ___
cat("\n\n")

## Which of the AUID tab samples can be identified in metadata.
cat("Searching...")

## Switched from grep() to match() b/c some sample names have ".", but if we use
## grep(..., fixed=T) then a match of the sample name *within* the meta SAMPLEID
## will pass.  Probably don't want to allow that. So match() seems the best
## option.
a2m <- match(aSampIds, mSampIds)  # Exact match. (Takes the first.)
names(a2m) <- aSampIds
cat("\n")
if (any(table(a2m)) > 1) {
    ## This doesn't happen.**  If it does, report the sample names.
    warning("Some sample names seem to be duplicated in the AUID table, because\n",
            "they map to the same metadata SAMPLEID.\n")
}
cat(sum(!is.na(a2m)), "sample names in the AUID table map (identically) to metadata SAMPLEIDs.\n")
cat(sum(is.na(a2m)),  "sample names in the AUID table do NOT map to a metadata SAMPLEID.\n")
if (sum(is.na(a2m))) {
    cat("The unmapped sample names are\n")
    print(names(a2m)[is.na(a2m)])
    cat("\n\nPossible reasons the sample names above cannot be found in metadata SAMPLEIDs:\n",
        "\t- Different naming convention used for the sample names that were seen by\n",
        "\t  the DADA2 pipeline vs. what is in metadata tables.  The former depends on\n",
        "\t  how you set up the input to organizeFastqs.R. (The final column defined the\n",
        "\t  sample.) The latter is probably from the SRA run number.\n\n",
        "\t- The metadata table with information for the sample(s) is missing, possibly\n",
        "\t  dropped by gatherMetadata.R because it is misformatted.\n\n",
        "\t- The sample name was numeric so R prepended an \"X\" to it during the pipeline.\n\n")
}
cat("Done checking mapping from AUID abundance table to metadata table.\n")
