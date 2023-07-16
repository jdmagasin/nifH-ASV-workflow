#!/usr/bin/env Rscript

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
##  -- Metatranscriptomic samples were dropped by gatherMetadata.R but could
##     have been run through the pipeline.
##

cat("Checking with match() whether AUID sample names (columns) can be looked up in metadata\n",
     "SAMPLEIDs.  For this script use everything after the \"___\" in the AUID column\n",
     "(i.e. drop the tags preceding the \"___\") as the sample name.\n\n")

mSampIds <- read.table('metadata.tsv', header=T)$SAMPLEID
stopifnot(!(table(mSampIds) > 1))

droppedSamps <- readLines('dropped_metatranscriptomic_samples.txt')

##aSampIds <- colnames(read.table('../GatherAsvs/asv2auid.abundances.tsv.gz', nrows=3))
aSampIds <- colnames(read.table('../FilterAuids/auid.abundances.filtered.tsv.gz', nrows=3))
stopifnot(table(aSampIds) == 1)  # uniqueness check
x <- table( sub('^.+___(.+)$','\\1', aSampIds) )
if (any(x > 1)) {
    cat("Some column names in the AUID abundance table have identical",
        "strings at their ends, after the '___':\n")
    print(names(which(x > 1)))
}
aSampIds <- names(which(x == 1))  # makes aSampIds be just the part after ___
cat("\n\n")

## Which of the AUID tab samples can be identified in metadata.  Ignore
## droppedSamps.
cat("Searching...")

## Switched from grep() to match() b/c some sample names have ".", but if we use
## grep(..., fixed=T) then a match of the sample name *within* the meta SAMPLEID
## will pass.  Probably don't want to allow that. So match() seems the best
## option.
keptSamps <- setdiff(aSampIds, droppedSamps)
a2m <- match(keptSamps, mSampIds)  # Exact match. (Takes the first.)
names(a2m) <- keptSamps
cat("\n")
if (any(table(a2m)) > 1) {
    ## This doesn't happen.  If it does, report the sample names.
    warning("Some sample names seem to be duplicated in the AUID table, because\n",
            "they map to the same metadata SAMPLEID.\n")
}
cat(sum(!is.na(a2m)), "sample names in the AUID table map (identically) to metadata SAMPLEIDs.\n")
cat(sum(is.na(a2m)),  "sample names in the AUID table do NOT map to a metadata SAMPLEID.\n")
if (sum(is.na(a2m))) {
    cat("The unmapped sample names are\n")
    print(names(a2m)[is.na(a2m)])
}
cat("Note that the checks just above *excluded* samples in dropped_metatranscriptomic_samples.txt.\n")

cat("\n\nPossible reasons the sample names above cannot be found in metadata SAMPLEIDs:\n",
    "\t- Different naming convention used for the sample names that were seen by\n",
    "\t  the DADA2 pipeline vs. what is in metadata tables.  The former depends on\n",
    "\t  how you set up the input to organizeFastqs.R. (The final column defined the\n",
    "\t  sample.) The latter is probably from the SRA run number.\n\n",
    "\t- The metadata table with information for the sample(s) is missing, possibly\n",
    "\t  dropped by gatherMetadata.R because it is misformatted.\n\n",
    "\t- The sample name was numeric so R prepended an \"X\" to it during the pipeline.\n\n")

cat("Done checking mapping from AUID abundance table to metadata table.\n")
