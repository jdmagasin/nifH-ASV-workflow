#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

##
## Creates one big abundance table, asv2auid.abundances.tsv, which has ASV counts
## for all the studies.  Inputs are the AUID abundance tables that were created
## for the studies by makeAUIDCountTable.R.
##

ofile <- 'asv2auid.abundances.tsv'
cat("Making sure",ofile,"does not exist already.\n")
stopifnot(!file.exists(ofile))
cat("Okay, will aggregate AUID abundance tables into one big table...\n")

##
## Load abundance tables. Then copy them into one massive matrix, with
## rows=AUIDs and columns=samples.
##
abundFiles <- list.files('CountTables', pattern='^counts\\.*', full=T)
tags <- sub('^counts\\.(.+)\\.tsv','\\1', basename(abundFiles))
abundTabs <- lapply(1:length(abundFiles),
                    function(i) {
                        ## Get counts, and prefix samples with the tag.
                        m <- as.matrix(read.table(abundFiles[i], sep="\t"))
                        colnames(m) <- paste0(tags[i],"___",colnames(m))  # see hack note below
                        m
                    })
## Make sure AUIDs and sample names are not duplicated.
auids <- unlist(lapply(abundTabs, rownames));  x <- unique(auids);  stopifnot(setequal(x,auids));  auids <- x
samps <- unlist(lapply(abundTabs, colnames));  x <- unique(samps);  stopifnot(setequal(x,samps));  samps <- x
stopifnot(setequal(sub('___.*','',samps),tags))  # ~Hack: rely on '___' not being in any tag
abunds <- matrix(0, nrow=length(auids), ncol=length(samps), dimnames=list(auids,samps))
sumReadsByTab <- 0  # sanity check
sumSampsByTab <- 0
for (i in 1:length(abundTabs)) {
    a <- rownames(abundTabs[[i]])
    s <- colnames(abundTabs[[i]])
    sumReadsByTab <- sumReadsByTab + sum(abundTabs[[i]])
    sumSampsByTab <- sumSampsByTab + ncol(abundTabs[[i]])
    abunds[a,s] <- abundTabs[[i]][a,s]
}
stopifnot(sumReadsByTab == sum(abunds))  # Whole should *equal* sum of parts:)
## Next check really should not be needed b/c when 'samps' is created we check
## for duplicates.  But let's do it anyway.
stopifnot(sumSampsByTab == ncol(abunds))

## Development debug: Check that the total counts for each AUID match what was
## calculated by prepareFastaForUchime3_denonovo.R and stuffed in the deflines.
if (FALSE) {
    library(ShortRead)
    VerifyAuidTots <- function() {
        x <- sub('size=','', as.character(id(readFasta('all_AUIDs_sizes.len280to340nt.fasta'))))
        x <- strsplit(x,';')
        auid2size <- sapply(x, function(v) as.numeric(v[1]))
        names(auid2size) <- sapply(x, function(v) v[2])
        stopifnot( rowSums(abunds)[names(auid2size)] == auid2size )
    }
    VerifyAuidTots()
}

write.table(abunds, ofile, sep="\t", quote=F)
cat("Wrote",ofile,"which has all AUIDs (rows) and columns of form\n",
    "     tag___sample\n",
    "where tag encodes metadata for a DADA2 run (e.g. the region and size fraction).\n\n")

cat("Number of AUIDs   (rows) =",nrow(abunds),"\n")
cat("Number of samples (cols) =",ncol(abunds),"\n")
cat("Total reads              =",sum(abunds),"\n")