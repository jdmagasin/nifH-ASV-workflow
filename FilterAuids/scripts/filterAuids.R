#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

usageStr = "
    filterAuids auidAbundTable.tsv{.gz}  auids.fasta  minLen maxLen  {contams}  {nifH-like dir}

  The two {optional} parameters are from the post-DADA2-pipeline checkers:
  
   - contams is the tabular BLAST output from
     Contaminants/check_nifH_contaminants.sh, e.g. something like blast.96id.out
     which has 96%id alignments to known contaminants.
     
   - nifH-like dir is the directory with results from
     CheckNifHLike/classifyNifH.sh. The directory includes the outputs
     positives.ids, negatives.ids, unsure.ids, and nohits.ids.

  Main output of this script is the filtered abundance table
  auid.abundances.filtered.tsv.gz
"

## Name of output, the filtered abundance table
FOUT = "auid.abundances.filtered.tsv.gz"

args <- commandArgs(T)
auidAbundTabName <- args[1]
auidFastaName    <- args[2]
lenFiltRange     <- as.numeric(args[3:4])
stopifnot(lenFiltRange[1] < lenFiltRange[2])
contamBlastOut   <- args[5] # optional
nifHlikeCheckDir <- args[6] # optional (#2 must be present to use this)


options(width = 120) # More columns for printing.


##------------------------------------------------------------------------------
##
## Parameters that affect filtering of AUIDs and samples
##

MIN_SAMPLE_SIZE = 500            # Drop samples with fewer than this many reads
AUID_MIN_SAMPS  = 2              # AUIDs must be in at least this many samples and have
AUID_MIN_READS  = 1              # at least this many reads in each of those samples
AUID_ONE_HIT_WONDER_READS = 1000 # Also keep AUIDs with >= this many reads >=1 samples

RETAIN_CONTAMINANTS = T  # If TRUE, do not drop contaminants. Just report them.
CONTAM_PCTID = 96        # AUIDs that align (nt) to known contaminants at this %id
CONTAM_ALEN  = 70        # and length will considered contaminants.

RETAIN_NON_NIFH_LIKE = F # If TRUE, do not drop non-NifH-like AUIDs.


##------------------------------------------------------------------------------
##
## Helper functions
##

dim2desc <- function(d) { paste(d[1],"AUIDs by",d[2],"samples") }

AnnounceSection <- function(title)
{
    cat("\n\n")
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    cat("##",title,"\n\n")
}


## Record changes in total reads and num ASVs.
gFilterStages <- data.frame()
RecordDeltasAtFilterStep <- function(stage,
                                     readsAfter, readsBefore,
                                     asvsAfter,  asvsBefore)
{
    df <- data.frame(Stage=stage,
                     ReadsAfter=readsAfter, ReadsBefore=readsBefore,
                     ReadsPctDrop=round(100*(readsBefore-readsAfter)/readsBefore, 1),
                     AsvsAfter=asvsAfter, AsvsBefore=asvsBefore,
                     AsvsPctDrop=round(100*(asvsBefore-asvsAfter)/asvsBefore, 1))
    gFilterStages <<- rbind(gFilterStages, df)  # Note the global assignment operator.
}

## Record number of reads in each sample at the end of the specified 'stage'.  Pass an
## abundance table (AUID rows x Sample cols)
gFilterStages.samples <- data.frame()
RecordSampleReadsAfterFilterStep <- function(stage, atab)
{
    df <- data.frame(colnames(atab), colSums(atab))
    colnames(df) <- c('Sample', stage)
    stopifnot(!(stage %in% colnames(gFilterStages.samples)))
    if (nrow(gFilterStages.samples) > 0) {
        ## One stage of FilterAuids can only drop samples as it progresses, so all.x=T, and
        ## make the counts be 0 for the dropped samples after the stage at which dropped.
        df <- merge(gFilterStages.samples, df, by='Sample', all.x = T)
        df[is.na(df)] <- 0
    }
    gFilterStages.samples <<- df
}


##------------------------------------------------------------------------------
##
## Load stuff
##
AnnounceSection("Setting up")
cat("Loading libraries...")
suppressMessages(library(ShortRead))  # Needed only for aseqs
cat("done.\n")

cat("Loading AUID abundance table...")
stopifnot(file.exists(auidAbundTabName))
atab <- read.table(auidAbundTabName, check.names=F) # do not check/fix sample names
cat("done!  Table is",dim2desc(dim(atab)),"\n")
RecordSampleReadsAfterFilterStep('Initial', atab)

cat("Loading FASTA...")
aseqs <- readFasta(auidFastaName)
id2len <- width(sread(aseqs))
names(id2len) <- sub(' +.*$','',as.character(id(aseqs)))
stopifnot(rownames(atab) %in% names(id2len))
cat("done!  Length range is",paste(range(id2len),collapse=" to "),"nt.\n")
rm(aseqs) # currenly only need lengths

contamAsvs <- NULL
if (!is.na(contamBlastOut)) {
    cat("Loading results from contaminants check. ")
    stopifnot(file.exists(contamBlastOut))
    if (file.size(contamBlastOut) == 0) {
        cat("There are no contaminants:  ", contamBlastOut, "is empty.\n")
    } else {
        contamAsvs <- unique(subset(contamAsvs, V3 >= CONTAM_PCTID & V4 >= CONTAM_ALEN)[,'V1'])
        cat(length(contamAsvs), "AUIDs align at >=", CONTAM_PCTID, "%id and >=",
            CONTAM_ALEN,"nt to contaminants.\n")
        if (length(contamAsvs) == 0) { contamAsvs <-  NULL }
    }
}
nifHlike <- NULL
if (!is.na(nifHlikeCheckDir)) {
    cat("Loading results from NifH-like classifications: ")
    stopifnot(dir.exists(nifHlikeCheckDir))
    nifHlike <- lapply(c('positives','negatives','unsure'), function(x) {
        readLines(file.path(nifHlikeCheckDir, paste0(x,'.ids')))
    })
    names(nifHlike) <- c('positives','negatives','unsure')
    cat(paste(sapply(nifHlike,length), names(nifHlike), collapse=", "),"\n")
}
cat("\n")


##------------------------------------------------------------------------------
##
## Report contaminants and discard if asked.
##
if (!is.null(contamAsvs)) {
    AnnounceSection("Checking for contaminants")
    cat("Looking over contaminants in each study.\n")
    ## Hack: GatherAsvs's makeBigAbundTable.R formats the column names in the abund
    ## table as <tags>___<sample name>.  Use the <tags> as the "study".
    studies <- sapply(strsplit(colnames(atab),'___'), '[[', 1)
    idx.C <- which(rownames(atab) %in% contamAsvs)           # Big C
    totsBefore <- c(reads=sum(atab), asvs=nrow(atab))
    for (sid in unique(studies)) {
        idx.s <- which(studies == sid)                       # This study
        x <- atab[,idx.s, drop = F]
        idx.c <- which(rowSums(x[idx.C,, drop = F]) > 0)     # Little c
        if (length(idx.c) > 0) {
            idx.c <- idx.C[idx.c]
            cat(paste0(round(100*sum(x[idx.c,]) / sum(x),1),"%"),
                "of all reads from", sid, "are in", length(idx.c),
                "AUIDs that might be contaminants.")
            if (!RETAIN_CONTAMINANTS) {
                cat("Removing them.")
                stopifnot(rownames(x) == rownames(atab))
                atab[idx.c,idx.s] <- 0
            }
            cat("\n")
        } else {
            cat(sid,"has no contaminants!\n")
        }
    }
    cat("\n")
    RecordDeltasAtFilterStep("Contaminants",
                             sum(atab),  totsBefore['reads'],
                             nrow(atab), totsBefore['asvs'])
    RecordSampleReadsAfterFilterStep("Contaminants", atab)
}


##------------------------------------------------------------------------------
##
## Drop small samples and rare AUIDs
##
AnnounceSection("Dropping small samples and rare AUIDs")

## Drop samples with too few reads, and then drop AUIDs that were only in the
## dropped samples (i.e. remove rows that are all 0).
DropSmallSamples <- function(tab, minReadsInSamp = MIN_SAMPLE_SIZE)
{
    cat("Threw out samples with <=", minReadsInSamp, "reads.\n",
        "This reduced the abundance table from\n",
        dim2desc(dim(tab)), "to ")
    tab <- tab[,colSums(tab) > minReadsInSamp]
    tab <- tab[rowSums(tab) > 0,]
    cat(dim2desc(dim(tab)),"\n")
    tab
}

DropRareAuids <- function(tab, minSamps=AUID_MIN_SAMPS,
                               minReadsEachSamp=AUID_MIN_READS,
                               oneHitWonderReads=AUID_ONE_HIT_WONDER_READS)
{
    cat("Threw out AUIDs that did not have at least", minReadsEachSamp, "reads",
        "in at least", minSamps, "samples,\n",
        "or at least", oneHitWonderReads, "reads in one sample.\n",
        "This reduced the abundance table from", dim2desc(dim(tab)),"to\n")
    idx <- which(apply(tab, 1, function (rowv) {
                     (sum(rowv >= minReadsEachSamp) >= minSamps) ||
                     (any(rowv > oneHitWonderReads))
                 }))
    tab <- tab[idx,]
    tab <- tab[,colSums(tab) > 0] # dropping AUIDs could make some samples all 0
    cat(dim2desc(dim(tab)),"\n")
    tab
}

totsBefore <- c(reads=sum(atab), asvs=nrow(atab))
atab <- DropSmallSamples(atab)
RecordDeltasAtFilterStep("Drop small samples",
                         sum(atab),  totsBefore['reads'],
                         nrow(atab), totsBefore['asvs'])
RecordSampleReadsAfterFilterStep("Drop small samples", atab)

totsBefore <- c(reads=sum(atab), asvs=nrow(atab))
atab <- DropRareAuids(atab)
RecordDeltasAtFilterStep("Drop rare ASVs",
                         sum(atab),  totsBefore['reads'],
                         nrow(atab), totsBefore['asvs'])
RecordSampleReadsAfterFilterStep("Drop rare ASVs", atab)


##------------------------------------------------------------------------------
##
## Report whether AUIDs are NifH-like for each sample.
## Discard the "negatives" at the end, if asked.
##
if (!is.null(nifHlike)) {
    AnnounceSection("Applying previous results from NifH-like checker tool")
    cat("\nLooking over AUIDs in each study to see whether they are NifH-like.\n")
    ## Hack: GatherAsvs's makeBigAbundTable.R formats the column names in the abund
    ## table as <tags>___<sample name>.  Use the <tags> as the "study".
    studies <- sapply(strsplit(colnames(atab),'___'), '[[', 1)
    for (sid in unique(studies)) {
        cat(paste0("\n>>> ",sid,":\n"))
        idx.s <- which(studies == sid) # this study
        x <- atab[,idx.s, drop = F]
        idx.T <- intersect(rownames(x), nifHlike[['positives']])
        if (all(colSums(x) == colSums(x[idx.T,, drop = F]))) {
            ## For every sample, total reads equals total reads from +'s.
            cat("  + All AUIDs are NifH-like!\n")
        } else {
            stopifnot(rownames(x) == rownames(atab))
            for (typ in c('positives','unsure','negatives')) {
                idx.T <- intersect(nifHlike[[typ]], rownames(atab))
                idx.t <- which(rowSums(x[idx.T,, drop = F]) > 0)
                sym <- paste0("  ",c(positives="+",unsure="?",negatives="-")[typ])
                if (length(idx.t) > 0) {
                    idx.t <- idx.T[idx.t]
                    cat(sym, paste0(round(100*sum(x[idx.t,]) / sum(x),1),"%"),
                        "of all reads are in", length(idx.t),
                        paste0("AUIDs that have '",typ,"' similarity to NifH.\n"))
                    lenStat <- summary(id2len[idx.t])
                    lenStat <- do.call(sprintf,
                                       c(fmt="median %s nt, range %s-%s nt, IQR %s-%s nt",
                                       as.list(lenStat)[c(3,1,6,2,5)]))
                    cat(paste0("    Lengths of '",typ,"' AUIDs: "), lenStat, "\n")
                    if (FALSE && !RETAIN_NON_NIFH_LIKE && typ=="negatives") {
                        ## DISABLED: Do for all samples after the loop.
                        cat("    Removing them.\n")
                        atab[idx.t,idx.s] <- 0
                    }
                    #cat("\n")
                } else {
                    cat(sym,"No AUIDs in",sid,"are", typ, "NifH-like.\n")
                }
            }
        }
    }
    cat("\n")
    if (!RETAIN_NON_NIFH_LIKE && length(nifHlike[['negatives']]) > 0) {
        ## Drop the "negatives."
        idx <- setdiff(rownames(atab), nifHlike[['negatives']])
        if (length(idx) < nrow(atab)) {
            cat("Removing",length(idx),"negative NifH-like AUIDs.\n")
            totsBefore <- c(reads=sum(atab), asvs=nrow(atab))
            atab <- atab[idx,]
            RecordDeltasAtFilterStep("Drop NifH negatives",
                                     sum(atab),  totsBefore['reads'],
                                     nrow(atab), totsBefore['asvs'])
            RecordSampleReadsAfterFilterStep("Drop NifH negatives", atab)
        }
    } else { cat("Not filtering the NifH-like \"negatives\".\n") }
    cat("Not filtering the NifH-like \"unsures\".\n")
}


##------------------------------------------------------------------------------
##
## Length-based filtering done last so we can see impacts of the
## approaches above.
##
AnnounceSection("Length-based filtering")

## Returns dataframe of per sample impacts of filtering: % reads kept, lost.
IfFilterByLens <- function(atab, lmin, lmax, wantIdsIN=FALSE)
{
    asvs <- rownames(atab)
    idx.in  <- intersect(asvs, names(which(  ((lmin < id2len) & (id2len < lmax)))))
    if (wantIdsIN) {
        return(idx.in)
    }
    idx.out <- intersect(asvs, names(which(! ((lmin < id2len) & (id2len < lmax)))))
    tots <- colSums(atab)
    x <- data.frame(IN  = colSums(atab[idx.in,])/tots,
                    OUT = colSums(atab[idx.out,])/tots)
    stopifnot( abs(x$IN + x$OUT -1) < 0.0001 )
    round(100*x, 2)
}

cat("Checking impacts of filtering with different length requirements.")
lenRanges <- data.frame(min=seq(240,290,10), max=seq(390,340,-10))
x <- apply(lenRanges, 1,
           function(lf) {
                  cat(".")
                  x <- IfFilterByLens(atab, lf[1], lf[2])
                  c("Num samps with >5% read loss"   = nrow(subset(x, OUT > 5)),
                    "Num samps with >10% read loss"  = nrow(subset(x, OUT > 10)),
                    "Num samps with >25% read loss"  = nrow(subset(x, OUT > 25)))
                 #  "Pct samps with >90% reads kept" = round(mean(subset(x, IN > 90)$IN),1),
                 #  "Pct samps with >95% reads kept" = round(mean(subset(x, IN > 95)$IN),1),
                 #  "Pct samps with >99% reads kept" = round(mean(subset(x, IN > 99)$IN),1))
           })
cat("done.\n")
colnames(x) <- apply(lenRanges,1,paste, collapse='-')
print(x)
cat("\n")

cat("Filtering at length range",paste(lenFiltRange, collapse='-'),"nt as you asked.\n")
totsBefore <- c(reads=sum(atab), asvs=nrow(atab))
x <- IfFilterByLens(atab, lenFiltRange[1], lenFiltRange[2])
atab <- atab[IfFilterByLens(atab, lenFiltRange[1], lenFiltRange[2], TRUE), ]
RecordDeltasAtFilterStep("Length based",
                         sum(atab),  totsBefore['reads'],
                         nrow(atab), totsBefore['asvs'])
RecordSampleReadsAfterFilterStep("Length based", atab)

## This file is somewhat redundant now that gFilterStages.samples gets written out
## at the end of the script.
f <- paste0("lenFilterImpactsToSamps.",paste(lenFiltRange, collapse='.'),".tsv")
write.table(x, f, sep="\t", quote=F)
cat("Wrote table",f,"of % reads kept/lost for each sample.\n")
rm(x,f,lenRanges)


##------------------------------------------------------------------------------
##
## Wrap up
##
AnnounceSection("Wrap-up")

cat("Impacts of each filter step on reads and ASVs:\n")
row.names(gFilterStages) <- NULL
print(gFilterStages)
cat("\n")
## Useful to have this for summary stats/plots.
write.table(gFilterStages, 'filterAuids_byStages.tsv', quote=FALSE, sep="\t")

## Save a table of all the counts retained in each sample at each stage.
write.table(gFilterStages.samples, 'filterAuids_samplesByStages.tsv', quote=F, sep="\t")

cat("Writing out", dim2desc(dim(atab)), "abundance table",FOUT,"...")
gcon <- gzfile(FOUT, "wb")
write.table(atab, gcon, quote=FALSE, sep="\t")
close(gcon)
cat("Done!\n")
