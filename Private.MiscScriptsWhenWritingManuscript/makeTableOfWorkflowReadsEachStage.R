#!/usr/bin/env Rscript

##
## Make tables that shows how many reads are retained through each stage of the post-DADA2-pipeline
## workflow: by sample, by study.  Also a multi-box plot faceted by study.
##

## Set this as in make_nifH_ASV_database.R.  By July 2023 we decided this was TRUE.
KEEP_SAMPS_MISSING_META_OR_CMAP_DATA <- TRUE

## Create descriptions for the different stages of the workflow.  For clarity do this here even
## though stageRename is only needed in the plotting code far below.  Note the 'hack' to
## conditionally include filtering steps hasSampMeta and hasEnvMeta.
stageRename <- c(ReadsPipeline                     = 'From pipeline',
                 ReadsGatherAsvs                   = 'GatherAsvs',
                 ReadsFilterAuids.SmallSamp        = 'Undersequenced samples',
                 ReadsFilterAuids.Rare             = 'Rare',
                 ReadsFilterAuids.NonNifH          = 'Not nifH-like',
                 ReadsFilterAuids.Length           = 'Too short or long',
                 ReadsWorkspaceStartup.notEmpty    = 'Empty samples',
                 hack_hasSampMeta                  = 'No metadata',
                 hack_hasEnvMeta                   = 'No CMAP data',
                 ReadsWorkspaceStartup.noAnnot     = 'No annotation')
if (KEEP_SAMPS_MISSING_META_OR_CMAP_DATA) {
    ## make_nifH_ASV_database.R didn't drop samples that lacked meta or CMAP data.
    stageRename <- stageRename[-grep('^hack', names(stageRename))]
} else {
    ## We DID drop samples that lacked meta or CMAP data so fixup those entries.
    names(stageRename) <- sub('^hack_', 'ReadsWorkspaceStartup.', names(stageRename))
}
stageRename <- rev(stageRename)  # plot stages top-to-bottom


options(width=120)
library(ggplot2)
library(reshape2)

workflowDir <- commandArgs(T)
if (length(workflowDir) != 1 || !dir.exists(workflowDir)) {
    stop("Pass a directory that contains the workflow stages, i.e. ",
         "with GatherAsvs, FilterAuids, and other subdirectories. ",
         "You may pass '..' or relative paths.\n")
}
# All workflow stages present?  (Ignoring Metadata and CMAP env data.)
stageNames <- c("GatherAsvs","FilterAuids","AnnotateAuids","GatherMetadata","WorkspaceStartup")
stopifnot( stageNames %in% list.files(workflowDir) )

StageDir  <- function(sn)    { normalizePath(file.path(workflowDir,sn),    mustWork=T)  }
StageFile <- function(sn,fn) { normalizePath(file.path(workflowDir,sn,fn), mustWork=T)  }


##################################################
## From pipeline.

## The input table for GatherAsvs points to the paired FASTA / abundance tables
gatherTabTsv <- StageFile('GatherAsvs','asvs.noChimera.fasta_table.tsv')
cat("Loading table of pipeline outputs", gatherTabTsv, "\n")
gatherTab <- read.table(gatherTabTsv)
## Replace FASTA name in col 2 with abundance table name.  Make all paths
## absolute, i.e.  convert those that are relative to GatherAsvs.
gatherTab[,2] <- sub('\\.fasta$','.tsv',gatherTab[,2])
colnames(gatherTab) <- c('Study', 'ASV_noChimera_tsv', 'tag')
relPathIdx <- grep('^[~/]', gatherTab$ASV_noChimera_tsv, invert=T)
if (length(relPathIdx) > 0) {
    gatherTab$ASV_noChimera_tsv[relPathIdx] <- file.path(dirname(gatherTabTsv),
                                                         gatherTab$ASV_noChimera_tsv[relPathIdx])
}
gatherTab$ASV_noChimera_tsv <- normalizePath(gatherTab$ASV_noChimera_tsv, mustWork=T)
rm(relPathIdx)
stopifnot(file.exists(gatherTab$ASV_noChimera_tsv))
## Format of the tag is <studyID>:<processingGroupNum>:<SRA study>


## Count total reads in each sample in an abundance table from the pipeline.
## Returned data frame rows are:  sampId, studyName, tag, totReads
CountReadsInPipelineAbundanceTable <- function(atabPath, studyName, tag, verbose=F, onlyTheseAuids = NULL)
{
    if (verbose) cat("Counting reads from", atabPath, "\n")
    readsBySample <- colSums(read.table(atabPath, check.names=T))
    ## We need a sample name that will match the sample names used in the
    ## workflow (created during GatherAsvs) which consist of the tag (slightly
    ## modified), ___, then the sample name.
    wfSampNams <- paste0(gsub(':','.',tag), "___", names(readsBySample))
    data.frame(SampleWF = wfSampNams, Sample = names(readsBySample),
               Study  = studyName, tag = tag, ReadsPipeline  = readsBySample,
               row.names=NULL)
}

## Initialize the workflow table with all the read counts from the pipeline.
x <- apply(gatherTab, 1, function(rvec) {
         CountReadsInPipelineAbundanceTable(rvec['ASV_noChimera_tsv'], rvec['Study'], rvec['tag'])
     })
workflowTable <- do.call(rbind, x)
rm(x)

## Helper which loads an abundance table from a workflow stage and appends its
## samples' total reads as a new column to the growing workflow table.
MergeAbundTable <- function(wftab, atabPath, colName)
{
    cat("Loading the", colName, "abundance table\n\t",atabPath,"\n...")
    atab <- read.table(atabPath, check.names=T)
    cat("done.  Now merging it to the big table.\n")
    x <- data.frame(SampleWF = colnames(atab), newCol = colSums(atab))
    stopifnot( x$newCol >= 0 ) # no NA's
    stopifnot(x$SampleWF %in% wftab$SampleWF)  # atab must be a subset of wft
    x <- merge(wftab, x, by='SampleWF', all=T)
    stopifnot(nrow(x) == nrow(wftab))
    ## If merge added any NA's in the new col, make them 0.
    x[is.na(x$newCol),'newCol'] <- 0
    ## Reads cannot increase from stage i to i+1
    stopifnot( colnames(x)[ncol(x)] == 'newCol' )
    stopifnot(is.numeric(x[,ncol(x)]) && is.numeric(x[,ncol(x)-1]))
    idx <- which(x[,ncol(x)] > x[,ncol(x)-1])
    if (length(idx) > 0) {
        cat("WARNING: These", length(idx),"samples somehow gained reads at stage",colName,":\n")
        print(x[idx,])
    }
    ##stopifnot( sum(x[,ncol(x)] > x[,ncol(x)-1], na.rm=T) == 0 ) ## na.rm should not be needed.
    colnames(x) <- c(colnames(x)[-ncol(x)],colName)
    x
}


##################################################
## GatherAsvs (end of):  Chimeras dropped
x <- MergeAbundTable(workflowTable,
                     StageFile('GatherAsvs','asv2auid.abundances.tsv.gz'),
                     'ReadsGatherAsvs')
workflowTable <- x


##################################################
## FilterAuids
## March 21, 2024, I added a table that records for each sample the reads retained at
## the end of each FilterAuids stage. So no longer neeeded to back-calculate, yay!
##   1. Report contaminants but do not remove them (b/c RETAIN_CONTAMINANTS is T)
##   2. Drop small samples
##   3. Drop rare AUIDs
##   4. Drop non-nifH-like AUIDs (b/c RETAIN_NON_NIFH_LIKE is F)
##   5. Length-based filtering
filtTab <- read.table(StageFile('FilterAuids','filterAuids_samplesByStages.tsv'),
                      sep="\t", header=T)
stopifnot(!is.na(filtTab)) # Expect 0's in stage i+1,2,3 if a sample was dropped at stage i
stopifnot(colnames(filtTab) == c("Sample", "Initial", "Contaminants", "Drop.small.samples",
                                 "Drop.rare.ASVs", "Drop.NifH.negatives", "Length.based"))
filtTab <- filtTab[,-c(2,3)]  # Initial is same as GatherAsvs, and we don't drop Contaminants.
## Shorten/simplify names
colnames(filtTab) <- c('SampleWF',
                       paste0('ReadsFilterAuids.', c('SmallSamp','Rare','NonNifH','Length')))

## Merge
x <- merge(workflowTable, filtTab, by='SampleWF', all=T)
x[is.na(x)] <- 0
workflowTable <- x
rm(filtTab, x)


##################################################
## AnnotateAuids <-- No filtering done here.


##################################################
## GatherMetadata <-- No filtering.
##   Identifies metaT samples. These used to get dropped in WorkspaceStartup but no longer.
metaTsamples <- readLines(StageFile('GatherMetadata','metatranscriptomic_samples.txt'))


##################################################
## CMAP <-- No filtering of AUIDs.


##################################################
## WorkspaceStartup (end of)

## See WorkspaceStartup's Makefile and make_nifH_ASV_database.R for filtering steps in order.
##  - Drop samples with 0 reads.
##  - Show effects of removing samples with no sample metadata, (conditionally)
##  - Show effects of removing samples with no environmental metadata. (conditionally)
##  - Show  effects of removing ASVs that have no annotation.
##

## Drop samples with 0 reads.  Nothing to do but copy the previous stage (which
## might have 0 reads).
workflowTable$ReadsWorkspaceStartup.notEmpty <- workflowTable[,ncol(workflowTable)]

## Drop specified samples from the passed workflow table.  Duplicate the rightmost column in wftab
## but name it colName and make the 'samps' have 0 total counts.
DropSamps <- function(wftab, samps, colName)
{
    ## Use grep so that we can find the sample even if R prefixed with an X to make it a valid
    ## (column) name.  Note that some of the 'samps' might have already been dropped from wftab for
    ## other reasons.
    idx <- unlist(sapply(samps, function(s) grep(s, wftab$Sample)))
    stopifnot( names(idx) %in% samps )
    idx <- as.vector(idx)  # stress that values, not names, are used next
    x <- wftab[,ncol(wftab)]  # Last stage's total
    stopifnot(idx %in% 1:length(x))
    x[idx] <- 0  # Sample removal means 0-ing the reads
    wftab[,colName] <- x
    wftab
}


if (!KEEP_SAMPS_MISSING_META_OR_CMAP_DATA) {
    ## Now remove samples that had no sample metadata
    dropme <- readLines(StageFile('WorkspaceStartup','sampsWithoutMetadata.txt'))
    x <- DropSamps(workflowTable, dropme, 'ReadsWorkspaceStartup.hasSampMeta')
    workflowTable <- x

    ## Now remove samples that had no environmental metadata
    dropme <- readLines(StageFile('WorkspaceStartup','sampsWithoutEnvdata.txt'))
    x <- DropSamps(workflowTable, dropme, 'ReadsWorkspaceStartup.hasEnvMeta')
    workflowTable <- x
}


## Get for each sample the number of reads in the final abundance table.
## Because dropping of ASVs with no annotation is the final filtering step in
## WorkspaceStartup, these read counts minus those in the rightmost column of
## workflowTable equal the reads lost due to lack of annotation.
ReadsInFinalAbundanceTable <- function(colName = "final_reads")
{
    load(StageFile('WorkspaceStartup','workspace.RData'))
    stopifnot(colnames(abundTab)[1] == 'AUID')
    sampTots <- colSums(abundTab[,-1])
    df <- data.frame(Sample = names(sampTots), tots = as.numeric(sampTots), row.names = NULL)
    colnames(df) <- c('Sample', colName)
    df
}

y <- ReadsInFinalAbundanceTable("ReadsWorkspaceStartup.noAnnot")
## The workflowTable might have samples not in the final abundance table so all.x = T.
workflowTable <- merge(workflowTable, y, by = 'Sample', all.x = T)
idx <- which(is.na(workflowTable$ReadsWorkspaceStartup.noAnnot))
workflowTable$ReadsWorkspaceStartup.noAnnot[idx] <- 0
rm(y,idx)


cat("Checking that all 'Reads' columns in the table are numeric.\n")
idx <- grep('^Reads',colnames(workflowTable))
stopifnot( apply(workflowTable[,idx], 2, is.numeric) )
cat("Checking for every sample that read counts decrease from each stage i to i+1 through the workflow.\n")
##stopifnot( apply(workflowTable[,idx], 1, function(v) all(order(v,decreasing=T)==1:length(v))) )
## There are a very few cases where reads increase by 1 or 2 from ReadsGatherAsvs to ReadsFilterAuids.Rare.
## I think this stems from the round() to calculate endOf3. Just report it.
idx <- which(!apply(workflowTable[,idx], 1, function(v) all(order(v,decreasing=T)==1:length(v))))
if (length(idx) > 0) {
    cat("These have some strange increases in reads.\n")
    print(workflowTable[idx,])
}

cat("\nHere are the total reads and % retained at each stage:\n")
x <- colSums(workflowTable[,-c(1:4)])
data.frame(Reads = x, Pct = 100*x/x[1])
rm(x)


cat("\nWriting out table of reads retained at each stage:  workflowTable.csv\n")
write.csv(format(workflowTable, digits=2), file='workflowTable.csv', row.names=F, quote=F)


## Also prepare Table 5 for the manuscript which has mean reads (over samples) retained for each
## study at each stage of the workflow, and some overall stats in the last three rows (across
## studies).
dat <- read.csv('workflowTable.csv')[,-c(1,2,4)]
dat$PctRetained <- 100 * dat$ReadsWorkspaceStartup.noAnnot / dat$ReadsPipeline
dat$PctRetained[is.nan(dat$PctRetained)] <- 0

## Top rows of table will have the per study stats:
## For each stage (column) take the mean but only using samples that have reads.  Alternatively
## could restrict to just the samples that make it into the ASV database, but some of them have 0
## reads because after unannotated ASVs were dropped empty samples were not scrubbed. Perhaps
## there's some value in having these empty samples in the DB (e.g. in the CMAP and meta data) so I
## have at this time decided not to scrub them.
studyMeans <- aggregate(.~Study, dat, function (v) mean(v[v>0]))

## Last three rows have mean, median, and sums taken across all samples:
overall <- matrix(c(apply(dat[,-1], 2, function (v) mean(v[v>0])),
                    apply(dat[,-1], 2, function (v) median(v[v>0])),
                    colSums(dat[,-1])),
                  nrow = 3, byrow=T,
                  dimnames = list(c('mean','median','sum'), colnames(dat)[-1]))
overall['sum','PctRetained'] <- NA

## Construct the final table
dat <- rbind(studyMeans,
             data.frame(Study = rownames(overall), overall))
rownames(dat) <- NULL
cat("\nWriting out table_5_for_manuscript.csv\n")
write.csv(format(dat, digits=2), file='table_5_for_manuscript.csv', row.names=F, quote=F)
rm(dat, overall, studyMeans)


##################################################
## Summary plot

cat("Preparing plot that shows for each study how the workflow progressed",
    "(reads retained at each stage).\n")

wft <- workflowTable

## Simplest if drop a few empty-from-pipes now.  'notEmpty' can still show >0 in plot.
idx <- which(wft$ReadsPipeline > 0)
cat(length(idx),"of",nrow(wft),"samples have >0 reads from the pipeline.  Using only them.\n")
wft <- wft[idx,]

## When WorkspaceStartup drops empty samples, of course that does not change the number of reads
## retained by any sample, as next line verifies.
stopifnot( wft$ReadsWorkspaceStartup.notEmpty - wft$ReadsFilterAuids.Length == 0 )
## Therefore, it is not helpful to show ReadsWorkspaceStartup.notEmpty in the plots.
wft <- wft[, colnames(wft) != 'ReadsWorkspaceStartup.notEmpty']
stageRename <- stageRename[names(stageRename) != 'ReadsWorkspaceStartup.notEmpty']

cat("Preparing data for ggplot...\n")
## Have reads at each stage.  Want the % of reads lost at each stage. Get it
## with diff().  Note that diff() removes the leftmost col, ReadsPipeline, which
## is good because it would always be 100%.
x <- as.matrix(wft[,rev(names(stageRename))])  # order for diff!
x <- t(apply(x, 1, function(v) { 100 * (-diff(v) / v['ReadsPipeline']) }))
dat <- cbind(wft[,c('Sample','Study')], x)
dat.m <- melt(dat, id=c('Sample','Study'))
dat.m$variable <- factor(dat.m$variable,
                         levels = names(stageRename), labels = as.character(stageRename))
colnames(dat.m) <- c('Sample','Study','Filter','value')
rm(x,dat)
## Above we reported a few cases where 1-2 reads were gained, probably a bug in this
## script which must use approx %'s to estimate length-based filtering.  Set them to 0.
dat.m[which(dat.m$value < 0),'value'] <- 0

## Which studies contribute the most data (after all filtering)?
tots <- aggregate(wft[,ncol(wft)], by=list(study=wft$Study), sum, na.rm=T)
colnames(tots) <- c('study','totreads')
tots <- tots[order(tots$totreads,decreasing=T),]
tots$pct <- round(100*tots$totreads/ sum(tots$totreads),1)
idx <- match(dat.m$Study, tots$study)
stopifnot(!is.na(idx))
## Bracketted numbers in facet labels show the contribution [as a %] of each study to
## the total reads that made it through FilterAuids (the last column in wft).
dat.m$Study <-factor(dat.m$Study,
                     levels = tots$study,
                     labels = paste0(tots$study,' [',tots$pct,'%]'))
cat("Here are the total read contributions by study:\n")
print(tots)
rm(tots,idx)

cat("Plotting...\n")

## The plot shows at each stage, for each data set (shown as jittered points)
## the % of _all_ reads from the pipeline that were discarded *AT* the indicated
## stage.  Stages happen sequentially (top-to-bottom rows in plot) so each stage
## should be interpreted as the *additional* reads that are lost.  (One cannot
## tell if a read that was lost due to being in a chimera would have later been
## lost for lack of CMAP data.)
## Very thin violins so the jitter is important.
sfb <- scale_fill_brewer(palette = 'Set3', type = 'qual')
scb <- scale_color_brewer(palette = 'Set3', type = 'qual')
g.box <- ggplot(dat.m, aes(Filter, value, fill=Filter)) +
           geom_violin() + geom_jitter(aes(color=Filter), size=0.5, alpha=0.5) +
           coord_flip() + sfb + scb +
           theme_bw() + labs(x=NULL,y='% pipeline reads lost at stage') +
           guides(color = 'none') +  # color legend (jitter) is redundant
           guides(fill = guide_legend(reverse = TRUE)) +
           facet_wrap(~Study)
ggsave(filename='workflowReadsRetainedEachStage.png', plot=g.box, width=11.5, height=7, units='in', dpi=144)

## Also make an overall plot (all samples across all studies).  Drop legend and
## labels since will manually layout to use those from plot above.
g.box <- ggplot(dat.m, aes(Filter, value, fill=Filter)) +
           geom_violin() + coord_flip() + sfb +
           theme_bw() + theme(axis.text.y = element_blank()) +
           labs(x = element_blank(), y = element_blank()) +
           guides(fill = 'none')   # Use legend of main plot
ggsave(filename='workflowReadsRetainedEachStage_overall_violin.png', plot=g.box,
       width=11.5/6, height=7/4, units='in', dpi=144)

## Box plot might look better.
g.box <- ggplot(dat.m, aes(Filter, value, fill=Filter)) +
           geom_boxplot(outlier.size=0.5, outlier.alpha=0.5) + coord_flip() + sfb +
           theme_bw() + theme(axis.text.y = element_blank()) +
           labs(x = element_blank(), y = element_blank()) +
           guides(fill = 'none')   # Use legend of main plot
ggsave(filename='workflowReadsRetainedEachStage_overall_box.png', plot=g.box,
       width=11.5/6, height=7/4, units='in', dpi=144)

cat("Done. See workflowReadsRetainedEachStage.png and the 'overall' versions.\n",
    "Note that stage ReadsWorkspaceStartup.notEmpty is not shown because it discards\n",
    "empty samples and does not impact the read counts.\n")
quit(save='no')
