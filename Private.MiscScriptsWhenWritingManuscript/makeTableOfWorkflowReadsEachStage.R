#!/usr/bin/env Rscript

## Modified because now we changed the post-pipeline to not drop metaT samples,
## and because I want to run this script without yet having WorkspaceStartup
## ready for GitHub.  Also fixed a bug (in this script) in tracking of
## FilterAuids %'s.

##
## Make tables that shows how many reads are retained through each stage of the post-DADA2-pipeline
## workflow: by sample, by study.  Also a multi-box plot faceted by study.
##

options(width=120)
library(ggplot2)
library(reshape2)

workflowDir <- commandArgs(T)
if (length(workflowDir) != 1 || !dir.exists(workflowDir)) {
    stop("Pass a directory that has contains the workflow stages, i.e. ",
         "with GatherAsvs, FilterAuids, and other subdirectories. ",
         "You may pass '..' or relative paths.\n")
}
# All workflow stages present?  (Ignoring Metadata and CMAP env data.)
stageNames <- c("GatherAsvs","FilterAuids","AnnotateAuids","GatherMetadata") #,"WorkspaceStartup")
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
##   See FilterAuids/scripts/filterAuids.R for the order of filters.
##   1. Report contaminants but do not remove them (b/c RETAIN_CONTAMINANTS is T)
##   2. Drop small samples and rare AUIDs
##   3. Drop non-nifH-like AUIDs (b/c RETAIN_NON_NIFH_LIKE is F)
##   4. Length-based filtering
## Unfortunately, FilterAuids does not record reads retained for each sample and
## stage. (filterAuids_byStages.tsv is overall.)  It does output
## lenFilterImpactsToSamps.280.360.tsv, which I can used to back calculate the
## num reads at the end of step #3.  And I can identify non-nifH AUIDs using
## FilterAuids/Check_NifH_like/negatives.ids.

filtAbunds <- read.table(StageFile('FilterAuids','auid.abundances.filtered.tsv.gz'))
lenFiltTab <- read.table(StageFile('FilterAuids','lenFilterImpactsToSamps.280.360.tsv'))
stopifnot(rownames(lenFiltTab) %in% workflowTable$SampleWF)      # Reverse is not true
stopifnot(setequal(rownames(lenFiltTab), colnames(filtAbunds)))  # Identical samples
negIds <- readLines(StageFile('FilterAuids', file.path('Check_NifH_like','negatives.ids')))

x <- data.frame(endOf4 = colSums(filtAbunds)[rownames(lenFiltTab)])  # tot reads each sample
x$endOf3 <- round(x$endOf4 / (lenFiltTab$IN/100))  # % retained in #4 tells us num at end of #3

## No great way to get to reads at end of #2 unless re-read this table.
## That gets us the reads in non-nifH-like AUIDs, which we sum by sample
nonNifH <- read.table(StageFile('GatherAsvs','asv2auid.abundances.tsv.gz'))[negIds,rownames(x)]
stopifnot(!is.na(nonNifH))
nonNifH <- colSums(nonNifH)  # by sample
stopifnot( names(nonNifH) == rownames(x) )
x$endOf2 <- x$endOf3 + nonNifH  # add back the non-nifH-like reads that got dropped in #3

## Improve the names, and set up column order before the merge.
x <- x[,c('endOf2','endOf3','endOf4')]
colnames(x) <- sub('endOf4', 'ReadsFilterAuids.Length',
                   sub('endOf3', 'ReadsFilterAuids.NonNifH',
                       sub('endOf2', 'ReadsFilterAuids.Rare', colnames(x))))
## Check strictly non-increasing
x[is.na(x)] <- 0
stopifnot(sum(x[,1] < x[,2]) == 0)
stopifnot(sum(x[,2] < x[,3]) == 0)

## Merge
x$SampleWF <- rownames(x)
x <- merge(workflowTable, x, by='SampleWF', all=T)
x[is.na(x)] <- 0
workflowTable <- x

## Clean up
rm(filtAbunds, lenFiltTab, negIds, nonNifH, x)


##################################################
## AnnotateAuids <-- No filtering.


##################################################
## GatherMetadata <-- No filtering.
##   Identifies metaT samples which get dropped in WorkspaceStartup
##droppedMetaT <- readLines(StageFile('GatherMetadata','dropped_metatranscriptomic_samples.txt'))
## 25 July 2023: They don't get dropped.
metaTsamples <- readLines(StageFile('GatherMetadata','metatranscriptomic_samples.txt'))


##################################################
## CMAP <-- No filtering.


##################################################
## WorkspaceStartup (end of)

##  - Skip: It combines samples from the same date/lat/lon (but not replicates).
##    Ignore because I don't want to merge rows (samples), and combining
##    preserve the total number of reads in the combined samples (by default --
##    see docs for combineSamples.R).
##  - Drop samples with 0 reads.
##  - Show effects of removing metaT samples
##  - Show effects of removing samples with no sample metadata,
##  - ... or  environmental metadata.
## See Makefile and makeWorkspace.R to get proper order to show.
##

## Drop samples with 0 reads.  Nothing to do but copy the previous stage (which
## might have 0 reads).
if (F) {
## 25 July 2023: Leave out until have this.
workflowTable$ReadsWorkspaceStartup.notEmpty <- workflowTable[,ncol(workflowTable)]
}

## Drop specified samples, e.g. samps can be metatranscriptomic.
DropSamps <- function(wftab, samps, colName)
{
    ## Use grep so that we can find the sample even if R prefixed with an X to
    ## make it a valid (column) name.  Note that some of the 'samps' might have
    ## already been dropped from wftab for other reasons.
    idx <- unlist(sapply(samps, function(s) grep(s, wftab$Sample)))
    stopifnot( names(idx) %in% samps )
    idx <- as.vector(idx)  # stress that values, not names, are used next
    x <- wftab[,ncol(wftab)]  # Last stage's total
    stopifnot(idx %in% 1:length(x))
    x[idx] <- 0  # Sample removal means 0-ing the reads
    wftab[,colName] <- x
    wftab
}

## 25 July 2023: Keep the metaT samples in the table.
if (F) {
x <- DropSamps(workflowTable, droppedMetaT, 'ReadsWorkspaceStartup.genes')
workflowTable <- x
}

## Now remove samples that had no sample metadata
## 25 July 2023: Cannot do this because WorkspaceStartup/sampsWithoutMetadata.txt does
## not yet exist in the GitHub version of the post-pipeline.
if (F) {
dropme <- readLines(StageFile('WorkspaceStartup','sampsWithoutMetadata.txt'))
x <- DropSamps(workflowTable, dropme, 'ReadsWorkspaceStartup.hasSampMeta')
workflowTable <- x
}

## Now remove samples that had no environmental metadata
## 25 July 2023: Same issue as above
if (F) {
dropme <- readLines(StageFile('WorkspaceStartup','sampsWithoutEnvdata.txt'))
x <- DropSamps(workflowTable, dropme, 'ReadsWorkspaceStartup.hasEnvMeta')
workflowTable <- x
}

## FIXME: <-- old. still true?
## The total number of reads in the whole data set (sum of the last column just
## appended) is 18456 more than then sum of the abundTab in the
## ../WorkspaceStartup/globalNcdWorkspace.RData.  That's a small fraction but
## I'd like to know why.

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
cat("If you want to read this table in R, you should specify column classes:\n",
    " wftab <- read.csv('workflowTable.csv', colClasses=c( rep(NA,4), rep('integer',8)))\n",
    "or else some of the Reads count columns might be characters.  Not sure why.\n")
##save.image('workflowReadsAtEachStage.RData')


##################################################
## Summary plot

cat("Preparing plot that shows for each study how the workflow progressed",
    "(reads retained at each stage).\n")

## In the plot use more intuitive, shorter names than in workflowTable. For now...
stageRename <- c(ReadsPipeline                     = 'From pipeline',
                 ReadsGatherAsvs                   = 'GatherAsvs',
                 ReadsFilterAuids.Rare             = 'Rare',
                 ReadsFilterAuids.NonNifH          = 'Not nifH-like',
                 ReadsFilterAuids.Length           = 'Too short or long')
## 25 July 2023: Leave out until renabled above.
#                 ReadsWorkspaceStartup.notEmpty    = 'Empty samples',
#                 ReadsWorkspaceStartup.genes       = 'mRNA samples',
#                 ReadsWorkspaceStartup.hasSampMeta = 'No metadata',
#                 ReadsWorkspaceStartup.hasEnvMeta  = 'No CMAP data')
stageRename <- rev(stageRename)  # plot stages top-to-bottom

wft <- workflowTable
## Simplest if drop a few empty-from-pipes now.  'notEmpty' can still show >0 in plot.
idx <- which(wft$ReadsPipeline > 0)
cat(length(idx),"of",nrow(wft),"samples have >0 reads from the pipeline.  Using only them.\n")
wft <- wft[idx,]

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
## 25 July 2023: Oops, instead use last stage (column)
##tots <- aggregate(wft$ReadsWorkspaceStartup.hasEnvMeta, by=list(study=wft$Study), sum, na.rm=T)
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

cat("Done. See workflowReadsRetainedEachStage.png and the 'overall' versions.\n")
quit(save='no')
