#!/usr/bin/env Rscript

## NOT DONE.  Adapting 'Reads' namesake to count AUIDs.

## Make tables that shows how many ASVs/AUIDs are retained through each stage of the
## post-DADA2-pipeline workflow: by sample, by study.  Also a multi-box plot faceted by study.

## Set this as in make_nifH_ASV_database.R.  By July 2023 we decided this was TRUE.
KEEP_SAMPS_MISSING_META_OR_CMAP_DATA <- TRUE

## Create descriptions for the different stages of the workflow.  For clarity do this here even
## though stageRename is only needed in the plotting code far below.  Note the 'hack' to
## conditionally include filtering steps hasSampMeta and hasEnvMeta.
stageRename <- c(ASVsPipeline                     = 'From pipeline',
                 ASVsGatherAsvs                   = 'GatherAsvs',
                 ASVsFilterAuids.RareAndNonNifH   = 'Rare, then not nifH-like',  # don't know the rare ASVs unfortunately.
                 ASVsFilterAuids.Length           = 'Too short or long',
                 ASVsWorkspaceStartup.notEmpty    = 'Empty samples',
                 hack_hasSampMeta                 = 'No metadata',
                 hack_hasEnvMeta                  = 'No CMAP data',
                 ASVsWorkspaceStartup.noAnnot     = 'No annotation')
if (KEEP_SAMPS_MISSING_META_OR_CMAP_DATA) {
    ## make_nifH_ASV_database.R didn't drop samples that lacked meta or CMAP data.
    stageRename <- stageRename[-grep('^hack', names(stageRename))]
} else {
    ## We DID drop samples that lacked meta or CMAP data so fixup those entries.
    names(stageRename) <- sub('^hack_', 'ASVsWorkspaceStartup.', names(stageRename))
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


## Given an abundance table (ASV rows X sample columns, count for each
## sample the number of detected ASVs.  
CountDetectedASVsInSamples <- function(atab, minDet = 1)
{
    colSums(atab >= minDet, na.rm=T)
}


## Count total ASVs in each sample in an abundance table from the pipeline.
## Returned data frame rows are:  sampId, studyName, tag, totASVs
CountASVsInPipelineAbundanceTable <- function(atabPath, studyName, tag, verbose=F, onlyTheseAuids = NULL)
{
    if (verbose) cat("Counting ASVs detected in each sample from", atabPath, "\n")
    asvsBySample <- CountDetectedASVsInSamples(read.table(atabPath, check.names=T))
    ## We need a sample name that will match the sample names used in the
    ## workflow (created during GatherAsvs) which consist of the tag (slightly
    ## modified), ___, then the sample name.
    wfSampNams <- paste0(gsub(':','.',tag), "___", names(asvsBySample))
    data.frame(SampleWF = wfSampNams, Sample = names(asvsBySample),
               Study  = studyName, tag = tag, ASVsPipeline  = asvsBySample,
               row.names=NULL)
}

## Initialize the workflow table with all the read counts from the pipeline.
x <- apply(gatherTab, 1, function(rvec) {
         CountASVsInPipelineAbundanceTable(rvec['ASV_noChimera_tsv'], rvec['Study'], rvec['tag'])
     })
workflowTable <- do.call(rbind, x)
rm(x)


## Helper which loads an abundance table from a workflow stage and appends its
## samples' total ASVs as a new column to the growing workflow table.
MergeAbundTable <- function(wftab, atabPath, colName)
{
    cat("Loading the", colName, "abundance table\n\t",atabPath,"\n...")
    atab <- read.table(atabPath, check.names=T)
    cat("done.  Now merging its detected ASVs to the big table.\n")
    x <- data.frame(SampleWF = colnames(atab),
                    newCol = CountDetectedASVsInSamples(atab))
    stopifnot( x$newCol >= 0 ) # no NA's
    stopifnot(x$SampleWF %in% wftab$SampleWF)  # atab must be a subset of wft
    x <- merge(wftab, x, by='SampleWF', all=T)
    stopifnot(nrow(x) == nrow(wftab))
    ## If merge added any NA's in the new col, make them 0.
    x[is.na(x$newCol),'newCol'] <- 0
    ## ASV counts cannot increase from stage i to i+1
    stopifnot( colnames(x)[ncol(x)] == 'newCol' )
    stopifnot(is.numeric(x[,ncol(x)]) && is.numeric(x[,ncol(x)-1]))
    idx <- which(x[,ncol(x)] > x[,ncol(x)-1])
    if (length(idx) > 0) {
        cat("WARNING: These", length(idx),"samples somehow gained ASVs at stage",colName,":\n")
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
                     'ASVsGatherAsvs')
workflowTable <- x


##################################################
## FilterAuids
##   See FilterAuids/scripts/filterAuids.R for the order of filters.
##   1. Report contaminants but do not remove them (b/c RETAIN_CONTAMINANTS is T)
##  *2. Drop small samples and rare AUIDs
##  *3. Drop non-nifH-like AUIDs (b/c RETAIN_NON_NIFH_LIKE is F)
##   4. Length-based filtering
## Unfortunately, FilterAuids does not record reads or ASVs retained for each sample and
## stage. (filterAuids_byStages.tsv is overall.)  Unlike in the 'reads' script I cannot
## use lenFilterImpactsToSamps.280.360.tsv to calculate #2 alone.

filtAbunds <- read.table(StageFile('FilterAuids','auid.abundances.filtered.tsv.gz'))
## Cannot *use* this table b/c it records only read %'s, not ASVs.
lenFiltTab <- read.table(StageFile('FilterAuids','lenFilterImpactsToSamps.280.360.tsv'))
stopifnot(rownames(lenFiltTab) %in% workflowTable$SampleWF)      # Reverse is not true
stopifnot(setequal(rownames(lenFiltTab), colnames(filtAbunds)))  # Identical samples
negIds <- readLines(StageFile('FilterAuids', file.path('Check_NifH_like','negatives.ids')))

##OLD##x <- data.frame(endOf4 = colSums(filtAbunds)[rownames(lenFiltTab)])  # tot reads each sample
asvsBySample <- CountDetectedASVsInSamples(filtAbunds)
x <- data.frame(endOf4 = asvsBySample[rownames(lenFiltTab)])  # tot ASVs each sample

##OLD##x$endOf3 <- round(x$endOf4 / (lenFiltTab$IN/100))  # % retained in #4 tells us num at end of #3
## Darn, I cannot use lenFiltTab b/c it's %'s are of reads, not ASVs.
## Instead just get this from the non-NifH results just below.

## To get reads at end of #3.
nonNifH <- read.table(StageFile('GatherAsvs','asv2auid.abundances.tsv.gz'))[negIds,rownames(x)]
stopifnot(!is.na(nonNifH))
nonNifH <- CountDetectedASVsInSamples(nonNifH)
stopifnot( names(nonNifH) == rownames(x) )
x$endOf3 <- x$endOf4 + nonNifH  # add back non-nifH-like ASVs to get to end of 3.

## Improve the names, and set up column order before the merge.
x <- x[,c('endOf3','endOf4')]
colnames(x) <- sub('endOf4', 'ASVsFilterAuids.Length',
                   sub('endOf3', 'ASVsFilterAuids.RareAndNonNifH', colnames(x)))
## Check strictly non-increasing
x[is.na(x)] <- 0
stopifnot(sum(x[,1] < x[,2]) == 0)


## Merge
x$SampleWF <- rownames(x)
x <- merge(workflowTable, x, by='SampleWF', all=T)
x[is.na(x)] <- 0
workflowTable <- x

## Clean up
rm(filtAbunds, lenFiltTab, negIds, nonNifH, x)


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
workflowTable$ASVsWorkspaceStartup.notEmpty <- workflowTable[,ncol(workflowTable)]

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
    x <- DropSamps(workflowTable, dropme, 'ASVsWorkspaceStartup.hasSampMeta')
    workflowTable <- x

    ## Now remove samples that had no environmental metadata
    dropme <- readLines(StageFile('WorkspaceStartup','sampsWithoutEnvdata.txt'))
    x <- DropSamps(workflowTable, dropme, 'ASVsWorkspaceStartup.hasEnvMeta')
    workflowTable <- x
}


## Get for each sample the number of detected ASVs in the final abundance table.
## Because dropping of ASVs with no annotation is the final filtering step in
## WorkspaceStartup, these ASV counts minus those in the rightmost column of
## workflowTable equal the ASVs lost due to lack of annotation.
ASVsInFinalAbundanceTable <- function(colName = "final_reads")
{
    load(StageFile('WorkspaceStartup','workspace.RData'))
    stopifnot(colnames(abundTab)[1] == 'AUID')
    asvsBySample <- CountDetectedASVsInSamples(abundTab[,-1])
    df <- data.frame(Sample = names(asvsBySample), tots = as.numeric(asvsBySample), row.names = NULL)
    colnames(df) <- c('Sample', colName)
    df
}

cat("Getting ASVs in final abundance table.\n")
y <- ASVsInFinalAbundanceTable("ASVsWorkspaceStartup.noAnnot")
## The workflowTable might have samples not in the final abundance table so all.x = T.
workflowTable <- merge(workflowTable, y, by = 'Sample', all.x = T)
idx <- which(is.na(workflowTable$ASVsWorkspaceStartup.noAnnot))
workflowTable$ASVsWorkspaceStartup.noAnnot[idx] <- 0
rm(y,idx)


cat("Checking that all 'ASVs' columns in the table are numeric.\n")
idx <- grep('^ASVs',colnames(workflowTable))
stopifnot( apply(workflowTable[,idx], 2, is.numeric) )
cat("Checking for every sample that ASV counts decrease from each stage i to i+1 through the workflow.\n")
##stopifnot( apply(workflowTable[,idx], 1, function(v) all(order(v,decreasing=T)==1:length(v))) )
## Check for any cases where ASVs increase from one workflow stage to the next.
idx <- which(!apply(workflowTable[,idx], 1, function(v) all(order(v,decreasing=T)==1:length(v))))
if (length(idx) > 0) {
    cat("These have some strange increases in ASVs.\n")
    print(workflowTable[idx,])
}

cat("\nHere are the mean ASVs per sample at the end of each stage (count and %):\n")
x <- round(colMeans(workflowTable[,-c(1:4)]))
data.frame(ASVs.mean   = x,
           ASVs.median = apply(workflowTable[,-c(1:4)], 2, median),
           ASVs.sd     = round(apply(workflowTable[,-c(1:4)], 2, sd)),
           ASVs.mean.pct = round(100*x/x[1],2))
rm(x)


cat("\nWriting out table of reads retained at each stage:  workflowTable.csv\n")
write.csv(format(workflowTable, digits=2), file='workflowTable_ASVs.csv', row.names=F, quote=F)
##save.image('workflowASVsAtEachStage.RData')


##################################################
## Summary plot

cat("Preparing plot that shows for each study how the workflow progressed",
    "(reads retained at each stage).\n")

wft <- workflowTable

## Simplest if drop a few empty-from-pipes now.  'notEmpty' can still show >0 in plot.
idx <- which(wft$ASVsPipeline > 0)
cat(length(idx),"of",nrow(wft),"samples have >0 reads from the pipeline.  Using only them.\n")
wft <- wft[idx,]

## When WorkspaceStartup drops empty samples, of course that does not change the number of ASVs
## retained by any sample, as next line verifies.
stopifnot( wft$ASVsWorkspaceStartup.notEmpty - wft$ASVsFilterAuids.Length == 0 )
## Therefore, it is not helpful to show ASVsWorkspaceStartup.notEmpty in the plots.
wft <- wft[, colnames(wft) != 'ASVsWorkspaceStartup.notEmpty']
stageRename <- stageRename[names(stageRename) != 'ASVsWorkspaceStartup.notEmpty']

cat("Preparing data for ggplot...\n")
## Have ASVs at each stage.  Want the % of ASVs lost at each stage. Get it
## with diff().  Note that diff() removes the leftmost col, ASVsPipeline, which
## is good because it would always be 100%.
x <- as.matrix(wft[,rev(names(stageRename))])  # order for diff!
x <- t(apply(x, 1, function(v) { 100 * (-diff(v) / v['ASVsPipeline']) }))
dat <- cbind(wft[,c('Sample','Study')], x)
dat.m <- melt(dat, id=c('Sample','Study'))
dat.m$variable <- factor(dat.m$variable,
                         levels = names(stageRename), labels = as.character(stageRename))
colnames(dat.m) <- c('Sample','Study','Filter','value')
rm(x,dat)
## Above we verified that there are no ASVs were gained from stage i to i+1. But if there
## were any such cases, squash them in the plot to 0.
idx <- which(dat.m$value < 0)
if (length(idx) > 0) {
    cat("What?!  I thought gained ASVs were impossible and I checked, but now I see",length(idx),"gain cases.\n")
    dat.m[idx,'value'] <- 0
}


## The 'reads' script orders studies by proportion of reads contributed.  Cannot
## (easily) do similar for ASVs because AUIDs are not tracked in this script.
## Disabled code below -- it's wrong:
if (FALSE) {
tots <- aggregate(wft[,ncol(wft)], by=list(study=wft$Study), sum, na.rm=T)
colnames(tots) <- c('study','totASVs')
tots <- tots[order(tots$totASVs,decreasing=T),]
tots$pct <- round(100*tots$totASVs/ sum(tots$totASVs),1)
idx <- match(dat.m$Study, tots$study)
stopifnot(!is.na(idx))
## Bracketted numbers in facet labels show the contribution [as a %] of each study to
## the total ASVs that made it through FilterAuids (the last column in wft).
dat.m$Study <-factor(dat.m$Study,
                     levels = tots$study,
                     labels = paste0(tots$study,' [',tots$pct,'%]'))
cat("Here are the total read contributions by study:\n")
print(tots)
rm(tots,idx)
}

cat("Plotting...\n")

## The plot shows at each stage, for each data set (shown as jittered points)
## the % of _all_ ASVs from the pipeline that were discarded *AT* the indicated
## stage.  Stages happen sequentially (top-to-bottom rows in plot) so each stage
## should be interpreted as the *additional* ASVs that are lost.  (One cannot
## tell if an ASV that was lost due to being in a chimera would have later been
## lost for some other reason.)
## Very thin violins so the jitter is important.
sfb <- scale_fill_brewer(palette = 'Set3', type = 'qual')
scb <- scale_color_brewer(palette = 'Set3', type = 'qual')
g.box <- ggplot(dat.m, aes(Filter, value, fill=Filter)) +
           geom_violin() + geom_jitter(aes(color=Filter), size=0.5, alpha=0.5) +
           coord_flip() + sfb + scb +
           theme_bw() + labs(x=NULL,y='% pipeline ASVs lost at stage') +
           guides(color = 'none') +  # color legend (jitter) is redundant
           guides(fill = guide_legend(reverse = TRUE)) +
           facet_wrap(~Study)
ggsave(filename='workflowASVsRetainedEachStage.png', plot=g.box, width=11.5, height=7, units='in', dpi=144)

## Also make an overall plot (all samples across all studies).  Drop legend and
## labels since will manually layout to use those from plot above.
g.box <- ggplot(dat.m, aes(Filter, value, fill=Filter)) +
           geom_violin() + coord_flip() + sfb +
           theme_bw() + theme(axis.text.y = element_blank()) +
           labs(x = element_blank(), y = element_blank()) +
           guides(fill = 'none')   # Use legend of main plot
ggsave(filename='workflowASVsRetainedEachStage_overall_violin.png', plot=g.box,
       width=11.5/6, height=7/4, units='in', dpi=144)

## Box plot might look better.
g.box <- ggplot(dat.m, aes(Filter, value, fill=Filter)) +
           geom_boxplot(outlier.size=0.5, outlier.alpha=0.5) + coord_flip() + sfb +
           theme_bw() + theme(axis.text.y = element_blank()) +
           labs(x = element_blank(), y = element_blank()) +
           guides(fill = 'none')   # Use legend of main plot
ggsave(filename='workflowASVsRetainedEachStage_overall_box.png', plot=g.box,
       width=11.5/6, height=7/4, units='in', dpi=144)

cat("Done. See workflowASVsRetainedEachStage.png and the 'overall' versions.\n",
    "Note that stage ASVsWorkspaceStartup.notEmpty is not shown because it discards\n",
    "empty samples and does not impact the read counts.\n")
quit(save='no')
