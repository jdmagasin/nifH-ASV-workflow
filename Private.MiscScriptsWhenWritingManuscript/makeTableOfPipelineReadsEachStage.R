#!/usr/bin/env Rscript

##
## Make tables that shows how many reads are retained through each stage of the DADA2 pipeline: by
## sample, by study.  Also a multi-box plot faceted by study.
##
## Adapted from version in ../../../MiscScriptsWhenWritingManuscript/
##

options(width=120)
library(ggplot2)
library(reshape2)

## Take the workflow dir so that we can load GatherASVs's main tsv, which points to
## the pipeline files we need (readCountsDuringProcessing*.csv).
workflowRootDir = commandArgs(T)[1]
cat("Will use pipeline outputs rooted at",workflowRootDir,"\n")
stopifnot(dir.exists(workflowRootDir))

gatherTabTsv <- file.path(workflowRootDir, "GatherAsvs", "asvs.noChimera.fasta_table.tsv")
cat("Loading", gatherTabTsv, "to get to the pipeline reads at each stage.\n")
gatherTab <- read.table(gatherTabTsv)
colnames(gatherTab) <- c('Study', 'ASV_noChimera_FASTA', 'tag')
## Format of the tag is <studyID>:<processingGroupNum>:<SRA study>

gatherTab$ASV_noChimera_FASTA <- sapply(gatherTab$ASV_noChimera_FASTA, function (fas) {
    if (startsWith(fas, '../')) {
        ## 'fas' is relative to GatherAsvs. Make it relative to the workflowRootDir.
        fas <- file.path(workflowRootDir, sub('^\\.\\./', '', fas))
    }
    fas
})

if (FALSE) {
    ## DISABLED because barely used and I want to focus on the pipeline.

    ## Load workspace to get the abundance table so that we can add a column for the reads retained
    ## *after* the pipeline, e.g. after uchime3_denovo and any other post-pipeline removal of ASVs or
    ## entire samples.  The sample names in smetaTab don't easily match up to the Sample names that end
    ## up in pipe.samps so currently the column is NOT added.
    cat("Loading the workspace created by WorkspaceStartup.\n")
    load(file.path(workflowRootDir, "WorkspaceStartup", "workspace.RData"))
    sampsInAbundTab <- sub('^X','',colnames(abundTab))
}

## For each FASTA, get the pipeline stages file which is 2 dirs up from each asvs.noChimera.fasta.
## Except for the data sets with mixed-orientation seq libs, which need special handling.
cat("\nGathering the readCountsDuringProcessing* csv files for every pipeline run\n",
    "(processing group), except not for those that had mixed-orientation seq libs.\n")
readsOverStages <- sapply(gatherTab$ASV_noChimera_FASTA, function (fas) {
    tgtPat <- 'readCountsDuringProcessing.*.csv$'
    ## Data sets with mixed-orientation seq libs have different dir structure. For
    ## all others, can back up 2 dirs to find the readCountsDuringProcessing csvs.
    mixLibs <- paste0('Shiozaki_',c('2017','2018GBC'))
    x <- dirname(dirname(fas))
    if (grepl('\\/Dada2PipeOutput\\/',x)) {  # Normal data set
        x <- list.files(x, tgtPat, full.names=T)
        stopifnot(length(x) == 1)
    } else if (any(sapply(mixLibs, grepl, x=fas))) {
        ## Hack.  Works for the two Shiozaki m-o data sets.  This finds both the Normal and Swapped
        ## csvs.  For simplicity use just the Normal.  Means counts will be ~halved for these data
        ## sets -- note in figure caption (see reminder at end of this script).
        x <- dirname(fas)
        x <- list.files(x, tgtPat, full.names=T, recursive=T)
        want <- c('RNA','DNA')[grepl('combined\\.DNA',fas)+1]
        x <- x[grep(want,    x)]  # Get the D/RNA as specified by the gatherTab entry
        x <- x[grep('Normal',x)]
        stopifnot(length(x) == 1)
    } else { x <- NULL }
    x
})
cat("Found",length(readsOverStages),"readCountsDuringProcessing* csv's.\n")
## Next check assumes the loop gives just 1 csv for the mixed-orientation cases.
stopifnot(length(readsOverStages) == nrow(gatherTab))
## Nothing should be missing at this near-submission stage!
stopifnot(!sapply(readsOverStages,is.null))

## Load each reads-over-stages csv.  Then bind into one big data frame.  Note comments in
## ProcessOneReadCountCsv() about dropping the R2 counts so that the final data frame has counts for
## read *pairs*.
cat("Loading the reads at each stage for every pipeline run.\n")
stopifnot( length(unique(gatherTab$tag)) >= length(unique(readsOverStages)) )
names(readsOverStages) <- gatherTab$tag

## Helper to make a data frame for one processing group (readCountsDuringProcessing csv).
ProcessOneReadCountCsv <- function (csv, tag)
{
    d <- read.csv(csv)
    ## Each sample has R1 and R2 have counts in the Initial, Trimmed, and Filtered stages, but the
    ## counts ~halve at/after ASV creation.  Drop the R2 counts, and keep/treat the R1 counts as if
    ## they are for read *pairs*.  The (_|$) is critical or some sets will not have their R2's
    ## dropped and it will appear as if the Merged stage throws out ~half the reads ==> Jonathan
    ## wondering if the pipeline parameters need reevaluation, gahr!!
    d <- d[grep('_R2(_|$)', d$Sample, invert=T),]  
    
    ## The "Sample" names change from stages that use fastqs to later stages.  Define "root" sample
    ## names based on the Merged stage names. The roots lack _R1_, or _L001_R1_, or
    ## <station>_L001_R1, etc. Unfortunately the ways names changed pre/post-Merged isn't the same
    ## for all studies so below we match the 'snams' to the pre-Merged Sample names. Yuck.
    snams <- unique(subset(d, Stage=='Merged')$Sample)
    stopifnot(sum(d$Stage=='Initial') == length(snams))
    samp2root <- lapply(snams, function(s) grep(paste0('^',s), d$Sample))
    names(samp2root) <- snams
    ## Check that no 2 roots hit the same sample.
    stopifnot( table(unlist(samp2root)) == 1 )
    for (i in 1:length(samp2root)) { d[samp2root[[i]],'SampRoot'] <- names(samp2root)[i] }
    ## Verify that every sample has 6 stages.
    stopifnot( aggregate(Stage ~ SampRoot, d, length)$Stage == 6 )
    d$tag <- tag
    d$Study <- sub(':.*','',tag)
    d <- d[,colnames(d) != 'Sample']
    colnames(d) <- sub('SampRoot','Sample',colnames(d))
    ## Debug:  Can check each sample's read flow:  aggregate(Reads~Stage+Sample, d, sum)
    d
}

## Bind into one big data frame of per stage stats for read *pairs*.
## (x is a temp, see pipe.samps below.)
x <- lapply(names(readsOverStages), function(tag) {
    d <- NULL
    rosCsv <- readsOverStages[[tag]]
    if (length(rosCsv) > 0) {  # Shiozaki studies 0 due to mixed-orientation
        ##cat("Working on",rosCsv,"\n")
        d <- ProcessOneReadCountCsv(rosCsv,tag)
    }
    d
})
names(x) <- NULL
x <- do.call(rbind, x)
if (F) { ## Debug: 
    y <- aggregate(Reads ~ Stage + Study, x, sum)
    y[order(y$Study,y$Reads,decreasing=T),]
    rm(y)
}


cat("Reformatting the data into a table of sample read pair counts at each stage.\n")
## 'x' is set up for ggplot (just like the csvs when the pipeline saved them).  Cast it into data
## frame that has one row per sample and columns for each pipeline stage.
pipe.samps <- dcast(x, Study + Sample + tag ~ Stage, value.var='Reads', sum)
## Exclude InASVs because it is the same as Merged based on comment in
## countReadsAtStages.R (pipeline script) as well as checks in this script below.
x <- c('Study', 'Sample', 'tag', 'Initial', 'Trimmed', 'Filtered', 'Merged', 'InNonChimericASVs')
pipe.samps <- pipe.samps[,x]

## Verify that read pair counts strictly decrease
if (T) {
    x <- c('Initial', 'Trimmed', 'Filtered', 'Merged', 'InNonChimericASVs')
    x <- t(apply(pipe.samps[,x], 1, order, decreasing=T))
    stopifnot( unique(x) == 1:5 ) # 1..5 is the only pattern
}

pipe.samps$PctReadPairsRetained <- 100*(pipe.samps$InNonChimericASVs / pipe.samps$Initial)
x <- subset(pipe.samps, InNonChimericASVs > 0)

if (FALSE) { # Not super helpful
    cat("Here is a stem plot of the % read pairs retained (initial to non-bimera ASVs) for the\n",
        nrow(x),"samples run through the DADA2 pipeline with >0 reads in ASVs at the end.\n")
    print(stem(x$PctReadPairsRetained))
}

cat("\nStats on the % of paired reads retained (when *pipeline* finishes).\n",
    "Only includes",nrow(x),"samples that had >0 non-bimera ASVs.\n")
round(summary(x$PctReadPairsRetained), 3)

cat("\nStats on the number of paired reads InNonChimericASVs (i.e. pipeline finish):\n")
summary(x$InNonChimericASVs)

if (FALSE) {
    ## DISABLED because barely used and I want to focus on the pipeline.
    
    cat("\nFor comparison, stats on the counts in the abundance table. These are post-pipeline\n",
        "and include additional filtering (uchime3_denovo, etc.).  There are", ncol(abundTab)-1,
        "\nsamples in the abundance table:\n")
    summary(colSums(abundTab[,-1]))
    cat("\n")
}

fnam <- "Table_SreadsAtEachStage_samples.csv"
cat("Writing",fnam,"\n")  # with 2 signif digits, impacts PctReadPairsRetained
x <- pipe.samps[,colnames(pipe.samps) != 'tag']
write.csv(format(x, digits=2), file=fnam, row.names=F, quote=F)


# # # # # # # # # # # #
# Plot!
#

tab <- pipe.samps
##stopifnot( tab$Merged == tab$InASVs )           # Dropped InASVS (see above)
if (FALSE) {
    cat("Discarding",sum(tab$Initial < 1000),"samples with <1K initial reads since",
        "post-pipeline they would have\nhad <1K reads in ASVs and so would have been dropped.\n")
    tab <- subset(tab, Initial >= 1000)
}
x <- c('Study','Sample','Initial','Trimmed','Filtered','Merged','InNonChimericASVs')
x <- tab[,x[-c(1:2)]]                             # Just the numeric cols
x <- scale(t(x), center=F, scale = x$Initial)     # Proportions relative to Initial
x <- as.data.frame(t(x))
stopifnot( max(abs(x$InNonChimericASVs - tab$PctReadPairsRetained/100)) < 0.001 )  # check scale()
tab.propKept <- cbind(tab[,c('Study','Sample')], x)
stopifnot( tab.propKept$Initial == 1 )
## Replace Initial with the initial paired read counts (R1 as approx), and tack
## on the final read count.
tab.propKept$Initial <- tab$Initial
tab.propKept$Final_count <- tab$InNonChimericASVs
stages <- c('Initial_count','Trimmed','Filtered','Merged','InNonChimericASVs','Final_count')
colnames(tab.propKept) <- c('Study','Sample',stages)


## Another useful table which shows the each study's average proportion of reads
## retained at each stage, and initial and final read counts.
fnam <- "Table_SreadsAtEachStage_studies_medProportions.csv"
cat("\nAlso writing out",fnam,"\n",
    "This table has the median (over samples in each study) proportions of reads\n",
    "retained at each stage of the DADA2 pipeline. Median Initial and Final counts are\n",
    "also provided.\n")
x <- tab.propKept[,colnames(tab.propKept) != 'Sample']
nums <- aggregate(Initial_count ~ Study, x, length); colnames(nums) <- c('Study','NumSamps')
x <- aggregate(. ~ Study, x, median)
x <- merge(x, nums)
x <- x[,c('Study','NumSamps',setdiff(colnames(x),c('Study','NumSamps')))]
write.csv(format(x, digits=2), file=fnam, row.names=F, quote=F)

fnam <- "Table_SreadsAtEachStage_studies_counts.csv"
cat("\nAlso writing out",fnam,"\n",
    "which has the paired read counts at each stage for each study. (Not medians!)\n")
## Drop 'tag' and 'Sample'.
x <- pipe.samps[,c('Study',names(which(sapply(pipe.samps, is.numeric))))]
x <- aggregate( . ~ Study, x, sum)  # PctReadPairsRetained wrong but updated next
x$PctReadPairsRetained <- round(100*x[,'InNonChimericASVs']/x[,'Initial'],2)
x <- merge(x, nums)
x <- x[,c('Study','NumSamps',setdiff(colnames(x),c('Study','NumSamps')))]
write.csv(x, file=fnam, row.names=F, quote=F)


##
## Box plot panel (by study)
##
x <- c('Study','Sample',setdiff(stages,c('Initial_count','Final_count')))
dat <- melt(tab.propKept[,x], id=c('Study','Sample'))
colnames(dat) <- c('Study','Sample','Step','PropKept')
dat$Step  <- factor(dat$Step, levels = c("Trimmed","Filtered","Merged","InNonChimericASVs"),
                              labels = c("Trimmed [1]", "Filtered [4]", "Merged [7]", "Not bimera [9]"))
dat$Study <- factor(dat$Study, levels = unique(sort(dat$Study)))  # order studies (facets) alphabetically

g <- ggplot(dat, aes(x=Step, y=PropKept, group=Step, fill=Step)) +
     geom_boxplot() + facet_wrap(~Study) +
     theme_bw() +
     ##theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
     theme(axis.text.x = element_blank(), strip.text.y = element_text(size = 8, )) +  ## facet lab font sz
     labs(x="DADA2 pipeline step", y="Proportion reads retained")
ggsave('pipeReadsRetainedEachStep.png', plot=g, width=8.3, height=7, units='in', dpi=144)
cat("\nSaved pipeReadsRetainedEachStep.png\n")


## Also plot for all the samples together. Just drop the facets.  Make the plot
## the ~same size as each of the 5x4 facets so that we can manually place it on
## the right of the final Fig 3, in the same column as the legend.
g <- ggplot(dat, aes(x=Step, y=PropKept, group=Step, fill=Step)) +
     geom_boxplot(show.legend = FALSE) + 
     theme_bw() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
     labs(x = element_blank(), y = element_blank()) +
     guides(fill = 'none')  # Instead use legend of main plot
ggsave('pipeReadsRetainedEachStep_overall.png', plot=g,
       width = 8.3/6, height = 7/4, units='in', dpi=144)  # Should make plot similar size as facets
cat("\nSaved pipeReadsRetainedEachStep_overall.png\n")

cat("NOTE!  The studies", paste0('Shiozaki_',c('2017','2018GBC')), "are mixed-orientation.\n",
    "The stats reflect just the Normal pipepline runs, and since ~half of the reads in were\n",
    "for the wrong strand and thus rejected, the proportions of reads will be < 50%.\n")
##save.image('image.Rdata')
quit(save='no')