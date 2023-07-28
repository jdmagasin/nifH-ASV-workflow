#!/bin/env Rscript

##
## Given ASVs from the Normal and Swapped-primer runs of the pipeline, check the
## extent of ASVs shared and combine ASVs to produce one abundance table and
## FASTA.
##
## This script assumes there is one processing group with Normal and Swapped ASV
## tables that need to be combined.  If there are >1 processing groups, instead
## use multiCombineNormalAndSwapped.sh.
##

indexOfAsvTables <- as.numeric(commandArgs(T)[1])  ## See "HACK" below.
stopifnot(indexOfAsvTables > 0)
rundirs <- c('Normal','Swapped') # safer
stopifnot(sapply(rundirs, dir.exists))
names(rundirs) <- c('normal','swapped')

cat("Loading libraries...")
suppressMessages(library(ShortRead))
suppressMessages(library(vegan))
cat("\n")


GetAsvFastas <- function(dir, revcomp=FALSE) {
    fl <- sort(list.files(dir, "asvs.noChimera.fasta", full.names=T, recursive=T))
    fasl <- lapply(fl, readFasta)
    names(fasl) <- fl
    cat(paste0("Loaded ",length(fasl)," ASVs from fasta in ",dir,".\n"))
    if (revcomp) {
        cat("Reverse-complementing the",dir,"ASVs.\n")
        fasl <- lapply(fasl,
                       function(fas) {
                           ShortRead(id = id(fas),
                                     sread = reverseComplement(sread(fas)))
                       })
    }
    fasl
}

LoadAbundanceTables <- function(dir) {
    fl <- sort(list.files(dir, "asvs.noChimera.tsv", full.names=T, recursive=T))
    tabl <- lapply(fl, read.table)
    names(tabl) <- fl
    lapply(names(tabl), function (nam) {
        cat(nam, "has",nrow(tabl[[nam]]), "ASVs in", ncol(tabl[[nam]]), "samples",
            "and",sum(tabl[[nam]]),"total reads.\n")
    })
    tabl
}


asvs.normal  <- GetAsvFastas(rundirs['normal'])
asvs.swapped <- GetAsvFastas(rundirs['swapped'], TRUE)
## Better correspond
if (! all(sub('^Normal', '', dirname(names(asvs.normal))) ==
          sub('^Swapped','', dirname(names(asvs.swapped)))) ) {
    stop("The names of the Normal and Swapped ASV tables do",
         "not have identical paths except for 'Normal' or",
         "'Swapped' at the start. So I do not know which",
         "normal and swapped ASV tables correspond.")
} else {
    cat("\nGood Normal and Swapped ASV tables appear to correspond\n",
        "base on their file paths.\n")
}
atab.normal  <- LoadAbundanceTables(rundirs['normal'])
atab.swapped <- LoadAbundanceTables(rundirs['swapped'])


##
## HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
##
## This script was written assuming one normal and one swapped ASV table to
## combine -- except for the code above which loads all ASV tables.  I don't
## want to rewrite everything below to make it handle corresponding lists of ASV
## tables.  Instead use the index specified on the command line.
cat(">>>>>>> Will only combine the",indexOfAsvTables,"'th ASV tables. <<<<<<<\n")
stopifnot(indexOfAsvTables >= 1 && indexOfAsvTables <= length(atab.normal))
asvs.normal  <- asvs.normal[[indexOfAsvTables]]
asvs.swapped <- asvs.swapped[[indexOfAsvTables]]
atab.normal  <- atab.normal[[indexOfAsvTables]]
atab.swapped <- atab.swapped[[indexOfAsvTables]]


## Rename ASV ids as sequences so that we can compare abundances (in the same
## namespace).
stopifnot( match(rownames(atab.normal),  id(asvs.normal))  == 1:nrow(atab.normal) )
stopifnot( match(rownames(atab.swapped), id(asvs.swapped)) == 1:nrow(atab.swapped) )
rownames(atab.normal)  <- as.character(sread(asvs.normal))
rownames(atab.swapped) <- as.character(sread(asvs.swapped))


## What ASVs are identical between the two runs?
sharedAsvs <- intersect(sread(asvs.normal), sread(asvs.swapped))
cat("There are",length(sharedAsvs),"identical ASVs in Normal and Swapped.\n")
x <- as.character(sharedAsvs)
cat("and they account for",
    paste0(round(100*sum(rowSums(atab.normal[x,])) / sum(atab.normal), 1),"%"),
    "of total Normal reads and",
    paste0(round(100*sum(rowSums(atab.swapped[x,])) / sum(atab.swapped), 1),"%"),
    "of total Swapped reads.\n")
rm(x)


## Expand tables to have the same ASVs (rows) for later rbind'ing.
m.n <- atab.normal
m.s <- atab.swapped
allSeqs <- union(rownames(m.n),rownames(m.s))
m.n[setdiff(allSeqs,rownames(m.n)),] <- 0
m.s[setdiff(allSeqs,rownames(m.s)),] <- 0
m.n <- m.n[allSeqs,]
m.s <- m.s[allSeqs,]



x <- decostand(t(m.n), 'total')            # 'total' so normalize samples (rows)
rownames(x) <- paste0(rownames(x),'.n')
y <- decostand(t(m.s), 'total')
rownames(y) <- paste0(rownames(y),'.s')
dat <- rbind(x,y)
rm(x,y)


## Dendrogram showing whether Normal and Swapped sample pairs cluster.
cat("Hierarchically clustering to see whether Normal and Swapped\n",
    "abundance profiles for samples pair up. See cluster_normal_vs_swapped.png\n")
hc <- hclust(vegdist(dat,'bray'), 'centroid')
sink('/dev/null')
png(paste0('cluster_normal_vs_swapped.',indexOfAsvTables,'.png'),
     width=7, height=7, units='in', bg='white', res=144)
plot(hc, main='Samples rep\'d by Normal vs. Swapped ASVs',
     sub='Centroid clustering of Bray-Curtis\nsample dissimilarities', xlab='',
     cex=0.7)
dev.off()
sink()


## Make abundance table and FASTA.  Simply add counts from the two runs.
## Each run used ~half the reads.
stopifnot(rownames(m.n)==rownames(m.s))  # ASVs match up
stopifnot(colnames(m.n)==colnames(m.s))  # samples match up
dat <- m.n + m.s
## Now create new IDs for the ASVs that reflect the abundance in Normal and
## Swapped (simply added).
asvTots <- sort(rowSums(dat), decreasing=T)
asvIds <- paste0('ASV.0',1:length(asvTots))  # "0" prefix to signal that not from pipeline
rownames(dat) <- asvIds[match(rownames(dat), names(asvTots))]
dat <- dat[asvIds,]
fn <- paste0("asvs.noChimera.combined.",indexOfAsvTables,".tsv")
write.table(dat, fn, sep="\t", quote=F)
sr <- ShortRead(id = BStringSet(asvIds), sread = DNAStringSet(names(asvTots)))
fn <- paste0("asvs.noChimera.combined.",indexOfAsvTables,".fasta")
writeFasta(sr, fn)
cat("Created", paste0("asvs.noChimera.combined.",indexOfAsvTables,".{tsv,fasta}."),"\n")
##cat("You should rename outputs to reflect the processing group since '",
##    indexOfAsvTables, "' is not super informative.\n")

quit('no')
