#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin
##
## Make a big table of ASV annotation with respect to several references:
##  - ARB 2017:  Zehr Lab nifH ARB database as of 2017.
##  - Genomes879: Genomes from 879 diazotrophs
##  - Marine diazotrophs nifH DB: Small DB with cyano and noncyano diazotrophs, some assembled from
##    metagenomic data.
##  - CART:  NifH cluster classifier by Frank et 2016
##
## The first 3 tables are from BLAST searches; take the best hit (by E-value). For ARB2017 hits that
## lack useful annotation, try to find an equally good and informative hit.
##  -- Unfortunately this didn't help.  Rather than 'equally good' perhaps try hits with slightly
##     worse E-values...
##
usageStr <- "
makeAnnotationTable.R  <ARB2017 .tab>  <genome 879 .tab>  <marine diazo .tab>  <CART clusters .map>

Input:  Four annotation files listed above. The first three inputs are tabular BLAST results (tab-
        separated columns).  The CART input is a '.map' file but uses commas to separate columns.

Output: auids.annot.tsv   Table with annotation based on BLAST searches of ARB2017, Genomes879, and
                          the marine diazotroph DB, including %id and E-values.  NifH classification
                          by CART is also included.  NA appears if no annotation was found (at the
                          cut offs you used). For the ARB2017 and diazotroph DB hits, taxonomy
                          information is included.
"

##
## Get annotation file names
##
annotFiles <- commandArgs(T)[c(1:4)]
idx <- which(sapply(annotFiles, file.exists))
if (length(idx) != 4) {
    cat("\nError: Missing input file(s)",annotFiles[-idx],"\n")
    cat(usageStr)
    stop("Aborting.")
}
names(annotFiles) <- c('ARB2017','Genomes879','MarineDiazo','CART')


##
## Load results from annotation tools
##
cat("Loading annotation tables\n")
annot <- lapply(names(annotFiles), function(x) {
    fn <- annotFiles[x]
    cat("Loading",fn,"...")
    if (grepl('\\.map$',fn)) { x <- read.csv(fn) } else { x <- read.table(fn) }
    cat("okay, loaded", paste(dim(x),collapse=' by '),"table.\n")
    x
})
names(annot) <- names(annotFiles)

blastHeader <- c('AUID', 'subject','pid','alen','mismatch','gapopen',
                 'q.start','q.end','s.start','s.end','evalue','bits')
colnames(annot$ARB2017)     <- blastHeader
colnames(annot$Genomes879)  <- blastHeader
colnames(annot$MarineDiazo) <- blastHeader
colnames(annot$CART)        <- c('AUID', 'cluster','subcluster')

cat("Loading Genomes879 taxomomy map...")
g879taxPath <- file.path('./data','genome879_acc_taxstring.txt')
g879taxa <- read.table(g879taxPath, header=T, sep="\t")
## Next line shows that NA appears for no field in the G879 table, except 'kingdom'
## apply(g879taxa,2,function(x) sum(is.na(x)))
cat("done.\n")

cat("Loading marine cyano and NCD sequence information...")
cyanoNCDtaxPath <- file.path('./data','Marine_NCD_cyano_nifH',
                             'Marine_NCD_cyano_refsequences.txt')
cyanoNCDtaxa <- read.table(cyanoNCDtaxPath, header=T, sep="\t")
cat("done.\n")

bestHits <- list()

##
## Get top hits to ARB2017. Then if necessary replace with equally top but
## better annotated hits.
##
## match() gets the top b/c blast results are sorte by E-value for each query.
auids <- unique(annot$ARB2017$AUID)
bestHitsIdx <- match(auids, annot$ARB2017$AUID)

taxa <- sub('^.+\\|','',annot$ARB2017[bestHitsIdx,'subject'])  # strip accessions
useless <- c('^uncultured_marine','^uncultured_bacterium','^uncultured_nitrogen-fixing',
             '^uncultured_microorganism','^uncultured_soil','^metagenome')
idx <- unique(unlist(lapply(useless, function(pat) which(regexpr(pat,taxa)!=-1))))

## For all the best hits with useless annotation, try to find an equally good hit
## that is not useless.  This code is a little complicated.  The function takes
## an index into annot$ARB2017, uses it to create an index of the next 50 hits,
## which it searches: must be from the target AUID, have same E-value, and a
## subject that doesn't start with "uncultured". Function returns an index into
## annot$ARB2017.
alts <- sapply(bestHitsIdx[idx],
       function(i) {
           ## Get this AUID's top E-value.
           auid <- annot$ARB2017[i,'auid']
           eval <- annot$ARB2017[i,'evalue']
           ## Consider the next 50 hits.  Little risk of this running off table's end.
           altsIdx <- seq(i+1,i+51)
           altsIdx <- altsIdx[ which((annot$ARB2017[altsIdx,'AUID'] == auid) &
                                     (annot$ARB2017[altsIdx,'evalue'] <= eval)) ]
           if (length(altsIdx) > 0) {
               altsIdx <- altsIdx[ grep('^uncultured',annot$ARB2017[altsIdx,'subject'],invert=T) ]
           }
           if (length(altsIdx) > 0) { i <- altsIdx[1] } # Found a replacement(s)
           i
       })

cat("There were",length(idx),"AUIDs that had top hits to ARB that had unhelpful annotation (e.g. uncultured_microorganism).\n",
    "Tried to find hits with equally good E-value but helpful annotation.  Was able to improve the annotation for",
    sum(bestHitsIdx[idx] != alts),"AUIDs.\n")
bestHitsIdx[idx] <- alts
stopifnot(length(bestHitsIdx) == length(auids))

bestHits$ARB2017 <- annot$ARB2017[bestHitsIdx, c('AUID','subject','pid','alen','evalue')]


##
## Genomes879 don't need alternatives -- we always know the gennome we hit:)
##
auids <- unique(annot$Genomes879$AUID)
bestHitsIdx <- match(auids, annot$Genomes879$AUID)
bestHits$Genomes879 <- annot$Genomes879[bestHitsIdx,c('AUID','subject','pid','alen','evalue')]

## Merge the blast annotations
df <- merge(x = bestHits$ARB2017, y = bestHits$Genomes879, by='AUID', all=T)
colnames(df) <- c('AUID',
                  'ARB2017.id',   'ARB2017.pctId',   'ARB2017.alen',   'ARB2017.evalue',
                  'Genomes879.id','Genomes879.pctId','Genomes879.alen','Genomes879.evalue')
stopifnot(setequal(df$AUID, c(bestHits$ARB2017$AUID, bestHits$Genomes879$AUID)))

## For top Genomes879 hits, replace the id (which is an ARB name apparently)
## with the genus_species.  Also tack on taxonomic info.
cat("Looking up more informative names and taxa for the hits to Genomes879.\n")
ghits <- as.character(df$Genomes879.id)
idx <- match(ghits, g879taxa$arb_name)
ghits[!is.na(idx)] <- g879taxa[idx[!is.na(idx)],'genus_species']
stopifnot(!g879taxa$genus_species %in% c(NA,'')) # So we did not make ghits worse.
df$Genomes879.id <- ghits
## Now the taxa, in a single column. Separate levels with ';'.  Use 'unknown' for
## levels that are NA (only kingdom).
x <- apply(g879taxa[idx[!is.na(idx)],c('kingdom','phylum','class','order','family','genus')],1,
          function(x) { x[is.na(x)] <- 'unknown';  paste(x,collapse=';') })
df[!is.na(idx),'Genomes879.tax'] <- x

## Genomes879 has NifH cluster or subcluster.  Here we could temporarily tack on
## the cluster and then compare it to what CART assigned.  But I don't know if
## the subcluster for the best hit in G879 should be the same as the CART
## assigned subcluster.  Often I see mismatches.
## For *cluster* mismatches are rare e.g. 63 mismatches between CART vs the best
## G879 hit, but 3707 matches.  The mismatches are often ~catch-all cluster "4"
## from CART but cluster 3 for the G879 hit.
## df[!is.na(idx),'TEMPCLUST'] <- sub('','',g879taxa[idx[!is.na(idx)],'nifHCluster'])


##
## Marine cyano/NCD nifH DB
##
auids <- unique(annot$MarineDiazo$AUID)
bestHitsIdx <- match(auids, annot$MarineDiazo$AUID)
bestHits$MarineDiazo <- annot$MarineDiazo[bestHitsIdx,c('AUID','subject','pid','alen','evalue')]

## Tack on 'description' from the annotation companion to the FASTA
mdat <- bestHits$MarineDiazo
idx <- match(mdat$subject, cyanoNCDtaxa$seqID)
stopifnot(!is.na(idx))  # If trips, then the FASTA and annotation are out of sync.
mdat$description <- cyanoNCDtaxa[idx,'description']

## Merge with the ARB and G879 blast annotations
colnames(mdat) <- c('AUID', paste0('MarineDiazo.', colnames(mdat)[-1]))
df <- merge(x = df, y = mdat, by='AUID', all=T)
stopifnot(setequal(df$AUID, c(bestHits$ARB2017$AUID, bestHits$Genomes879$AUID,
                              bestHits$MarineDiazo$AUID)))

##
## And now add CART annotation.
##
df <- merge(x = df, y = annot$CART, by='AUID', all=T)
stopifnot(setequal(df$AUID, c(bestHits$ARB2017$AUID, bestHits$Genomes879$AUID,
                              bestHits$MarineDiazo$AUID, annot$CART$AUID)))

ofile <- "auids.annot.tsv"
cat("Writing out the big annotation table",ofile,"\n")
write.table(df,ofile,quote=F,sep="\t",row.names=F)
