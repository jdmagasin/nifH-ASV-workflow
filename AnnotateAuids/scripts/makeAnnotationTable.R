#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin
##
## Make a big table of ASV annotation with respect to several references:
##  1- ARB 2017:  Zehr Lab nifH ARB database as of 2017.
##  2- Genome879: Genomes from 879 diazotrophs
##  3- Marine diazotrophs nifH DB: Small DB with cyano and noncyano diazotrophs, some assembled from
##     metagenomic data.
##  4- UCYN-A oligos: A small DB of UCYN-A oligotypes
##  5- CART:  NifH cluster classifier by Frank et 2016
##  6- Conserved cysteines check:  Results from a check for C's usually found in the highly conserved
##     Switch I and II regions of NifH.  AMP just after the Switch II C is also part of the checked
##     patterns.
##
## The first 4 tables are from BLAST searches; take the best hit (by E-value). For ARB2017 hits that
## lack useful annotation, try to find an equally good and informative hit.
##  -- Unfortunately this didn't help.  Rather than 'equally good' perhaps try hits with slightly
##     worse E-values...
##
usageStr <- "
makeAnnotationTable.R  <ARB2017 .tab>  <genome 879 .tab>  <marine diazo .tab>  <ucyna_oligos .tab>
                       <CART clusters .map>  <CCAMP check .csv>

Input:  Six annotation files listed above. The first four inputs are tabular BLAST results (tab-
        separated columns).  The CART input is a '.map' file but uses commas to separate columns.
        The CCAMP input is a table from check_CCAMP.R.

Output: auids.annot.raw.tsv   Table with annotation based on BLAST searches of ARB2017, Genome879,
                              and the marine diazotroph DB, including %id and E-values.  NifH
                              classification by CART is also included, as well as results from the
                              conserved cysteines check.  NA appears if no annotation was found (at
                              the cut offs you used). For the ARB2017 and diazotroph DB hits,
                              taxonomy information is included.
"

##
## Get annotation file names
##
annotFiles <- commandArgs(T)[c(1:6)]
idx <- which(sapply(annotFiles, file.exists))
if (length(idx) != 6) {
    cat("\nError: Missing input file(s)", annotFiles[-idx], "\n")
    cat(usageStr)
    stop("Aborting.")
}
names(annotFiles) <- c("ARB2017", "Genome879", "MarineDiazo", "UCYNAoligos", "CART", "CCAMP")

##
## Load results from annotation tools
##
cat("Loading annotation tables\n")
annot <- lapply(names(annotFiles), function(x) {
    fn <- annotFiles[x]
    cat("Loading", fn, "...")
    if (grepl("\\.(map|csv)$", fn)) {
        x <- read.csv(fn)
    } else {
        x <- read.table(fn)
    }
    cat("okay, loaded", paste(dim(x), collapse = " by "), "table.\n")
    x
})
names(annot) <- names(annotFiles)

blastHeader <- c(
    "AUID", "subject", "pid", "alen", "mismatch", "gapopen",
    "q.start", "q.end", "s.start", "s.end", "evalue", "bits"
)
colnames(annot$ARB2017) <- blastHeader
colnames(annot$Genome879) <- blastHeader
colnames(annot$MarineDiazo) <- blastHeader
colnames(annot$UCYNAoligos) <- blastHeader
colnames(annot$CART) <- c("AUID", "cluster", "subcluster")
colnames(annot$CCAMP) <- c("AUID", "orf.len", "hasCC.len", "hasCCAMP.len")

cat("Loading Genome879 taxomomy map...")
g879taxPath <- file.path("./data", "genome879_acc_taxstring.txt")
g879taxa <- read.table(g879taxPath, header = T, sep = "\t")
## Next line shows that NA appears for no field in the G879 table, except 'kingdom'
## apply(g879taxa,2,function(x) sum(is.na(x)))
cat("done.\n")

cat("Loading marine cyano and NCD sequence information...")
cyanoNCDtaxPath <- file.path(
    "./data", "Marine_NCD_cyano_nifH",
    "Marine_NCD_cyano_refsequences.txt"
)
cyanoNCDtaxa <- read.table(cyanoNCDtaxPath, header = T, sep = "\t")
cat("done.\n")


cat("Loading UCYN-A oligos sequence information...")
ucynaOLIGOtaxPath <- file.path(
    "./data", "UCYNA_oligoreps",
    "UCYNA_oligoreps.txt"
)
ucynaOLIGOtaxa <- read.table(ucynaOLIGOtaxPath, header = T, sep = "\t")
cat("done.\n")


bestHits <- list()

##
## Get top hits to ARB2017. Then if necessary replace with equally top but
## better annotated hits.
##
## match() gets the top b/c blast results are sorte by E-value for each query.
auids <- unique(annot$ARB2017$AUID)
bestHitsIdx <- match(auids, annot$ARB2017$AUID)

taxa <- sub("^.+\\|", "", annot$ARB2017[bestHitsIdx, "subject"]) # strip accessions
useless <- c(
    "^uncultured_marine", "^uncultured_bacterium", "^uncultured_nitrogen-fixing",
    "^uncultured_microorganism", "^uncultured_soil", "^metagenome"
)
idx <- unique(unlist(lapply(useless, function(pat) which(regexpr(pat, taxa) != -1))))

## For all the best hits with useless annotation, try to find an equally good hit
## that is not useless.  This code is a little complicated.  The function takes
## an index into annot$ARB2017, uses it to create an index of the next 50 hits,
## which it searches: must be from the target AUID, have same E-value, and a
## subject that doesn't start with "uncultured". Function returns an index into
## annot$ARB2017.
alts <- sapply(
    bestHitsIdx[idx],
    function(i) {
        ## Get this AUID's top E-value.
        auid <- annot$ARB2017[i, "auid"]
        eval <- annot$ARB2017[i, "evalue"]
        ## Consider the next 50 hits.  Little risk of this running off table's end.
        altsIdx <- seq(i + 1, i + 51)
        altsIdx <- altsIdx[which((annot$ARB2017[altsIdx, "AUID"] == auid) &
            (annot$ARB2017[altsIdx, "evalue"] <= eval))]
        if (length(altsIdx) > 0) {
            altsIdx <- altsIdx[grep("^uncultured", annot$ARB2017[altsIdx, "subject"], invert = T)]
        }
        if (length(altsIdx) > 0) {
            i <- altsIdx[1]
        } # Found a replacement(s)
                i
    }
)

cat(
    "\nThere were", length(idx), "AUIDs that had top hits in ARB with vague annotation (e.g.\n",
    "uncultured_microorganism).  I tried to find alternative hits with equally good E-values and\n",
    "more useful annotation.  I was able to do this for", sum(bestHitsIdx[idx] != alts), "AUIDs.\n"
)
bestHitsIdx[idx] <- alts
stopifnot(length(bestHitsIdx) == length(auids))

bestHits$ARB2017 <- annot$ARB2017[bestHitsIdx, c("AUID", "subject", "pid", "alen", "evalue")]


##
## Genome879 don't need alternatives -- we always know the gennome we hit:)
##
auids <- unique(annot$Genome879$AUID)
bestHitsIdx <- match(auids, annot$Genome879$AUID)
bestHits$Genome879 <- annot$Genome879[bestHitsIdx, c("AUID", "subject", "pid", "alen", "evalue")]

## Merge the blast annotations
df <- merge(x = bestHits$ARB2017, y = bestHits$Genome879, by = "AUID", all = T)
colnames(df) <- c(
    "AUID",
    "ARB2017.id", "ARB2017.pctId", "ARB2017.alen", "ARB2017.evalue",
    "Genome879.id", "Genome879.pctId", "Genome879.alen", "Genome879.evalue"
)
stopifnot(setequal(df$AUID, c(bestHits$ARB2017$AUID, bestHits$Genome879$AUID)))

## For top Genome879 hits, replace the id (which is an ARB name apparently)
## with the genus_species.  Also tack on taxonomic info.
cat("\nLooking up more informative names and taxa for the hits to Genome879.\n")
ghits <- as.character(df$Genome879.id)
idx <- match(ghits, g879taxa$arb_name)
ghits[!is.na(idx)] <- g879taxa[idx[!is.na(idx)], "genus_species"]
stopifnot(!g879taxa$genus_species %in% c(NA, "")) # So we did not make ghits worse.
df$Genome879.id <- ghits
## Now the taxa, in a single column. Separate levels with ';'.  Use 'unknown' for
## levels that are NA (only kingdom).
x <- apply(
    g879taxa[idx[!is.na(idx)], c("kingdom", "phylum", "class", "order", "family", "genus")], 1,
    function(x) {
        x[is.na(x)] <- "unknown"
        paste(x, collapse = ";")
    }
)
df[!is.na(idx), "Genome879.tax"] <- x


## Genome879 has NifH cluster or subcluster.  Here we could temporarily tack on
## the cluster and then compare it to what CART assigned.  But I don't know if
## the subcluster for the best hit in G879 should be the same as the CART
## assigned subcluster.  Often I see mismatches.
## For *cluster* mismatches are rare e.g. 63 mismatches between CART vs the best
## G879 hit, but 3707 matches.  The mismatches are often ~catch-all cluster "4"
## from CART but cluster 3 for the G879 hit.
## df[!is.na(idx),'TEMPCLUST'] <- sub('','',g879taxa[idx[!is.na(idx)],'nifHCluster'])


#############
##
## Marine cyano/NCD nifH DB
##
auids <- unique(annot$MarineDiazo$AUID)
bestHitsIdx <- match(auids, annot$MarineDiazo$AUID)
bestHits$MarineDiazo <- annot$MarineDiazo[bestHitsIdx, c("AUID", "subject", "pid", "alen", "evalue")]

## Tack on 'description' from the annotation companion to the FASTA
mdat <- bestHits$MarineDiazo
idx <- match(mdat$subject, cyanoNCDtaxa$seqID)
stopifnot(!is.na(idx)) # If trips, then the FASTA and annotation are out of sync.
mdat$description <- cyanoNCDtaxa[idx, "description"]

## Merge with the ARB and G879 blast annotations
colnames(mdat) <- c("AUID", paste0("MarineDiazo.", colnames(mdat)[-1]))
df <- merge(x = df, y = mdat, by = "AUID", all = T)
stopifnot(setequal(df$AUID, c(
    bestHits$ARB2017$AUID, bestHits$Genome879$AUID,
    bestHits$MarineDiazo$AUID
)))


###############
##
## UCYN-A oligos
##
auids <- unique(annot$UCYNAoligos$AUID)
bestHitsIdx <- match(auids, annot$UCYNAoligos$AUID)
bestHits$UCYNAoligos <- annot$UCYNAoligos[bestHitsIdx, c("AUID", "subject", "pid", "alen", "evalue")]

## Tack on 'description' from the annotation companion to the FASTA
mdat <- bestHits$UCYNAoligos
idx <- match(mdat$subject, ucynaOLIGOtaxa$seqID)
stopifnot(!is.na(idx)) # If trips, then the FASTA and annotation are out of sync.
mdat$description <- ucynaOLIGOtaxa[idx, "description"]

## Merge with the ARB and G879 blast annotations
colnames(mdat) <- c("AUID", paste0("UCYNAoligos.", colnames(mdat)[-1]))
df <- merge(x = df, y = mdat, by = "AUID", all = T)
stopifnot(setequal(df$AUID, c(
    bestHits$ARB2017$AUID, bestHits$Genome879$AUID, bestHits$MarineDiazo$AUID,
    bestHits$UCYNAoligos$AUID
)))


###################
##
## And now add CART annotation.
##
df <- merge(x = df, y = annot$CART, by = "AUID", all = T)
stopifnot(setequal(df$AUID, c(
    bestHits$ARB2017$AUID, bestHits$Genome879$AUID,
    bestHits$MarineDiazo$AUID, bestHits$UCYNAoligos$AUID,
    annot$CART$AUID
)))
x <- table(df$AUID)
if (any(x > 1)) {
    ## We might get here if AUID had multiple ORFs --> multiple nifH clusters (likely one valid and the
    ## other "ERROR" as returned by CART).
    cat("There should be one row of annotation per AUID but that is not the case. The following\n",
        "AUIDs have multiple annotations:\n"
    print(x[x > 1])
    stop("Aborting. Please check your CART results/log for AUIDs with multiple nifH clusters",
         "(due to multiple open reading frame predictions) and contact the authors for help.")
}
rm(x)

#############
##
## Conserved cysteines check (for C's that coordinate the 4Fe-4S cluster)
## (1) For AUIDs that have annotation from above sources, include the Switch I/II pattern checks.
## (2) For AUIDs without annotation, add them to 'df' only if the full CCAMP pattern is there.
##

auids1 <- df$AUID                                     # AUIDs with annotation (1)
auids2 <- subset(annot$CCAMP, hasCCAMP.len > 0)$AUID  # For (2), AUIDs with full pattern...
auids2 <- setdiff(auids, auids1)                      # ...but no annotation
df <- merge(x = df, y = subset(annot$CCAMP, AUID %in% c(auids1, auids2)), by = "AUID", all = T)
stopifnot(nrow(df) == length(auids1) + length(auids2))

cat("\nThere are", length(auids1), "AUIDs with annotation from the nifH reference DBs or CART.\n")
if (length(auids2) > 0) {
    n <- length(auids2)
    cat("There are", n, "with no DB or CART annotation but which appear to have cysteines to\n",
        "coordinate the 4Fe-4S cluster and AMP after the N-terminal cysteine.  These AUIDs will\n",
        "be included in the annnotation table.\n")
} else {
    ## Since annot$CCAMP has a result for each AUID, we can calc AUIDs that lacked annotation.
    n <- nrow(annot$CCAMP) - length(auids1)
    cat("All", n, "AUIDs that lack annotation from the DBs or CART also appear to lack paired\n",
        "cysteines (and \"AMP\") needed to coordinate the 4Fe-4S cluster.\n")
}


ofile <- "auids.annot.raw.tsv"
cat("\nWriting out the big annotation table", ofile, "with", nrow(df), "AUIDS.\n")
write.table(df, ofile, quote = F, sep = "\t", row.names = F)
