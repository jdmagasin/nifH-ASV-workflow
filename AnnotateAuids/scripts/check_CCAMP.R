#!/bin/env Rscript

## Copyright (C) 2024 Jonathan D. Magasin
##
## Check a FASTA with AUID amino sequences for a pair of highly conserved cysteines in NifH that
## coordinate the 4Fe-4S cluster.  Sequences without the cysteines are unlikely to be NifH.  Also
## check for AMP which usually follows the C in Switch II.
##
## Usage:
##     check_CCAMP.R  <ORF FASTA>
## where the FASTA has amino sequences on a single line.
##
## In each sequence two patterns are searched for:
##   C [30-42 residues] C
##   C [30-42 residues] C [1-3 residues] AMP
## Figure 1 in Schlessman et al. 1998 (DOI: 10.1006/jmbi.1998.1898) shows the two cysteines that
## coordinate the 4Fe-4S cluster, separated by 34 residues in A.vinelandii.  The only other cysteine
## we expect in our AUIDs (amplicons between nif2 and nif1 primers) is 12 residues nearer the amino
## terminus.  The patterns exclude this cysteine by allowing at most 34 + 8 = 42 residues between
## the two coordinating cysteines.
## 
## The AMP in Switch II is usually present (usually in FAMPIRE), but we have seen some ASVs with
## substitutions.
##

##
## Output: ccamp.csv, a table with three columns:
##    AUID:          AUIDs extracted from FASTA definition lines.  AUIDs have the form AUID.<digits>
##                   so.  Example:  A FragGeneScan defintion line that starts ">AUID.12345_2_325_+"
##                   would have "AUID.12345" in the AUID column.
##    orf.len:       Number of amino acids in the ORF for this AUID.
##    hasCC.len:     Number of amino acids matched in first pattern, or 0 if pattern absent.
##    hasCCAMP.len:  Number of amino acids matched in second pattern, or 0 if pattern absent.
##
## If an AUID appears in the FASTA multiple times (because it has multiple ORFs), then only the
## longest match will be recorded.
##

## Regular expressions used to search AUID ORFs.
ccPat    <- 'C.{30,42}C'
ccampPat <- 'C.{30,42}C.{1,3}AMP'

## Parse commmand line arguments
orfFasta <- commandArgs(T)[1]
stopifnot(file.exists(orfFasta))

## Read amino fasta and verify that sequences are on single lines.
lines <- readLines(orfFasta)
lines <- lines[lines != '']  # Drop blank lines
defLines <- grep('^>', lines)
cat("Loaded", orfFasta, "which has", length(defLines),"sequences.\n")
stopifnot(length(defLines) == length(lines)/2)

## Extract AUIDs from the definition lines.
auids <- sub('^>(AUID\\.[0-9]+).*$', '\\1', lines[defLines])
orfs <- lines[-defLines]
rm(lines, defLines)

## Search each sequence for the two patterns.  For simplicity search for them independently.
## Technically it's possible that the CC and CCAMP patterns could be found in different parts of the
## sequence, in which case hasCCAMP.len and hasCC.len might not differ by 4-6 as expected.  Can
## handle this unlikely issue but the code will become a little messy/complicated.
mcc <- regexpr(ccPat, orfs)
mccamp <- regexpr(ccampPat, orfs)
df <- data.frame(AUID         = auids,
                 orf.len      = sapply(orfs, nchar),
                 hasCC.len    = attr(mcc,    'match.length'),
                 hasCCAMP.len = attr(mccamp, 'match.length'))
## If pattern absent, the match.length is -1.
df$hasCC.len[df$hasCC.len == -1] <- 0
df$hasCCAMP.len[df$hasCCAMP.len == -1] <- 0
cat("Pattern", ccPat,    "was detected in", sum(df$hasCC.len > 0),    "of", nrow(df), "sequences.\n")
cat("Pattern", ccampPat, "was detected in", sum(df$hasCCAMP.len > 0), "of", nrow(df), "sequences.\n")
stopifnot(mcc[df$hasCC.len == 0] == -1)
stopifnot(mccamp[df$hasCCAMP.len == 0] == -1)

## If an AUID has multiple ORFs, report only the longest match.  Here use match() to find the first
## occurrence of each AUID, which is the one wanted thanks to order().
matchLens <- apply(df[,c('hasCC.len','hasCCAMP.len')], 1, max)  # For each AUID use longer match
df <- df[order(matchLens, decreasing=T),]
idx <- match(unique(df$AUID), df$AUID)
if (length(idx) < nrow(df)) {
    cat(nrow(df) - length(idx), "AUIDs contain multiple occurrences of the patterns.  ",
        "Will keep only the longest pattern matches for each AUID.\n")
    df <- df[idx,]
}

write.csv(df, 'ccamp.csv', row.names=F, quote=F)
cat("Saved ccamp.csv.\n")
