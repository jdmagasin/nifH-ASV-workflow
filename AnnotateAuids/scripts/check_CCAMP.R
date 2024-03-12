#!/bin/env Rscript

## Copyright (C) 2024 Jonathan D. Magasin
##
## Check a FASTA with AUID amino sequences for a pair of cysteines, found in the highly conserved
## Switch I and II regions of NifH to coordinate the 4Fe-4S cluster.  Sequences without the
## cysteines are unlikely to be NifH.
##
## Usage:
##     check_CCAMP.R  <ORF FASTA>
## where the FASTA has amino sequences on a single line.
##
## In each sequence the following pattern is searched for: C [30-45 residues] C [1-3 residues] AMP
## The AMP usually occurs in Switch II often within FAMPIRE.  For example see Figure in Schlessman
## et al. 1998 (DOI: 10.1006/jmbi.1998.1898)
##
## Output: ccamp.csv, a table with three columns:
##     AUID:      AUIDs extracted from FASTA definition lines.  AUIDs have the form AUID.<digits>
##                so.  Example:  A FragGeneScan defintion line that starts ">AUID.12345_2_325_+"
##                would have "AUID.12345" in the AUID column.
##     orfLen:    Number of amino acids in the ORF for this AUID.
##     patLen:    Number of amino acids in the matched pattern.  0 if the pattern is absent.
##
## If an AUID appears in the FASTA multiple times (because it has multiple ORFs), then only the
## longest match will be recorded.
##

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

## Search each sequence for the C-C-AMP pattern.
ccampPat <- 'C.{30,45}C.{1,3}AMP'
m <- regexpr(ccampPat, orfs)
df <- data.frame(AUID = auids,
                 orfLen = sapply(orfs, nchar),
                 patLen = attr(m, 'match.length'))
## If pattern absent, the match.length is -1.
df$patLen[df$patLen == -1] <- 0
cat("The pattern", ccampPat, "was detected in", sum(df$patLen > 0), "of", nrow(df), "sequences.\n")
stopifnot( m[df$patLen == 0] == -1 )

## If an AUID has multiple ORFs, report only the longest match.  Here use match() to find the first
## occurrence of each AUID, which is the one wanted thankns to order().
df <- df[order(df$patLen, decreasing=T),]
idx <- match(unique(df$AUID), df$AUID)
if (length(idx) < nrow(df)) {
    cat(nrow(df) - length(idx), "AUIDs have multiple matches to the pattern.  ",
        "Will keep only the longest matches.\n")
    df <- df[idx,]
}

write.csv(df, 'ccamp.csv', row.names=F, quote=F)
cat("Saved ccamp.csv.\n")
