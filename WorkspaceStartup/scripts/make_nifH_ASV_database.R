#!/usr/bin/env Rscript

##
## Make the nifH ASV database (nifH_ASV_database.tgz) and an associated workspace.RData that
## can be used to jumpstart subsequent analysis.  The following kinds of information are in
## the database:
##  AUIDs:    Sequences, abundance counts, annotation
##  Samples:  Associated metadata and environmental data
##
## Main objects saved:
##     auids        ASV sequences.  1:1 with abundTab
##     abundTab     ASV abundance counts for all retained samples.  Every count is an
##                  integer so the table will work with the R vegan package.  (Normalization
##                  used to make the input abundance tables non-integer.  The code for
##                  converting to integer is retained.)
##     relabundTab  ASV relative abundances for all retained samples
##     cmapTab      Environmental data for samples in abundTab
##     metaTab      Metadata for samples in abundTab
##     annotTab     Annotation for ASVS in abundTab
##
##     asvCyanos    IDs of all cyanobacteria (best hit in Genomes879)
##     asvNCDs      IDs of all non-cyanobacteria (best hit in Genomes879)
##     GetTaxa()    Function used to find ASVs from a specified taxonomic level.
##                  This was used to create asvCyanos and asvNCDs. Documentation
##                  within the function (as comments).
##
## Usage*:
##   make_nifH_ASV_database.R  asvAbunds.tsv  asv.fas  asvAnnot.tsv  sampMeta.tsv  cmapData.csv
##
## *Intended to be called by the Makefile. Not a tool.
##

##------------------------------------------------------------------------------
##
## Start up
##
options(width = 110)  # Assume 110 cols for wide-ish table printing.
args <- commandArgs(T)
if (FALSE) {
    ## DEBUGGING
    args <- c('../FilterAuids/auid.abundances.filtered.tsv.gz',
              '../FilterAuids/auid.filtered.fasta',
              '../GatherMetadata/metadata.tsv',
              '../AnnotateAuids/auids.annot.tsv',
              '../CMAP/CMAP_data.csv.gz')
}
stopifnot(length(args) == 5)
names(args) <- c('abundTabTsv', 'fasta', 'annotTsv', 'metaTsv', 'envarCsv')
stopifnot(file.exists(args))

cat("Loading libraries...")
suppressMessages(library(tidyverse))  # FIXME: add to Installation/environment_nifH_ASV_workflow.yml
cat("done.\n")


##------------------------------------------------------------------------------
##
## Load tables
##

## Tidyverse is adverse to rownames. read_tsv() and friends do not properly read
## table with rownames: They put the row names into column 1, so for the
## abundance table we get AUID labels as the values for the first sample, and
## the final sample is garbage (all "0\t0").  I see no params in read_tsv() to
## work around this nor solutions online.  Better to use read.table() and
## convert to a tibble (assuming we prefer tibbles).

cat("Loading the abundance table...")
abundTab <- read.table(args['abundTabTsv'], stringsAsFactors=T, header=T, row.names=1)
abundTab$AUID <- rownames(abundTab)
## Rename columns to use just the sample ID and make AUID column 1
abundTab <- as_tibble(abundTab) %>% 
    rename_with(~ str_replace(.x, '^.+___','')) %>%
    select(AUID, everything())
cat("done. ",nrow(abundTab),"ASVs X",ncol(abundTab)-1,"working samples.\n")


cat("Loading the CMAP environmental variables...")
cmapTab <- read_csv(args['envarCsv'], col_types=cols())  # use all cols and avoids msgs
cat("done. ",nrow(cmapTab),"samples X",ncol(cmapTab)-1,
    "CMAP environmental variables.\n")


cat("Loading sample metadata...")
metaTab <- as_tibble(read.table(args['metaTsv'], stringsAsFactors=T, header=T,sep="\t"))
cat("done. ",nrow(metaTab),"samples X",ncol(metaTab)-1,"metadata variables.\n")

x <- names(which(abundTab %>% select(-AUID) %>% colSums() == 0))
if (length(x) > 0) {
    cat("Dropping samples that have 0 reads.  ")
    abundTab <- abundTab %>% select(!contains(x))
    cat(length(x), "samples dropped:\n")
    strwrap(x)
}


cat("Shrinking CMAP and sample metadata tables to have just the samples in the abundance table.\n")
sampIds <- setdiff(colnames(abundTab),'AUID')
cmapTab <- cmapTab %>% filter(SAMPLEID %in% sampIds)
missingIdx <- idx <- which(! sampIds %in% cmapTab$SAMPLEID)
cat(length(idx),"samples in the abundance table have no environmental data",
    "(sampsWithoutEnvdata.txt)\n")
writeLines(sampIds[idx],'sampsWithoutEnvdata.txt')

metaTab <- metaTab %>% filter(SAMPLEID %in% sampIds)
idx <- which(! sampIds %in% metaTab$SAMPLEID)
cat(length(idx),"samples in the abundance table have no sample metadata",
    "(sampsWithoutMetadata.txt)\n")
writeLines(sampIds[idx],'sampsWithoutMetadata.txt')
missingIdx <- union(idx, missingIdx)


if (length(missingIdx) > 0) {
    cat("Dropping", length(missingIdx), "samples from the ASV abundance table",
        "that lack environmental and/or metadata.\n")
    abundTab <- abundTab %>% select(!contains(sampIds[missingIdx]))
}


cat("Loading annotation...")
## Load annotation for ASVs in the abundance table, and split the Genomes879
## taxa string.
annotTab <- read_tsv(args['annotTsv'], col_types=cols()) %>%
  filter(AUID %in% abundTab$AUID) %>%
  separate(col = Genomes879.tax,
           sep = ';', paste0("Genomes879.",c('k','p','c','o','f','g')))
cat("done.\n")


## Very simple FASTA reader which assumes valid FASTA created by FilterAuids
## (which used ShortRead).  Assumes one line for nucleic acids.
## Create vector 'auids' which has values that are the AUID sequences and names
## that are the AUID identifiers (AUID.<num>).
cat("Loading ASV sequences...")
fas <- readLines(args['fasta'])
fas <- fas[fas != '']  # Drop empty lines
## Checks that lines alternate AUID, then sequence.
idx.defs <- grep('^>AUID',fas)
idx.seqs <- grep('^>AUID',fas,invert=T)
if ((length(idx.defs) != length(idx.seqs)) || !all(idx.defs + 1 == idx.seqs)) {
    stop("The FASTA",args['fasta'],"seems to not have alternating definition",
          "and sequence lines.")
}
auids <- tibble(AUID     = sub('^>(AUID.[^ ]+) .*$','\\1', fas[idx.defs]),
                sequence = fas[idx.seqs])
cat("Loaded",nrow(auids),"ASVs.  ")
stopifnot(setequal(auids$AUID, abundTab$AUID))
cat("ASVs in FASTA and abundance table are 1:1.\n")
rm(idx.defs, idx.seqs)


##------------------------------------------------------------------------------
##
## Define some helpful functions
##

GetTaxa <- function(lev, taxa, invert=FALSE)
{
    ##
    ## Return a character vector of ASVs that have taxonomic level 'lev' that is
    ## equal to the specified 'taxa'.  Pass 'lev' an element in {k,p,c,o,f,g}.
    ## For example to get all cyanobacterial and non-cyanobacterial diazotrophs:
    ##    asvCyanos <- GetTaxa('p','Cyanobacteria')
    ##    asvNCDs   <- GetTaxa('p','Cyanobacteria', invert=T)
    ##
    ## Taxomomic information is based on the ASV's best hit to Genomes879, which
    ## is recorded in the 'annotTab'.
    ##
    stopifnot(lev %in% c('k','p','c','o','f','g'))
    g879Taxon <- paste0("Genomes879.", lev)
    if (invert) {
        ss <- annotTab %>% filter(! .data[[g879Taxon]] %in% taxa )
    } else {
        ss <- annotTab %>% filter( .data[[g879Taxon]] %in% taxa )
    }
    ## as_vector names to the vector entries (to "AUID1","AUID2",etc.) which
    ## is silly/redundant for our purpose. Strip the names with as.character.
    ss %>% select(AUID) %>% as_vector %>% as.character
}


##------------------------------------------------------------------------------
##
## Useful ASV lists
##
asvCyanos <- GetTaxa('p','Cyanobacteria')
asvNCDs   <- GetTaxa('p','Cyanobacteria', invert=T)


##------------------------------------------------------------------------------
##
## Make a relative abundance table
##

relabundTab <- abundTab %>% mutate_at(vars(-AUID),  ~ ./sum(.))
## Verify all columns sum to ~1, and all cols are type double.
stopifnot( abs((relabundTab %>% select(-AUID) %>% colSums()) -1) < 1e-9 )
stopifnot( sapply(relabundTab[,-1],class) == 'numeric' )

##------------------------------------------------------------------------------
##
## Ensure that abundTab is integer-only.
##
## table(sapply(abundTab, class))  # Shows that there are integers and doubles.
cat("\nForcing abundTab to have only integers to be vegan-friendly.\n")
tot.initial <- abundTab %>% column_to_rownames('AUID') %>% as.matrix() %>% sum()
int_if_need <- function(v) {
    if (!is.integer(v)) { v <- as.integer(round(v)) }
    v
}
abundTab <- abundTab %>% mutate(across(-AUID, int_if_need))
tot.final <- abundTab %>% column_to_rownames('AUID') %>% as.matrix() %>% sum()
x <- abs(tot.final - tot.initial)
if (x > 0) {
    cat("This changed the total number of reads from", tot.initial, "to", tot.final,"\n",
        "a change of", x, "reads.\n")
} else {
    cat("After this step there are still", tot.initial, "reads.\n")
}


save(auids, abundTab, relabundTab, annotTab,
     metaTab, cmapTab, 
     asvCyanos, asvNCDs,
     GetTaxa,
     file="workspace.RData")
cat("Saved workspace.RData\n")
quit(save="no")
