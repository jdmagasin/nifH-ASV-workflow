#!/usr/bin/env Rscript

##
## Make the nifH ASV database (nifH_ASV_database.tgz) and an associated workspace.RData that
## can be used to jumpstart subsequent analysis.  The following kinds of information are in
## the database:
##    AUIDs:    Sequences, abundance counts, annotation
##    Samples:  Associated metadata and environmental data
## and are described in workspaceObjectDescriptions below.
##
## Usage*:
##   make_nifH_ASV_database.R  asvAbunds.tsv  asv.fas  asvAnnot.tsv  sampMeta.tsv  cmapData.csv
##
## *Intended to be called by the Makefile. Not a tool.
##

workspaceObjectDescriptions <- "
  These objects comprise the nifH ASV database:
     asvSeqs      ASV sequences, 1:1 with abundTab. Each ASV has an ID of the form 'AUID.<number>'.
                  The number only ensures uniqueness.  It does not indicate the ASV's abundance.
     abundTab     ASV abundance counts for all retained samples.  Every count is an integer so the
                  able will work with the R vegan package.
     relabundTab  ASV relative abundances for all retained samples.
     cmapTab      Environmental data for samples in abundTab.
     metaTab      Metadata for samples in abundTab.
     annotTab     Annotation for ASVS in abundTab.

  The workspace.RData includes the above as well as the following:
     asvCyanos    Cyanobacteria ASVs based on best hit in Genomes879.
     asvNCDs      Non-cyanobacgteria ASVs based on best hit in Genomes879.
     GetTaxa()    Find ASVs with a specified taxonomic level (kingdom to genus) based on the ASV's
                  best hit in Genomes879.  Documentation is within the function comments.
"		  


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
    ## fixme: If going to do this, perhaps should drop these samples from metaTab.
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
asvSeqs <- tibble(AUID     = sub('^>(AUID.[^ ]+) .*$','\\1', fas[idx.defs]),
                  sequence = fas[idx.seqs])
cat("Loaded",nrow(asvSeqs),"ASVs.  ")
stopifnot(setequal(asvSeqs$AUID, abundTab$AUID))
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

##
## Make sure tables are coherent: AUIDs and samples work across tables.
## Compare to abundTab.
##
cat("Checking consistency of data tables (e.g. same ASVs and samples).\n")
stopifnot(setdiff(colnames(abundTab), metaTab$SAMPLEID) == 'AUID')  # Metadata for all samps
stopifnot(setdiff(colnames(abundTab), cmapTab$SAMPLEID) == 'AUID')  # CMAP for all samps
stopifnot(colnames(abundTab) == colnames(relabundTab))              # samples are 1:1
stopifnot(abundTab$AUID == relabundTab$AUID)                        # ASVs are 1:1 in abund tables
stopifnot(setequal(asvSeqs$AUID,abundTab$AUID))                     # ASVs are 1:1 with fasta
stopifnot(asvCyanos %in% abundTab$AUID)
stopifnot(asvNCDs %in% abundTab$AUID)

save(asvSeqs, abundTab, relabundTab, annotTab,
     metaTab, cmapTab, 
     asvCyanos, asvNCDs,
     GetTaxa,
     workspaceObjectDescriptions,
     file="workspace.RData")
cat("\nSaved workspace.RData which contains the following:")
cat(workspaceObjectDescriptions,"\n")
cat("The workspace objects can be used with R tidyverse, or without. No R packages are required.\n\n")

## Create the nifH ASV database
wdir <- 'nifH_ASV_database'
dir.create(wdir)
x <- file.copy(args['fasta'], file.path(wdir,'asvSeqs.fasta'))  # FASTA with original deflines, not asvSeqs
write_csv(abundTab,           file.path(wdir,'abundTab.csv'))
write_csv(relabundTab,        file.path(wdir,'relabundTab.csv'))
write_csv(annotTab,           file.path(wdir,'annotTab.csv'))
write_csv(metaTab,            file.path(wdir,'metaTab.csv'))
write_csv(cmapTab,            file.path(wdir,'cmapTab.csv'))
writeLines(workspaceObjectDescriptions, file.path(wdir,'manifest.txt'))
cat("Wrote files comprising the nifH ASV database.\n")

quit(save="no")
