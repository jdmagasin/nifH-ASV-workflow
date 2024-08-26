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
##   make_nifH_ASV_database.R  <asvAbunds.tsv>  <asv.fas>  <asvAnnot.tsv>  <sampMeta.tsv>  <cmapData.csv>
##
## *Intended to be called by the Makefile. Not a tool.
##

workspaceObjectDescriptions <- "
  These objects comprise the nifH ASV database:
     asvSeqs      ASV sequences, 1:1 with abundTab. Each ASV has an ID of the form 'AUID.<number>'.
                  The number only ensures uniqueness.  It does not indicate the ASV's abundance.
     abundTab     ASV abundance counts for all retained samples. Samples with 0 reads are excluded,
                  for example, those dropped by FilterAuids, or by WorkspaceStartup when ASVs
                  with no annotation get removed.  Every count is an integer so the table can be
                  used with the R vegan package.
     relabundTab  ASV relative abundances for all retained samples.
     annotTab     Annotation for ASVs in abundTab.
     metaTab      Metadata for samples in abundTab.
     cmapTab      Environmental data for samples in abundTab.

  The workspace.RData includes the above as well as 
     GetTaxaGenome879()    Find ASVs with a specified taxonomic level (kingdom to genus) based on
                           the ASV's best hit in the Genome 879 database.  For low levels or
                           specific ASVs, also consider the quality of the Genome 879 hit which is
                           provided in the annotTab fields Genome879.{pctId,alen,evalue}.  In the
                           annotTab see also the primary_id column since that can draw from any of
                           the resources used during AnnotateAuids.
"


## ------------------------------------------------------------------------------
##
## Start up
##

options(width = 110) # Assume 110 cols for wide-ish table printing.

## 16 Nov 2023: For now, do not drop samples from the abundance table that are missing
## CMAP or meta-data. They could still be dropped if they have 0 total reads because
## the CMAP and meta-data tables must remain in sync with abundance tables.
KEEP_SAMPS_MISSING_META_OR_CMAP_DATA <- TRUE

args <- commandArgs(T)
if (FALSE) {
    ## DEBUGGING
    args <- c(
        "../FilterAuids/auid.abundances.filtered.tsv.gz",
        "../FilterAuids/auid.filtered.fasta",
        "../GatherMetadata/metadata.csv",
        "../AnnotateAuids/auids.annot.tsv",
        "../CMAP/data/CMAP_metadata.csv.gz"
    )
}
stopifnot(length(args) == 5)
names(args) <- c("abundTabTsv", "fasta", "annotTsv", "metaCsv", "envarCsv")
stopifnot(file.exists(args))

cat("Loading libraries...")
suppressMessages(library(tidyverse))
cat("done.\n")


## ------------------------------------------------------------------------------
##
## Load tables
##

## Tidyverse is adverse to rownames. read_tsv() and friends do not properly read table with
## rownames: They put the row names into column 1, so for the abundance table we get AUID labels as
## the values for the first sample, and the final sample is garbage (all "0\t0").  I see no params
## in read_tsv() to work around this nor solutions online.  Better to use read.table() and convert
## to a tibble (assuming we prefer tibbles).

cat("Loading the abundance table...")
abundTab <- read.table(args["abundTabTsv"], stringsAsFactors = T, header = T, row.names = 1, check.names = F)
abundTab$AUID <- rownames(abundTab)
## Rename columns to use just the sample ID and make AUID column 1
abundTab <- as_tibble(abundTab) %>%
    rename_with(~ str_replace(.x, "^.+___", "")) %>%
    select(AUID, everything())
cat("done. ", nrow(abundTab), "ASVs X", ncol(abundTab) - 1, "working samples.\n")
cat("Note that", length(which(abundTab %>% select(-AUID) %>% colSums() == 0)),
    "samples have 0 total reads. They will be dropped after ASV annotations are checked.\n")


cat("Loading sample metadata...")
metaTab <- as_tibble(read.table(args["metaCsv"], stringsAsFactors = T, header = T, sep = ","))
cat("done. ", nrow(metaTab), "samples X", ncol(metaTab) - 1, "metadata variables.\n")
## Drop "_transcriptomic" which GatherAsvs appended to transcriptomic samples.
x <- sub("_transcriptomic$", "", metaTab$SAMPLEID)
## Currently we cannot handle DNA and cDNA sequencing runs with same sample ID. Would have to figure
## out which column in the abundance table is for DNA vs. cDNA (and allowing just one for each?)
stopifnot(table(x) == 1)
metaTab$SAMPLEID <- x


cat("Loading the CMAP environmental variables...")
## Use read.table() because it is ~instantaneous (vs. seconds) and avoids renaming of each column
## with ...<colNumber>.  Drop "_transcriptomic" as did for metaTab.  Drop columns that start with
## "FIXME." (e.g. FIXME.depth) which were added by prepareMetadataForCmap.R to flag samples that
## required (manual) fixing of illegal values, e.g. "surface" for the "depth".  Such cases should
## have been fixed.  If not then the CMAP stage would have dropped such samples with a warning.
cmapTab <- as_tibble(read.table(args["envarCsv"], stringsAsFactors = T, header = T, sep=",")) %>%
  mutate(SAMPLEID = str_remove(SAMPLEID, "_transcriptomic$")) %>%
  select(!starts_with("FIXME."))
cat(
    "done. ", nrow(cmapTab), "samples X", ncol(cmapTab) - 1,
    "CMAP environmental variables.\n\n"
)


cat("Loading annotation...")
## Load annotation for ASVs in the abundance table.  Drop ASVs with "unknown"
## primary_id unless they have cysteines that might coordinate the 4Fe-4S and
## the AMP (often in FAMPIRE) in the Switch II region.  Split the Genome879 taxa
## string.  Sort by AUID number (consistent with filter_by_annotTab_auid()).
annotTab <- read_tsv(args["annotTsv"], col_types = cols()) %>%
    filter((!primary_id %in% c("unknown", "unknownERROR")) | (hasCCAMP.len > 0)) %>%
    separate(
      col = Genome879.tax,
      sep = ";", paste0("Genome879.", c("k", "p", "c", "o", "f", "g"))
  ) %>%
  arrange(as.numeric(str_extract(AUID,"\\d+$")))
cat("done.\n")

## AnnotateAuids checks for AUIDs with multiple rows, but good to check again.
annotCount <- table(annotTab$AUID)
if (any(annotCount) > 1) {
    cat("The following AUIDs have multiple rows in the annotation table:\n")
    print(annotCount[annotCount > 1])
    stop("Aborting. Please check the log from AnnotateAuids and contact the authors for help.")
}
rm(annotCount)


## Very simple FASTA reader which assumes valid FASTA created by FilterAuids
## (which used ShortRead).  Assumes one line for nucleic acids.
## Create vector 'auids' which has values that are the AUID sequences and names
## that are the AUID identifiers (AUID.<num>).
cat("Loading AUID sequences...")
fas <- readLines(args["fasta"])
fas <- fas[fas != ""] # Drop empty lines
## Checks that lines alternate AUID, then sequence.
idx.defs <- grep("^>AUID", fas)
idx.seqs <- grep("^>AUID", fas, invert = T)
if ((length(idx.defs) != length(idx.seqs)) || !all(idx.defs + 1 == idx.seqs)) {
    stop(
        "The FASTA", args["fasta"], "seems to not have alternating definition",
        "and sequence lines."
    )
}
asvSeqs <- tibble(
    AUID = sub("^>(AUID.[^ ]+) .*$", "\\1", fas[idx.defs]),
    sequence = fas[idx.seqs]
)
cat("Loaded", nrow(asvSeqs), "AUIDs.\n")
rm(idx.defs, idx.seqs)


## ------------------------------------------------------------------------------
##
## Define some helpful functions
##

GetTaxaGenome879 <- function(lev, taxa, invert = FALSE) {
    ##
    ## Return a character vector of ASVs that have Genome 879 taxonomic level
    ## 'lev' that is equal to the specified 'taxa'.  Pass 'lev' an element in
    ## {k,p,c,o,f,g}.  For example to get all ASVs that are cyanobacterial and
    ## non-cyanobacterial diazotrophs based on their best Genome 879 hits:
    ##    asvCyanos.G879 <- GetTaxaGenome879('p','Cyanobacteria')
    ##    asvNCDs.G879   <- GetTaxaGenome879('p','Cyanobacteria', invert=T)
    ##
    stopifnot(lev %in% c("k", "p", "c", "o", "f", "g"))
    g879Taxon <- paste0("Genome879.", lev)
    if (invert) {
        ss <- annotTab %>% filter(!.data[[g879Taxon]] %in% taxa)
    } else {
        ss <- annotTab %>% filter(.data[[g879Taxon]] %in% taxa)
    }
    ## as_vector names to the vector entries (to "AUID1","AUID2",etc.) which
    ## is silly/redundant for our purpose. Strip the names with as.character.
    ss %>%
        select(AUID) %>%
        as_vector() %>%
        as.character()
}


## ------------------------------------------------------------------------------
##
## Filter AUIDs from the nifH database objects that did not receive an annotation
## during AnnotateAuids.  Annotation could be from the searched DBs, CART, or
## the presence of the the C-C-AMP pattern (conserved cysteines to coordinate
## the 4Fe-4S cluster).
cat("\nFiltering out AUIDs that had no annotation from workflow stage AnnotateAuids:\n")

## Create a function that filters based on annotation AUIDs
## Returns a filtered df
filter_by_annotTab_auid <- function(df) {
  ## Make a key of the AUIDs from the annotation table to filter input df.
  annotTab_auid_key <- annotTab %>%
    select(AUID) %>%
    pull()
  # filter input df with auid key and sort by auid number
  filt_df <- df %>%
    filter(AUID %in% annotTab_auid_key)  %>%
    arrange(as.numeric(str_extract(AUID,"\\d+$")))
  # Create print statement for number of rows removed
  rows_removed = nrow(df) - nrow(filt_df)
  cat(paste0("  Number of rows removed from: '", deparse(substitute(df)), "':\t", rows_removed,"\n"))

  return(filt_df)
}

# Filter abundance table and sequence file, removing unannotated ASVs
abundTab <- filter_by_annotTab_auid(abundTab)
asvSeqs <- filter_by_annotTab_auid(asvSeqs)
cat("After filtering out unannotated ASVs, there are a total of",
    length(which(abundTab %>% select(-AUID) %>% colSums() == 0)),
    "samples that have 0 total reads. They will be dropped shortly.\n")

# Verify that all objects and their associated files have the same AUIDs in the same order.
stopifnot(identical(asvSeqs$AUID, abundTab$AUID) &&
          identical(asvSeqs$AUID, annotTab$AUID))

# Convert filtered ASV sequences to FASTA format
stopifnot(setequal(asvSeqs$AUID, abundTab$AUID))
cat("ASVs in FASTA and abundance table are 1:1.\n")
asvSeqs_fasta <- paste0(">", asvSeqs$AUID, "\n", asvSeqs$sequence, collapse = "\n")


## ------------------------------------------------------------------------------
##
## Ensure that abundTab is integer-only.
##
## table(sapply(abundTab, class))  # Shows that there are integers and doubles.
cat("\nForcing abundTab to have only integers to be vegan-friendly.\n")
tot.initial <- abundTab %>%
    column_to_rownames("AUID") %>%
    as.matrix() %>%
    sum()
int_if_need <- function(v) {
    if (!is.integer(v)) {
        v <- as.integer(round(v))
    }
    v
}
abundTab <- abundTab %>% mutate(across(-AUID, int_if_need))
tot.final <- abundTab %>%
    column_to_rownames("AUID") %>%
    as.matrix() %>%
    sum()
x <- abs(tot.final - tot.initial)
if (x > 0) {
    cat(
        "This changed the total number of reads from", tot.initial, "to", tot.final, "\n",
        "a change of", x, "reads.\n"
    )
} else {
    cat("After this step there are still", tot.initial, "reads.\n")
}


##
## Removal of ASVs lacking annotation could (and does) renders some samples empty.  Also, a few
## samples are empty when received by WorkspaceStartup since FilterAuids does not do a final pass
## after its ASV filtering.  So here we remove any completely empty samples.
##
x <- names(which(abundTab %>% select(-AUID) %>% colSums() == 0))
if (length(x) > 0) {
    cat("Dropping samples that have 0 reads.  ")
    abundTab    <- abundTab %>% select(!contains(x))
    cat(length(x), "samples dropped:\n")
    cat(strwrap(x),"\n\n")
}

cat("Shrinking CMAP and sample metadata tables to have just the samples in the abundance table.\n")
sampIds <- setdiff(colnames(abundTab), "AUID")
cmapTab <- cmapTab %>% filter(SAMPLEID %in% sampIds)
missingIdx <- idx <- which(!sampIds %in% cmapTab$SAMPLEID)
cat(
    "  ", length(idx), "samples in the abundance table have no environmental data",
    "(sampsWithoutEnvdata.txt)\n"
)
writeLines(sampIds[idx], "sampsWithoutEnvdata.txt")

metaTab <- metaTab %>% filter(SAMPLEID %in% sampIds)
idx <- which(!sampIds %in% metaTab$SAMPLEID)
cat(
    "  ", length(idx), "samples in the abundance table have no sample metadata",
    "(sampsWithoutMetadata.txt)\n\n"
)
writeLines(sampIds[idx], "sampsWithoutMetadata.txt")
missingIdx <- union(idx, missingIdx)


if (!KEEP_SAMPS_MISSING_META_OR_CMAP_DATA && length(missingIdx) > 0) {
    cat(
        "Dropping", length(missingIdx), "samples from the ASV abundance table",
        "that lack environmental and/or metadata.\n\n"
    )
    abundTab <- abundTab %>% select(!contains(sampIds[missingIdx]))
    ## fixme: If going to do this, perhaps should drop these samples from metaTab.
}


## ------------------------------------------------------------------------------
##
## Make a relative abundance table
##

relabundTab <- abundTab %>% mutate_at(vars(-AUID), ~ . / sum(.))
## Verify all columns sum to ~1, and all cols are type double.
stopifnot(abs((relabundTab %>% select(-AUID) %>% colSums()) - 1) < 1e-9)
stopifnot(sapply(relabundTab[, -1], class) == "numeric")


##
## Make sure tables are coherent: AUIDs and samples work across tables.
## Compare to abundTab.
##
cat("Checking consistency of data tables (e.g. same AUIDs and samples).\n")
if (!KEEP_SAMPS_MISSING_META_OR_CMAP_DATA) {
    stopifnot(setdiff(colnames(abundTab), metaTab$SAMPLEID) == "AUID") # Metadata for all samps
    stopifnot(setdiff(colnames(abundTab), cmapTab$SAMPLEID) == "AUID") # CMAP for all samps
}
stopifnot(colnames(abundTab) == colnames(relabundTab)) # Samps 1:1 in abund tables
stopifnot(identical(asvSeqs$AUID, abundTab$AUID))      # ASVs 1:1 in abund tables and FASTA
stopifnot(identical(asvSeqs$AUID, relabundTab$AUID))
stopifnot(identical(asvSeqs$AUID, annotTab$AUID))

save(asvSeqs, abundTab, relabundTab, annotTab,
    metaTab, cmapTab,
    GetTaxaGenome879,
    workspaceObjectDescriptions,
    file = "workspace.RData"
)
cat("\nSaved workspace.RData which contains the following:")
cat(workspaceObjectDescriptions, "\n")
cat("The workspace objects can be used with R tidyverse, or without. No R packages are required.\n\n")

## Create the nifH ASV database
wdir <- "nifH_ASV_database"
dir.create(wdir)
# Filtered versions of each nifH ASV database object
writeLines(asvSeqs_fasta,     file.path(wdir,"asvSeqs.fasta"))  
write_csv(abundTab,           file.path(wdir,'abundTab.csv'))
write_csv(relabundTab,        file.path(wdir,'relabundTab.csv'))
write_csv(annotTab,           file.path(wdir,'annotTab.csv'))
write_csv(metaTab,            file.path(wdir,'metaTab.csv'))
write_csv(cmapTab,            file.path(wdir,'cmapTab.csv'))
writeLines(workspaceObjectDescriptions, file.path(wdir,'manifest.txt'))
cat("Wrote files comprising the nifH ASV database.\n")

quit(save = "no")
