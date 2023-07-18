#!/usr/bin/env Rscript

##
## Scrape all wanted information available from all available metadata files
## listed in studyID_metadata.csv.  Some values are converted into standard
## formats e.g. Collection_Date.
##
## Important: In the first column of Every metadata file there must be unique
## sample identifiers.  The column can be named anything but will be renamed
## SAMPLEID by this script.  
##
##     A SAMPLEID must correspond to a specific point in time/space. This is
##     required to look up environmental data for the SAMPLEID.  Date, lat/lon,
##     and depth recorded for each SAMPLEID in the metadata files is probably
##     sufficient for looking up envdata.
##
## SAMPLEIDs must appear only once in a metadata file.  If a source metadata
## file has a duplicated SAMPLEID, you will have to fix it by:
##   - commenting out the duplicates by starting their lines with '#', or
##   - renaming the duplicates by adding _dup1, _dup2, ... to each.
## Some of our metadata files have duplicated SAMPLEIDs for metagenomic and
## metatranscriptomic samples (detected in field LibrarySource).  To handle
## this the script appends _transcriptomic to the sample IDs.  Unfortunately
## the modified sample IDs cannot be used as-is to index into other tables.
##

## SAMPLEIDs can appear in multiple metadata files. If metadata files share some
## columns (e.g. Collection_Date), this script will check simply for consistent
## values across metadata files.  A very few fields are put into standard
## formats before this script does the comparison (e.g. Collection_Date).
## However, for most fields non-identical values between metadata files will
## cause the script to report the conflict and abort.
##

cat("Reading studyID_metadata.csv...")
sidTab <- read.csv('studyID_metadata.csv', row.names=NULL, header=F, comment.char='#')
colnames(sidTab) <- c('studyid','filepath')
cat(nrow(sidTab),"files are listed.\n")

## Warn if any metadata files cannot be read.
x <- file.access(sidTab$filepath, 4)
if (any(x == -1)) {
    cat("The following tables cannot be read. Permission problem or the file is missing.\n")
    cat(paste(names(which(x==-1)),collapse="\n"),"\n")
    stop("Aborting.")
}
rm(x)


## Columns to scrape from metadata files. Names are the R-friendly final names
## and values are the regular expressions to use to identify the column in a
## metadata table. (The regexps are used in a case insensitive grep, after
## anchoring -- see "grep for the fields we want" loop below.)
MakeWantCols <- function()
{
    wantCols <- c(Run = 'Run',  SRA_Study = 'SRA[ _]*Study',
                  BioProject = 'BioProject', BioSample = 'BioSample',
                  ## Allow an alternative sample name but give it a distinct
                  ## field because this simple script is not supposed to
                  ## reconcile duplicate/ambiguous fields.
                  Sample_Name = 'Sample[ _]*Name', Sample_Name_2 = 'Sample[ _]*Name[ _]2',
                  Region = 'Region',
                  Collection_Date = 'Collection[ _]*Date',
                  Time = 'Time', Depth = 'Depth', Station = 'Station',
                  geo_loc_name = 'geo_loc_name',
                  Lat_Lon = 'Lat_Lon', Lat = 'Lat', Lon = 'Lon',
                  Latitude = 'Latitude', Longitude = 'Longitude',
                  LibraryName = 'LibraryName', LibrarySelection = 'LibrarySelection',
                  LibrarySource = 'LibrarySource',
                  Size_fraction = 'Size_fraction',
                  FilterSize = 'FilterSize',
                  Sample_Type = 'Sample[\\._]{0,1}Type')
      ## Anchor each pattern
      x <- names(wantCols)
      wantCols <- paste0('^',wantCols,'$')
      names(wantCols) <- x
      wantCols
}
wantCols <- MakeWantCols()
cat("Will use the first column of each metadata file as a SAMPLEID\n",
    "and will scan for these columns (using regular expressions):\n")
print(names(wantCols))
cat("\n")


## Fixup dates to be in format e.g. 2021-1-17.  Some examples of the formats
## fixed up:   1/7/21  1.7.21  1/07/21  1/17/2021
## Not bullet proof.
Fixup_Collection_Date <- function(dates)
{
    reg.ymd <- '20[0-9]{2}-[0-9]{1,2}-[0-9]{1,2}'               # Wanted: 2016-5-19
    reg.mdy <- '[0-9]{1,2}[\\/\\.][0-9]{1,2}[\\/\\.][0-9]{2,4}' # Fixup these

    ## First drop leading 0's for stuff that otherwise is in desired format.
    idx <- grep(reg.ymd, dates)
    if (length(idx) > 0) {
        x <- as.character(lapply(strsplit(dates[idx],'-',T),
                 function(v) {
                     y <- as.numeric(v[1])
                     m <- as.numeric(v[2])
                     d <- as.numeric(v[3])
                     v <- paste0(v,collapse='-')
                     if (y <= 2022 && m <= 12 && d <= 31) {
                         v <- paste0(y,'-',m,'-',d)
                     }
                     v
                 }))
        dates[idx] <- x
    }

    idx <- grep(reg.mdy, dates)
    if (length(idx) > 0) {
        x <- as.character(lapply(strsplit(dates[idx],'[\\/\\.]',F),
                 function(v) {
                     if (nchar(v[3]) == 2) {
                         y <- as.numeric(paste0('20',v[3]))
                     } else if (nchar(v[3]) == 4) {
                         y <- as.numeric(v[3])
                     } else stop("Weird year",v)
                     m <- as.numeric(v[1])
                     d <- as.numeric(v[2])
                     v <- paste0(v,collapse='-') # might change the separator
                     if (y <= 2022 && m <= 12 && d <= 31) {
                         v <- paste0(y,'-',m,'-',d)
                     }
                     v
                 }))
        dates[idx] <- x
    }
    dates
}


## HACKS!
##  1. Rename transcriptomic samples by adding a "_transcriptomic" to the end of their
##     SAMPLEID.  Some studies use the same sample ID for their metaT and metaG samples
##     so this is a simple way to ensure unique SAMPLEIDs.  The modified SAMPLEIDs will
##     not index into other tables unfortunately, but our focus is on metagenomes.
##
##  2. Change GENOMIC --> METAGENOMIC in LibrarySource
##
## If need to add more hacks, break into separate functions.
##
## Return fixed up mdat.
##
gMetaT <- c()
RecordMetaT <- function(x) {
    gMetaT <<- union(gMetaT, as.character(x))
}
HackHandleProblemSamples <- function(mdat)
{
    ## Modify SAMPLEIDs of transcriptomic, and record the original IDs.
    idx <- grep('LibrarySource', colnames(mdat), ignore.case=T)
    if (length(idx) > 0) {
        stopifnot(length(idx)==1)
        idx <- grep('transcriptomic', mdat[,idx], ignore.case=T)
        if (length(idx) > 0) {
            cat("\tWARNING: Found",length(idx),"transcriptomic samples. ",
                "Will append _transcriptomic to their sample IDs.\n")
            RecordMetaT(mdat[idx,1])
            mdat[idx,1] <- paste0(mdat[idx,1], '_transcriptomic')
        }
    }

    ## GENOMIC --> METAGENOMIC
    idx <- grep('LibrarySource',colnames(mdat),ignore.case=T)
    if (length(idx) > 0) {
        stopifnot(length(idx)==1)
        idx2 <- grep('^genomic$',mdat[,idx],ignore.case=T)
        if (length(idx2) > 0) {
            cat("\tWARNING: Translated 'genomic' -> 'metagenomic' for",
                length(idx2), "samples.\n")
            mdat[idx2,idx] <- 'METAGENOMIC'
            stopifnot(nrow(mdat) > 0)
        }
    }
    mdat
}


## Read csv or tsv metadata.
ReadMetaTable <- function(fp)
{
    naStr = c('','NA','N/A',"na","n/a")
    if (grepl('\\.csv$',fp)) {
        mdat <- read.csv(fp, check.names=F, row.names=NULL,
                         colClasses = "character", na.strings = naStr,
                         comment.char='#')
    } else if (grepl('\\.(txt|tsv)$',fp)) {
        mdat <- read.table(fp, sep="\t", header=T, check.names=F, row.names=NULL,
                           colClasses = "character", na.strings = naStr)
    } else {
        cat("Ignoring",fp,"\n")
        return(NULL)
    }
    if (!all(dim(mdat) > c(1,1))) {
        stop("What kind of table is ",fp,"?  As loaded it has ~no rows/cols.\n")
    }
    mdat <- HackHandleProblemSamples(mdat)
    colnames(mdat) <- c('SAMPLEID', colnames(mdat)[-1])
    mdat
}


## Load all metadata.  As when we prepared the OSM abstract, the first column
## should be a SAMPLEID and must not be duplicated within any file.  If we
## detect duplication, then probably it is not a SAMPLEID (e.g. it could be
## station number, date, ...), so reject the entire file.
metaList <- lapply(sidTab$filepath, function(fname) {
    cat("Reading",fname,"\n")
    x <- ReadMetaTable(fname)
    if (!is.null(x)) {
        if (length(unique(x$SAMPLEID)) < nrow(x)) {
            cat("REJECTING",fname,"\n")
            cat("\tThe first column of", fname, "\n",
                "\tis supposed to function as a SAMPLEID.\n")
            cat("\tHowever, there are", nrow(x),"rows in the table but only",
                length(unique(x$SAMPLEID)), "distinct values in the first column.\n")
            cat("\tPlease correct the metadata file so that its first column can function\n")
            cat("\tas a sample ID.  A sample ID may appear at most once per metadata file.\n")
            x <- NULL
        }
    }
    x
})
names(metaList) <- sidTab$filepath
idx <- which(!sapply(metaList, is.null))
cat("\nWill use",length(idx),"metadata tables of the",length(metaList),
    "in studyID_metadata.csv\n")
metaList <- metaList[idx]
sidTab <- sidTab[idx,]
stopifnot(sidTab$filepath == names(metaList))
rm(idx)
##DEBUG: Check that each table's SAMPLEID's look okay.
##lapply(metaList, function(m) m[1,'SAMPLEID'])

cat("\n\nSaving metatranscriptomic_samples.txt which has the unmodified sample",
    "IDs for metatranscriptomic samples.\n")
writeLines(sort(gMetaT), "metatranscriptomic_samples.txt")

## grep for the fields we want, looping over the tables just loaded.
x <- lapply(names(metaList),
    function(fnam) {
        mdat <- metaList[[fnam]]
        studyid <- sidTab$studyid[match(fnam,sidTab$filepath)]
        cnams <- sapply(wantCols, function(re) grep(re, colnames(mdat), ignore.case=T))
        idx <- which(sapply(cnams, length) > 1)
        if (length(idx) > 0) {
            cat("Multiple columns match regexps for", paste(names(cnams)[idx]),
                "in", fnam,"\n")
            stop("Aborting because I don't know which column to use.")
        }
        idx <- cnams[which(sapply(cnams, length) == 1)]  # Unambiguous columns in metadata
        mdat <- data.frame(SAMPLEID=mdat$SAMPLEID, StudyID=studyid, mdat[,as.numeric(idx)])
        colnames(mdat) <- c('SAMPLEID','StudyID', names(idx))
        if ('Collection_Date' %in% colnames(mdat)) {
            mdat$Collection_Date <- Fixup_Collection_Date(mdat$Collection_Date)
        }
        mdat
    })
names(x) <- names(metaList)
metaList <- x
rm(x)


## Merge all tables. For each SAMPLEID, add all fields seen across all tables.
## This allows info for a given SAMPLEID to be in multiple input tables.
##   Caveat: If a SAMPLEID appears with and without _transcriptomic, we assume
##   there two sequencing runs.  Some studies might have multiple tables, one
##   that labels a sample as genomic (in LibrarySource) and the other which
##   labels it as transcriptomic.  We cannot detect this annotation error from
##   real cases where a study has transcriptomic and genomic sequencing runs
##   that use the same sample ID.
metaTab <- metaList[[1]]
if (length(metaList) >= 2) {
    for (i in 2:length(metaList)) {
        cat("Merging",names(metaList)[i],"...\n")
        tmp <- merge(x = metaTab, y = metaList[[i]], by = 'SAMPLEID', all = TRUE,
                     suffixes = c(':LEFT',':RIGHT'))
        ## Collapse shared columns, with no tolerance for conflicting values from
        ## the two tables.
        cnams <- colnames(tmp)
        shared <- unique(sub(':(LEFT|RIGHT)$','', cnams[grep(':(LEFT|RIGHT)$',cnams)]))
        for (sc in shared) {
            ##cat("    Reconciling",sc,"\n")
            lvals <- tmp[,paste0(sc,':LEFT')]
            rvals <- tmp[,paste0(sc,':RIGHT')]
            ## It L and R have values for a sample, the values must be identical (or NA).
            if (!all(lvals==rvals, na.rm=T)) {
                idx <- which(lvals != rvals)
                cat(paste0("!!! ",length(idx)," conflicts when merging ", sc,".\n",
                           "!!! Here the first few conflicts. The right column is from ",
                           names(metaList)[i],".\n"))
                df <- data.frame(SAMPLEID   = tmp$SAMPLEID[idx],
                                 alreadySaw = lvals[idx],
                                 nowSeeing  = rvals[idx])
                print(df[1:min(nrow(df),10),])
                stop("Aborting due to conflicting metadata values.")
            }
            vals <- lvals
            rdef <- !is.na(rvals)
            vals[rdef] <- rvals[rdef]
            tmp <- tmp[, -grep(paste0("^",sc,":(LEFT|RIGHT)$"),colnames(tmp))]
            tmp[,sc] <- vals
        }
        metaTab <- tmp
        stopifnot(!grepl(":(LEFT|RIGHT)$",colnames(metaTab)))
    }
}


## This should reasonablly order the rows. (Run probably renamed to SAMPLEID.)
idx <- order(metaTab$SAMPLEID, metaTab$StudyID, metaTab$Collection_Date)
metaTab <- metaTab[idx,]
## Also order the cols
cnams <- c('SAMPLEID','StudyID', names(wantCols))
cnams <- cnams[cnams %in% colnames(metaTab)]
metaTab <- metaTab[,cnams]


cat("\nMissing data check:  For each column here are the percentages of NA's\n",
    "in each of the", nrow(metaTab), "total rows of the metadata table.\n")
print(round(100*apply(metaTab, 2, function(v) sum(is.na(v)))/nrow(metaTab)))

write.table(metaTab, 'metadata.tsv', sep="\t", row.names=F)
cat("\nWrote out metadata.tsv which has",nrow(metaTab),"samples [rows]",
    "and",ncol(metaTab),"fields. Done!\n\n")
quit('no')
