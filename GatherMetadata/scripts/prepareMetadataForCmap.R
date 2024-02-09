#!/usr/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

##
## Prepare metadata so that it can be used to query CMAP.
##  - Standardize / correct some fields:  Collection_Date*
##
## *Actually Collection_Date gets fixed up by gatherMetadata.R. If we do not do
##  it there, then two rows for a given sample will appear to be in conflict due
##  to differently formatted Collection_Date's.
##

metadata <- read.table('metadata.tsv', header=T, check.names=F, stringsAsFactors=T)
sampsNeeded <- unique(readLines('samples_need_CMAP.txt'))

## gatherMetadata.R tacks on "_transcriptomic" to RNA sample names that will appear
## in metadata.tsv, but _transcriptomic is not in the abundance table column names
## and thus it is not in sampsNeeded. Drop _transcriptomic for the following check.
simpleSampID <- sub('_transcriptomic$','',metadata$SAMPLEID)
x <- setdiff(sampsNeeded, simpleSampID)
if (length(x) > 0) {
    warning("There is no metadata for the following ",length(x), " samples so CMAP ",
            "data will not be obtained: ", paste(x, collapse=','))
}
metadata <- metadata[simpleSampID %in% sampsNeeded, ]
stopifnot(nrow(metadata) > 0)  # bad samples_need_CMAP.txt?

## These are the fields that are carefully checked. (A few more are checked only
## to see that they always are NA.)
fieldsChecked <- c('Collection_Date','Depth','Lat_Lon','Latitude','Longitude')
cat("Will carefully check the following metadata fields:\n",
    "\t",paste(fieldsChecked, collapse=', '),"\n")

Announce <- function(field) {
    cat("\n\n")
    cat("##--------------------------------\n")
    cat("## Checking", field,"\n")
    cat("##\n")
}


needsFixin <- list()


##--------------------------------
##
## Collection_Date: Copied this from gatherAsvs.R. Don't really need to *fix*
## Collection_Dates but this script would seem incomplete if it did not try (and
## if fixup is ever taken out of gatherAsvs.R...)
##
## However, LATER (once we have lat and lon) we do need to convert
## to UTC and specific format for CMAP:
##     https://cmap.readthedocs.io/en/latest/user_guide/API_ref/pycmap_api/data_retrieval/pycmap_query.html
##
Announce('Collection_Date')

## Fixup dates to be in format e.g. 2021-1-17.  Some examples of the formats
## fixed up:   1/7/21  1.7.21  1/07/21  1/17/2021
## Not bullet proof.  Also, does not make 2-digit month or day as required
## by CMAP, but the UTC conversion happens to do that.
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


coldates <- Fixup_Collection_Date(as.character(metadata$Collection_Date))
## Inelegant.  Use reg.ymd from the function to report failed conversions.
reg.ymd <- '20[0-9]{2}-[0-9]{1,2}-[0-9]{1,2}'
idx <- grep(reg.ymd, coldates, invert=T)
if (length(idx) > 0) {
    cat(length(idx), "of the",length(coldates), "dates still need to be fixed.\n")
    cat("Here are unique date values that could not be converted to format YYYY-MM-DD:\n")
    x <- unique(coldates[idx])
    cat("\t",paste(x,collapse=', '),"\n")
    cat("These dates are associated with the following studies:\n")
    cat("\t",
        paste(as.character(unique(subset(metadata, Collection_Date %in% x)$StudyID)),
              collapse=', '),
        "\n")
    rm(x)
} else {
    cat("Good news -- all Collection_Dates are in format YYYY-MM-DD")
}

## Commit the corrections.
stopifnot(length(coldates) == nrow(metadata))
metadata$Collection_Date <- factor(coldates)


## Verify that what needs fixing is as determined above. Then record
## the rows in our "needs fixin'" mask.
stopifnot(idx == grep(reg.ymd, metadata$Collection_Date, invert=T))
needsFixin$Collection_Date <- idx
rm(idx,coldates)


##--------------------------------
##
## Depth
##
Announce('Depth')
cat("Fixing up depths.  Removing trailing 'm' but otherwise leaving alone.\n")
x <- sub(' *m *$', '', as.character(metadata$Depth))
metadata$Depth <- factor(x)

x <- suppressWarnings(as.numeric(x))
cat(sum(is.na(x)), "Depths have problems that need manual fixing.\n")
cat("Here are the bad Depth values:\n")
cat("\t",paste(as.character(unique(metadata[is.na(x),'Depth'])), collapse=', '),"\n")
cat("\nThose values appear in the following studies:\n")
cat("\t",paste(as.character(unique(metadata[is.na(x),'StudyID'])), collapse=", "),"\n")

## Save problem cases.
needsFixin$Depth <- which(is.na(x))


##--------------------------------
##
## Lat_Lon
##   Insist on form "36.251130 S 153.263230 E" or "36.251130 N 153.263230 W"
##   because that's what I mostly see e.g. from the SRA.
##   However, CMAP requires degrees East and North:
##       https://cmap.readthedocs.io/en/latest/faq_and_contributing/file_structure.html#lat
##
Announce('Lat_Lon')
cat("Checking Lat_Lon format. Although CMAP requires separate lat and lon\n",
    "in degrees N and E, we can use Lat_Lon to make those fields.\n")
reg.lat_lon <- paste0('^[[:digit:]]+\\.*[[:digit:]]*',' +[NS]+ +',
                       '[[:digit:]]+\\.*[[:digit:]]*',' +[WE]+$')
latlon <- as.character(metadata$Lat_Lon)
idx <- grep(reg.lat_lon, latlon, invert=T)
cat("Found",length(idx),"rows with Lat_Lon that does not follow the\n",
    "format \"<posNum> [NS] <posNum> [WE]\".  These are the misformatted values:\n")
cat("\t",paste(unique(latlon[idx]),collapse=', '),"\n\n")
cat("For some studies the fix is simply to add the compass direction after verifying\n",
    "that the degrees are correct for the study region.\n")
needsFixin$Lat_Lon <- idx


##--------------------------------
##
## Latitude and Longitude
##   Use these if available and Lat_Lon is missing.
##
Announce('Latitude and Longitude')
Make_Alt_Lat_Lon <- function()
{
    alt_lat_lon <- lapply(c('Latitude','Longitude'),
        function(v) {
            x <- as.character(metadata[,v])
            stopifnot(length(x) == nrow(metadata) && length(x) > 0)
            {
                ## Mo added some coords with NSWE so handle that.
                swreg <- '^[\\-]{0,1}[0-9]+\\.{0,1}[0-9]* *[SW]$'
                idx <- grep(swreg,x)  # Needs sign flip
                if (length(idx) > 0) {
                    annie <- as.numeric(sub(' *[SW]$','',x[idx]))
                    stopifnot(!is.na(annie)) # Fail if anything couldn't be coverted
                    x[idx] <- as.character(0 - annie)
                }
                ## Remove the directions for those that did not need a sign flip.
                nereg <- '^[\\-]{0,1}[0-9]+\\.{0,1}[0-9]* *[NE]$'
                idx <- grep(nereg,x)
                if (length(idx) > 0) { x[idx] <- sub(' *[NE]$','',x[idx]) }
            }
            idx <- which(!is.na(x))
            asNum <- suppressWarnings(as.numeric(x[idx]))  # NA's when cannot convert
            if (length(idx) > 0 && !is.numeric(asNum)) {
                cat("Found some",v,"that are are non-numeric. All must be",
                    "positive or negative decimal numbers.\n")
                x <- unique(x[idx][is.na(asNum)])
                cat("Some examples on num-numerics:", x[1:min(10, length(x))],"\n")
                stop("Aborting")
            }
            stopifnot(length(idx) == length(asNum))
            if (v == 'Latitude') {
                x[idx] <- paste(abs(asNum), c('N','S')[(asNum < 0)+1])  # N positive
            } else {
                x[idx] <- paste(abs(asNum), c('E','W')[(asNum < 0)+1])  # E positive
            }
            x[idx][is.na(asNum)] <- NA
            x
        })
    x <- paste(alt_lat_lon[[1]],alt_lat_lon[[2]])
    x[x=="NA NA"] <- NA
    stopifnot(length(x) == nrow(metadata))
    x
}
alt_lat_lon <- Make_Alt_Lat_Lon()

## Fixup *missing* values Lat_Lon (when we can, with alt_lat_lon). If there was
## a Lat_Lon specified but is has problems, this script will not try to figure
## out if the alt_lat_lon is a safe replacement.
latlon <- as.character(metadata$Lat_Lon)
idx <- needsFixin$Lat_Lon
idx <- idx[is.na(latlon[idx])]
if (length(idx) > 0) {
    cat("Will use", sum(table(alt_lat_lon)), "lat/lon coordinates constructed\n",
        "from Latitude and Longitude to set *missing* Lat_Lon.\n")
    stopifnot(is.na(latlon[idx]))
    stopifnot(length(latlon) == length(alt_lat_lon))
    latlon[idx] <- alt_lat_lon[idx]  # okay if overwrites NA with NA
    metadata$Lat_Lon <- factor(latlon)
    cat("Set",sum(!is.na(alt_lat_lon[idx])),
        "missing Lat_Lon using the corresponding Latitude and Longitude.\n")
    idx <- grep(reg.lat_lon, latlon, invert=T)
    needsFixin$Lat_Lon <- idx
    cat("There are now",length(idx),"rows that need Lat_Lon to be fixed.\n")
}


##--------------------------------
##
## CMAP lat and lon from Lat_Lon
##
Announce(' lat and lon for CMAP')

cat("This step takes the Lat_Lon from above and splits it into\n",
    "separate Lat and Lon columns with degrees N or E as required\n",
    "for CMAP queries.  If Lat or Lon columns already exist, will\n",
    "report any 'clobbers' that change a coordinate by > 0.001.\n")

for (x in c('Lat','Lon')) {
    if (x %in% colnames(metadata)) {
        idx <- which(!is.na(metadata[,x]) & !is.na(metadata$Lat_Lon))
        if (length(idx) > 0) {
            cat(x,"is defined",length(idx),"times in rows that have\n",
                "Lat_Lon defined. Will check in a moment whether",x,"\n",
                "will be clobbered by a substantially different value\n",
                "from Lat_Lon.\n")
        }
    }
}

## Only parse Lat_Lon with valid format!
vidx <- grep(reg.lat_lon, metadata$Lat_Lon)
x <- strsplit(as.character(metadata$Lat_Lon[vidx]), ' +')
lats <- as.numeric(sapply(x, '[[', 1)); stopifnot(!is.na(lats))
lats <- lats * c(+1,-1) [(sapply(x, '[[', 2) == 'S') +1  ]
lons <- as.numeric(sapply(x, '[[', 3)); stopifnot(!is.na(lons))
lons <- lons * c(+1,-1)[(sapply(x, '[[', 4) == 'W') +1]

## Helper to report clobbers.
check_clobber_coord <- function(oldvals,cnam)
{
    if (any(!is.na(oldvals))) {
        idx <- which(abs(oldvals - metadata[,cnam]) > 0.001)
        if (length(idx) > 0) {
            cat("The following",cnam,"'s got clobbered:\n")
            df <- data.frame(metadata$SAMPLEID[idx], oldals[idx],
                             metadata$Lat[idx])
            colnames(df) <- c('SAMPLEID',
                              paste0(cnam,'.original'),
                              paste0(cnam,'.Lat_Lon'))
            print(df)
        }
    }
}

## Create or overwrite Lat and Lon. Then check if clobbered.
oldvals <- NA
if ('Lat'  %in% colnames(metadata)) { oldvals <- as.numeric(metadata$Lat) }
metadata$Lat <- NA
metadata$Lat[vidx] <- lats
check_clobber_coord(oldvals,'Lat')

if ('Lon'  %in% colnames(metadata)) { oldvals <- as.numeric(metadata$Lon) }
metadata$Lon <- NA
metadata$Lon[vidx] <- lons
check_clobber_coord(oldvals,'Lon')


##--------------------------------
##
## Time
##   Not sure why this column is included if we are looking only at DNA.
##
Announce('Time')
if (all(is.na(metadata$Time))) {
    cat("The Time column is entirely NA.\n")
}  else {
    cat("Hmm, there are non-NA values in the Time column.\n",
        "Update script to check Time even if not use transcriptomic samples?.\n")
}


##--------------------------------
##
## UTC datetimes for local noon
##
Announce('Local noon --> UTC datetimes')

## Added this function after learning that CMAP requires UTC.
library(lutz)
MakeLocalNoonsAsUtc <- function(mtab)
{
    stopifnot(c('Lat','Lon','Collection_Date') %in% colnames(mtab))
    ## Get time zones (or NA).  The 'accurate' method requires additional
    ## libraries.  'fast' will cause a warning to be issued about possibly
    ## inaccurate tz's near unpopulated areas. However, for meta*g*enomic
    ## samples being off by 1 hr should have ~no impact on CMAP modeled data.
    ## For metatranscriptomic samples 'accurate' would be better.
    tzones <- tz_lookup_coords(mtab$Lat, mtab$Lon, method='fast')
    localNoons <- paste(mtab$Collection_Date,"12:00:00")
    ## Magic based on:
    ##  https://blog.revolutionanalytics.com/2009/06/converting-time-zones.html
    ## Nicely this also makes the months and days have 2 digits (leading 0).
    utcDateTimes <- sapply(1:length(localNoons), function (i) {
        x <- tzones[i]
        if (!is.na(x)) { x <- format(as.POSIXct(localNoons[i], x), tz="GMT") }
        x
    })
    utcDateTimes <- sub(' ','T',utcDateTimes)
    ## df <- data.frame(localNoons, tzones, utcDateTimes)  # for debugging
    utcDateTimes
}


cat("Adding column DateTimeUTC that has UTC time for noon at collection site and date.\n")
cat("This step converts from the timezone implied by lat and lon.\n")
cat("It tries to covert every Collection_Date.  If a Collection_Date needs fixing,\n",
    "I would not trust the conversion.\n")
metadata$DateTimeUTC <- MakeLocalNoonsAsUtc(metadata[,c('Lat','Lon','Collection_Date')])


##--------------------------------
##
## Flag rows with metadata that is missing (NA) or misformatted
##

## Flag rows that had a bad value other than NA (which is easy enough
## to just look for in the row itself).
tot <- c()
for (fc in fieldsChecked) {
    idx <- which(is.na(metadata[,fc]))
    idx <- setdiff( needsFixin[[fc]], idx )
    if (length(idx) > 0) {
        ## There are some actual bad values so add a FIXME column.
        metadata[,paste0('FIXME.',fc)] <- FALSE
        metadata[idx,paste0('FIXME.',fc)] <- TRUE
        tot <- union(tot,idx)
    }
}

cat("\n\nA total of", length(tot),"rows in the metadata table have\n",
    "one or more fields with a problem that needs manual fixing.  This\n",
    "excludes missing values (NA). See the rightmost \"FIXME\" columns\n",
    "to determine which rows/fields need attention.\n\n")

## Correct some of the field names to be what CMAP expects in queries:
##   https://cmap.readthedocs.io/en/latest/user_guide/API_ref/pycmap_api/data_retrieval/pycmap_query.html
## Simplest to do here e.g. in case CMAP ever changes the names.
cnams <- colnames(metadata)
map2qryfields <- c(Lat = 'lat', Lon = 'lon', Depth = 'depth', DateTimeUTC = 'time')
idx <- match(names(map2qryfields), cnams)
stopifnot(!is.na(idx))
cnams[idx] <- as.character(map2qryfields)
stopifnot(table(cnams) == 1)  # ensure no duplicated field namees
colnames(metadata) <- cnams

## Now restrict to just the fields CMAP needs to know (lat,lon,time,depth) and
## SAMPLEID, StudyID, and any of the "FIXME." fields.  The goal is to prevent
## CMAP from seeing/responding to any fields that we do not expect it to (due to
## API changes or bugs).
wanted <- c('SAMPLEID', 'StudyID', 'lat', 'lon', 'time', 'depth')
wanted <- c(wanted, colnames(metadata)[grep('^FIXME\\.',colnames(metadata))])
cat("Restricting CMAP query table metadata.cmap.tsv to these columns:\n", wanted, "\n")
stopifnot(wanted %in% colnames(metadata))
metadata <- metadata[,wanted]

cat("\nWriting metadata.cmap.tsv\n")
write.table(metadata, 'metadata.cmap.tsv', sep="\t", row.names=F)
cat("Done!\n")
quit('no')
