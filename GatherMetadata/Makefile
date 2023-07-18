##
## Gather and standardize sample metadata across all studies.  Sample metadata
## mainly means sequencing run tables obtained from the SRA.  Scrape for just
## the fields we need.
##
## Outputs:
##  - studyID_metadata.csv     Concatenated metadata tables from all studies
##
##  - metadata.tsv             Cleanup up version of previous.  Just the fields
##                             needed, standardized and quality checked.  Only
##                             metagenomic samples retained.
##
##  - metatranscriptomic_samples.txt   Transcriptomic sample IDs (original).
##                             In metadata.tsv the IDs have _transcriptomic
##                             appended to them to ensure that every sample
##                             ID in the table is unique.
##
##  - metadata.cmap.tsv        Modified version of metadata.tsv with fields
##                             formatted for use as CMAP queries.
##
## Usage:
##    make                # Make the four outputs
##    make clean          # Remove all outputs
##

SHELL=/usr/bin/bash

targets=studyID_metadata.csv metadata.tsv metadata.cmap.tsv metatranscriptomic_samples.txt

all: $(targets)

clean:
	@echo "Removing $(targets)"
	@rm -f $(targets)


## Use the curated/corrected metadata
studyID_metadata.csv: ../Data/StudyMetadata/*
	@echo ">>>>> Making a list of all the sample metadata tables <<<<<"
	@echo "This step makes $@"
	@find ../Data/StudyMetadata -type f > allMeta.tmp ; \
	rm -f $@ ; \
	while read path ; do  \
	    study="$$(basename $$(dirname $$path))" ; \
	    echo "$${study},$${path}" >> $@ ; \
	done < allMeta.tmp ; \
	rm allMeta.tmp
	@echo ">>>>> Finished listing sample metadata tables <<<<<"
	@echo


## gatherMetadata.R requires unique sample identifiers in the first column of a
## metadata file, or it will throw out the metadata file.  It also throws out
## metadata files that have multiple rows for a sample id. It drops
## transcriptomic samples (saving IDs to a file) and makes some value tweaks to
## avoid conflicts when merging tables.  SEE WARNING MESSAGES written to stdout
## by the script.
metadata.tsv metatranscriptomic_samples.txt: studyID_metadata.csv
	@echo ">>>>> Merging sample metadata tables <<<<<"
	@echo "This step makes $@"
	@scripts/gatherMetadata.R
	@echo ; \
	if [ -f ../GatherAsvs/asv2auid.abundances.tsv.gz ] ; then \
	    echo '>>>>> Check mapping of AUID samples to metadata <<<<<' ; \
	    scripts/checkMapAuidTab2MetaTab.R ; \
	fi
	@echo ">>>>> Finished merging sample metadata tables <<<<<"
	@echo


metadata.cmap.tsv: metadata.tsv
	@echo ">>>>> Standardizing metadata values when possible (dates, lat/lon) <<<<<"
	@echo "This step makes $@"
	@scripts/prepareMetadataForCmap.R
	@echo ">>>>> Finished standardizing metadata values <<<<<"
	@echo
