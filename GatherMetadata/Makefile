## Copyright (C) 2023 Jonathan D. Magasin
##
## Gather and standardize sample metadata across all studies.  Sample metadata
## mainly means sequencing run tables obtained from the SRA.  Scrape for just
## the fields we need.
##
## Outputs:
##  - studyID_metadata.csv     Concatenated metadata tables from all studies
##
##  - metadata.tsv             Cleanup up version of previous.  Just the fields
##                             needed, standardized and quality checked.  Some
##                             studies use the same sample ID for their DNA and
##                             cDNA sequencing runs. To ensure unique sample IDs
##                             in this table, "_transcriptomic" is appended to
##                             the ID for all metatranscriptomic samples.
##
##  - metatranscriptomic_samples.txt   Transcriptomic sample IDs without
##                             "_transcriptomic" appended as in metadata.tsv.
##
##  - samples_need_CMAP.txt    Samples that need CMAP environmental data.  This
##                             is created from the column names in the filtered
##                             abundance table from FilterAuids.
##
##  - metadata.cmap.tsv        Modified version of metadata.tsv with fields
##                             formatted for use as CMAP queries.
##
## Usage:
##    make                # Make the four outputs
##    make clean          # Remove all outputs
##

SHELL=/usr/bin/bash

targets=studyID_metadata.csv metadata.tsv metadata.cmap.tsv metatranscriptomic_samples.txt \
        samples_need_CMAP.txt

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


## Get a list (one file per line) of the samples in the abundance table.  We
## want only these samples in metadata.cmap.tsv so that the CMAP stage only asks
## for environmental data for samples with ASV data.
## Note that sed is used grab the sample ID by stripping of the tag preceding ___.
samples_need_CMAP.txt: ../FilterAuids/auid.abundances.filtered.tsv.gz
	@echo ">>>>> Determining which samples have ASVs and thus will need CMAP data. <<<<<"
	@zcat $< | head -n 1 | tr "\t" "\n" | sed -e 's/^..*___//' > $@
	@echo "There are" `cat $@ | wc -l` "samples that will need CMAP data."
	@echo ">>>>> Finished determining samples that need CMAP data. <<<<<"
	@echo

metadata.cmap.tsv: metadata.tsv samples_need_CMAP.txt
	@echo ">>>>> Standardizing metadata values when possible (dates, lat/lon) <<<<<"
	@echo "This step makes $@"
	@scripts/prepareMetadataForCmap.R
	@echo ">>>>> Finished standardizing metadata values <<<<<"
	@echo
