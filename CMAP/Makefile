## Copyright (C) 2023 Mo Morando
##
## FIXME: Assuming just one CMAP script called query_CMAP.py that takes as input only
## metadata.cmap.tsv and outputs CMAP_data.csv.gz.  (Or can gzip in rule below.)  Assuming
## CMAP_data.csv.gz is formatted similarly to 20220831CMAP_colocalized_nifH.csv.
##
## Use ../GatherMetadata/metadata.cmap.tsv to query CMAP for environmental data for each sample,
## using the sample's date (at local noon), latitude, longitude, and depth.  Produce CMAP_data.csv
## which has one row for each sample and columns that are the CMAP query fields and all the
## downloaded/colocated environmental variables.
##
## Usage:
##    make                # Default. Get CMAP data if metadata.cmap.tsv changed.
##    make clean          # Remove CMAP_data.csv.gz
##

SHELL=/usr/bin/bash

QUERY_TABLE=$(shell realpath ../GatherMetadata/metadata.cmap.tsv)
CMAP_TABLE=CMAP_data.csv.gz

$(CMAP_TABLE) : $(QUERY_TABLE)
	@echo ">>>>> For each sample getting colocated environmental data from CMAP <<<<<"
	@scripts/query_CMAP.py $(QUERY_TABLE)
	@echo ">>>>> Finished downloading CMAP data <<<<<"
	@echo

clean:
	rm -f $(CMAP_TABLE)
