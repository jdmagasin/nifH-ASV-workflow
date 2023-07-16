##
## Combine results from previous runs of the DADA2 nifH pipeline: Rename the ASVs to "AUIDs" (to
## resolve ASV ID collisions across pipeline runs) and combine the abundance tables.  Before
## combining the tables, remove chimera identified using uchime3 denovo.
##
## The 'quality' target forces the makefile to identify possible contaminants and to classify
## sequences as NifH or not by comparison to positive and negative GenBank sequences found by
## ARBitrator (Heller et al. 2015).
##
## Use of this makefile requires both the DADA2 nifH pipeline as well as the ancillary tools.
##
## Usage:
##    make                # Default.  Make AUID fasta and abundance table.
##    make clean          # Remove all outputs from this stage.
##

SHELL=/usr/bin/bash

## Handmade table of ASV files from pipeline runs.  See documentation in file.
TABLE_OF_ASV_FASTAS := asvs.noChimera.fasta_table.tsv

## From assignAUIDs2ASVs.R. Excludes asv2auid.abundances.tsv.gz which is created by
## makeBigAbundTable.R.
ASV2AUID_TSV_FASTA := asv2auid.tsv asv2auid.fasta


## Default targets.
all: $(ASV2AUID_TSV_FASTA) asv2auid.abundances.tsv.gz

clean:
	rm -f $(ASV2AUID_TSV_FASTA) asv2auid.abundances.tsv.gz
	rm -rf NonChimericAsvs Workdir.tmp CountTables


## Separately run uchime3_denovo on the ASVs from all the pipeline runs.  The
## output directory is NonChimericAsvs.
NonChimericAsvs: $(TABLE_OF_ASV_FASTAS)
	@echo "Running chimera filter (post-DADA2 nifH pipeline)"
	scripts/identifyChimeras.sh


## This target sets up for the next one, $(ASV2AUID_TSV_FASTA), by creating the
## working directory and getting the absolute paths to the ASV fasta files.
Workdir.tmp Workdir.tmp/absFastas.tsv : $(TABLE_OF_ASV_FASTAS)
	@echo Setting up Workdir.tmp which will be used for AUID assignments.
	@mkdir Workdir.tmp
	@cat $(TABLE_OF_ASV_FASTAS) | egrep -v '(^#|^$$)' \
	  | tr -s "\t" > Workdir.tmp/tsv.tmp
	@rm -f Workdir/absFastas.tsv
	@while read line ; do \
	  dsnam=`echo "$$line" | cut -f1`; \
	  dstag=`echo "$$line" | cut -f3 | tr -d "\n"` ; \
	  fasPath=`echo "$$line" | cut -f2`; \
	  if [ ! -f "$$fasPath" ] ; then \
	    echo "Cannot find $$fasPath" ; \
	    exit -1 ; \
	  fi ; \
	  fasPath=`realpath $$fasPath` ; \
	  echo -e "$${dsnam}\t$${fasPath}\t$${dstag}" >>  Workdir.tmp/absFastas.tsv ; \
	done < Workdir.tmp/tsv.tmp
	@rm Workdir.tmp/tsv.tmp


## Run assignAUIDs2ASVs.R on the chimera filtered FASTAs to create
## asv2auid.tsv. This takes a little set up because the required abundance
## tables are with the original ASV FASTAs but the de-chimera'd FASTAs are in
## NonChimericAsvs.  The while loop gathers the pairs and makes
## listOfAsvFastas.txt.  Copy the asv2auid.* but leave the Workdir with all
## symlinks because needed for making AUID count tables.
$(ASV2AUID_TSV_FASTA): Workdir.tmp/absFastas.tsv NonChimericAsvs
	@cd Workdir.tmp; \
	echo "Setting up listOfAsvFastas.txt." ; \
	echo "" > listOfAsvFastas.txt ; \
	while read line ; do \
	  dsnam=`echo "$$line" | cut -f1`; \
	  dstag=`echo "$$line" | cut -f3 | tr -d "\n"` ; \
	  dsnam="$${dsnam}_tag_$${dstag}" ; \
	  tsvpath=`echo "$$line" | cut -f2 | sed 's/\.fasta$$/.tsv/'` ; \
	  if [ ! -f "$$tsvpath" ] ; then \
	      echo "Cannot find $$tsvpath" ; \
	      exit -1 ; \
	  fi ; \
	  ncfasta="$${dsnam}.noChimera.fasta" ; \
	  abundTab="$${dsnam}.noChimera.tsv" ; \
	  ln -s "../NonChimericAsvs/$$ncfasta" . ; \
	  ln -s "$$tsvpath" "./$$abundTab" ; \
	  if [ `cat "$$abundTab" | wc -l` -eq 0 ] ; then \
	      echo "No rows in $$abundTab" ; \
	      exit -1 ; \
	  fi ; \
	  echo -e "$${ncfasta}\t$${dstag}" >> listOfAsvFastas.txt ; \
	done < absFastas.tsv ; \
	echo "Assigning AUIDs to the non-chimeric ASVs." ; \
	echo "Restricting to ASVs in length range [200,450]" ; \
	assignAUIDs2ASVs.R listOfAsvFastas.txt -fasta -length=200,450 ; \
	echo ; \
	echo "See NonChimericAsvs/chimera.log for results summaries for each study.\n" ; \
	echo
	@cp Workdir.tmp/asv2auid.* .


CountTables: $(ASV2AUID_TSV_FASTA) Workdir.tmp
	@echo "Making CountTables that use the AUIDs."
	@cd Workdir.tmp ; \
	makeAUIDCountTable.R asv2auid.tsv > log.counts.txt ; \
	mkdir CountTables ; \
	mv counts.*.tsv log.counts.txt CountTables
	@mv Workdir.tmp/CountTables .


asv2auid.abundances.tsv.gz: CountTables
	@echo "Combining CountTables into one big AUID abundance table $@"
	@scripts/makeBigAbundTable.R ; \
	gzip asv2auid.abundances.tsv
