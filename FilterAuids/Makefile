##
## First, check the AUIDs for possible contaminants and compare the AUIDs
## against positive and negative nifH examples in GenBank as of 2020, classified
## by ARBitrator (Heller et al. 2014).
##
## Then filter out AUIDs using scripts/filterAuids.R:
## - Discard AUIDs that are too rare (see AUID_MIN_SAMPS, AUID_MIN_READS, and
##   AUID_ONE_HIT_WONDER_READS in the script).
## - Discard AUIDs that are outside the length range defined in this makefile by
##   MINASVLEN and MAXASVLEN.
## - Discard samples with too few reads (MIN_SAMPLE_SIZE).
## - Discard non-nifH AUIDs (RETAIN_NON_NIFH_LIKE = F).
## - Assess but do not discard possible contaminants (RETAIN_CONTAMINANTS = T).
##
## This stage depends on specific outputs from GatherAsvs (noted just below).
## Recall that uchime3 denovo filtering is done in GatherAsvs.
##
## Usage:
##    make                # Default, described above.
##    make quality        # Only do the quality checks: contaminants, nifH-like.
##    make clean          # Remove main targets but not contaminant and nifH-like checks.
##    make superclean     # Remove main and quality targets
##

SHELL=/usr/bin/bash

## Required outputs from GatherAsvs stage
INITIAL_AUID_ABUND_TAB=$(shell realpath ../GatherAsvs/asv2auid.abundances.tsv.gz)
INITIAL_FASTA=$(shell realpath ../GatherAsvs/asv2auid.fasta)

## Nucleotide length range to filter out ASVs that are unlikely to be nifH.  In
## practice we mainly see nifH ASVs that are ~325 nt (if amplicons are from
## nested PCR that used nifH2 and nifH1 for the second round [Zehr and
## McReynolds, 1989], and sequencing was done with MiSeq 2 x 300nt).  The
## 280-360 range helps filter out PCR artefacts.
MINASVLEN=280
MAXASVLEN=360

FINAL_AUID_ABUND_TAB=auid.abundances.filtered.tsv.gz
FINAL_FASTA=auid.filtered.fasta
LENFILTIMPACTS=lenFilterImpactsToSamps.$(MINASVLEN).$(MAXASVLEN).tsv
## Main targets.  (Also need quality, see 'all' target.)
ALLTARGETS=$(FINAL_AUID_ABUND_TAB) $(FINAL_FASTA) $(LENFILTIMPACTS)

## Quality: contaminant check and nifH-like check.  Scripts are provided via
## symlinks in 'bin' when pipeline is installed.
CHECKCONTAMS := $(shell which check_nifH_contaminants.sh || exit -1)
CONTAM_CHECK=Check_contaminants/blast.96id.out
CLASSIFYNIFH := $(shell which classifyNifH.sh || exit -1)
NIFH_LIKE_CHECK=Check_NifH_like


all: $(ALLTARGETS) quality

quality: $(CONTAM_CHECK) $(NIFH_LIKE_CHECK)

clean:
	rm -f $(ALLTARGETS)

## These take awhile to create.
superclean:  clean
	rm -rf Check_contaminants Check_NifH_like


$(FINAL_AUID_ABUND_TAB) $(LENFILTIMPACTS) &: $(INITIAL_AUID_ABUND_TAB) $(INITIAL_FASTA) \
					     $(CONTAM_CHECK) $(NIFH_LIKE_CHECK)
	@echo ">>>>> Filtering AUID abundance table <<<<<"
	@scripts/filterAuids.R $(INITIAL_AUID_ABUND_TAB) \
	  $(INITIAL_FASTA) $(MINASVLEN) $(MAXASVLEN) \
	  $(CONTAM_CHECK) $(NIFH_LIKE_CHECK)
	@echo ">>>>> Finished filtering AUID abundance table <<<<<"
	@echo


$(FINAL_FASTA): $(INITIAL_FASTA) $(FINAL_AUID_ABUND_TAB)
	@echo ">>>>> Filtering FASTA based on filtered AUID abundance table <<<<<"
	@cat $(FINAL_AUID_ABUND_TAB) | gunzip | sed "1d" | cut -f1 > ids.tmp ; \
	extractFasta.pl ids.tmp $(INITIAL_FASTA) > $(FINAL_FASTA) ; \
	echo "Expected to extract" `cat ids.tmp | wc -l` "AUID sequences" ; \
	echo "Exctracted" `cat $(FINAL_FASTA) | grep -c '^>'` "AUID sequences." ; \
	rm ids.tmp
	@echo ">>>>> Finished filtering FASTA <<<<<"
	@echo


$(CONTAM_CHECK): $(INITIAL_FASTA)
	@echo ">>>>> Checking AUIDs for known contaminants (amino %id > 96) <<<<<"
	@mkdir Check_contaminants.tmp ; \
	cd Check_contaminants.tmp ; \
	$(CHECKCONTAMS) $(INITIAL_FASTA) 96 &> log.contaminants.txt
	@mv Check_contaminants.tmp Check_contaminants
	@echo ">>>>> Finished contaminants check <<<<<"
	@echo


## This one takes a few hours.
$(NIFH_LIKE_CHECK): $(INITIAL_FASTA)
	@echo ">>>>> Comparing ASVs to positive/negative nifH examples in GenBank <<<<<"
	@echo "This step compares the ASVS to positive and negative nifH examples classified by a 2020"
	@echo "run of ARBitrator (Heller et al. 2014). Filtering ASVs based on the relative strengths"
	@echo "of their hits to positive vs. negative examples is done later by filterAuids.R."
	@mkdir Check_NifH_like.tmp ; \
	cd Check_NifH_like.tmp ; \
	$(CLASSIFYNIFH) $(INITIAL_FASTA) &> log.classifyNifH.txt
	@mv Check_NifH_like.tmp Check_NifH_like
	@echo ">>>>> Finished ASV comparison to GenBank nifH <<<<<"
	@echo
