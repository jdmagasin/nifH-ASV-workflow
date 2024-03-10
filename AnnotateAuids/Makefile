## Copyright (C) 2023 Jonathan D. Magasin
##
## Annotate ASVs from stage FilterAuids using three tools that are included in the
## DADA2 nifH pipeline's "scripts.ancillary/Annotation" directory:
##  1. BlastnARB2017          - Blastn ASVs against the 2017 ARB nifH database 
##  2. BlastxGenome879        - Blastx ASVs against the Zehr Lab's database with
##                              879 diazotroph genomes.
##  3. NifHClustersFrank2016  - Determine nifH clusters for each ASV using the
##                              classification and regression tree (CART) 
##                              approach of Frank et al. 2016.
##
##  4. Marine_NCD_cyano_refsequences
##                            - 223 nifH sequences from selected
##                              cyanobacteria and NCDs to provide increased
##                              resolution.
##  5. UCYNA_oligoreps        - Forty-four UCNY-A oligotypes, thirty-eight
##                              from metagenome assemblies, to provide increased
##                              resolution of UCYN-A ASVs.
##
##
## This Makefile runs each of the annotation tools, then searches the cyano /
## NCD DB, the UCYN-A oligos DB, and then stitches together all the results 
## into auids.annot.tsv, using scripts/makeAnnotationTable.R. A consensus ID is 
## assigned for each AUID using information gathered across all the merged 
## databases using scripts/make_consensus_taxonomy.py and then outputs the 
## annotation table, auids.annot.consensus.tsv. This table is the final output 
## of the Makefile.
##
##
## Usage:
##    make                # Generate all annotation types.  
##                          Make auids.annot.consensus.tsv and 
##                          auids.annot.consensus.tsv
##    make ARB2017        # Do just the blastn against ARB 2017
##    make Genome879      # Do just the blastx against Genome879
##    make Clusters       # Do just the nifH cluster classifications
##    make NCD_cyano      # Do just the blastn against 267 NCD/cyano ref nifH
##    make UCYNA_oligos   # Do just the blastn against UCYN-A nifh oligotypes

##
##    make clean          # Remove auids.annot.tsv and auids.annot.consensus.tsv
##    make superclean     # Remove auids.annot.tsv, auids.annot.consensus.tsv,
##                          and all annotation output dirs.
##

SHELL=/usr/bin/bash

## Input: The fasta from FilterAuids
AUID_FASTA=$(shell realpath ../FilterAuids/auid.filtered.fasta)

## Output (main):  Table with all annottion
ANNOT_TAB=auids.annot.tsv

## Output (consensus taxonomy): Table with consensus taxonomy ID added
CONS_TAXONOMY_TAB=auids.annot.consensus.tsv

## The tools. All are sym-links in the pipeline's 'bin'
BLASTNARB2017=$(shell which blastnARB2017.sh)
BLASTXGENOME879=$(shell which blastxGenome879.sh)
ASSIGNNIFHCLUSTERS=$(shell which assignNifHclustersToNuclSeqs.sh)


## ARB2017
ARB2017_DIR=BlastnARB2017
ARB2017_LOG=log.BlastnARB2017.txt
ARB2017_PCTID=70
ARB2017_QRYCVG=90
ARB2017_DESC=$(ARB2017_PCTID)id.$(ARB2017_QRYCVG)cvg
ARB2017_OUTS=$(addprefix $(ARB2017_DIR)/, blastnNifH.$(ARB2017_DESC).tab blastnNifH.$(ARB2017_DESC).aln.gz)
ARB2017_TAB=$(word 1, $(ARB2017_OUTS))


## Genome879
GENOME879_DIR=BlastxGenome879
GENOME879_LOG=log.BlastxGenome879.txt
GENOME879_PCTID=50
GENOME879_QRYCVG=90
GENOME879_DESC=$(GENOME879_PCTID)pid.$(GENOME879_QRYCVG)cvg
GENOME879_OUTS=$(addprefix $(GENOME879_DIR)/, blastxGenome879.$(GENOME879_DESC).tab blastxGenome879.$(GENOME879_DESC).aln.gz)
GENOME879_TAB=$(word 1, $(GENOME879_OUTS))


## Marine NCD and cyanobacteria nifH reference sequences (~267 total)
NCDCYANO_DIR=Marine_NCD_cyano_nifH
NCDCYANO_LOG=log.marine_NCD_cyano.txt
NCDCYANO_FASTA=$(shell realpath data/Marine_NCD_cyano_nifH/Marine_NCD_cyano_refsequences.fasta.gz)
NCDCYANO_PCTID=97
NCDCYANO_QRYCVG=70
NCDCYANO_DESC=$(NCDCYANO_PCTID)pid.$(NCDCYANO_QRYCVG)cvg
NCDCYANO_TAB=$(addprefix $(NCDCYANO_DIR)/, blastn_NCD_cyano.$(NCDCYANO_DESC).tab)
NCDCYANO_ALN=$(addprefix $(NCDCYANO_DIR)/, blastn_NCD_cyano.$(NCDCYANO_DESC).aln.gz)
NCDCYANO_OUTS=$(NCDCYANO_TAB) $(NCDCYANO_ALN)


## UCYN-A nifH reference sequences (~267 total)
UCYNA_OLIGOS_DIR=UCYNA_oligoreps
UCYNA_OLIGOS_LOG=log.OLIGOS.txt
UCYNA_OLIGOS_FASTA=$(shell realpath data/UCYNA_oligoreps/UCYNA_oligoreps.fasta.gz)
UCYNA_OLIGOS_PCTID=97
UCYNA_OLIGOS_QRYCVG=70
UCYNA_OLIGOS_DESC=$(UCYNA_OLIGOS_PCTID)pid.$(UCYNA_OLIGOS_QRYCVG)cvg
UCYNA_OLIGOS_TAB=$(addprefix $(UCYNA_OLIGOS_DIR)/, blastn_OLIGOS.$(UCYNA_OLIGOS_DESC).tab)
UCYNA_OLIGOS_ALN=$(addprefix $(UCYNA_OLIGOS_DIR)/, blastn_OLIGOS.$(UCYNA_OLIGOS_DESC).aln.gz)
UCYNA_OLIGOS_OUTS=$(UCYNA_OLIGOS_TAB) $(UCYNA_OLIGOS_ALN)


## NifH clusters
NIFHCLUST_DIR=Classify_NifH_Clusters
NIFHCLUST_LOG=log.nifH_clusters.txt
NIFHCLUST_OUTS=$(addprefix $(NIFHCLUST_DIR)/, auid.filtered.clusters.map auid.filtered.withNifH_Clusters.fasta)
NIFHCLUST_MAP=$(word 1, $(NIFHCLUST_OUTS))


## Marine NCD and cyanobacteria nifH reference sequences (~267 total)
all: $(CONS_TAXONOMY_TAB)

clean:
	@echo "Removing $(ANNOT_TAB) $(CONS_TAXONOMY_TAB)"
	@if [ -f "$(ANNOT_TAB)" ] ; then rm "$(ANNOT_TAB)" ; fi
	@if [ -f "$(CONS_TAXONOMY_TAB)" ] ; then rm "$(CONS_TAXONOMY_TAB)" ; fi

superclean:
	@echo "Removing $(ANNOT_TAB), $(CONS_TAXONOMY_TAB), and all annotation working directories."
	@if [ -f "$(ANNOT_TAB)" ] ; then rm "$(ANNOT_TAB)" ; fi
	@if [ -f "$(CONS_TAXONOMY_TAB)" ] ; then rm "$(CONS_TAXONOMY_TAB)" ; fi
	@for x in $(ARB2017_DIR) $(GENOME879_DIR) $(NIFHCLUST_DIR) $(NCDCYANO_DIR) $(UCYNA_OLIGOS_DIR); do \
	    if [ -d "$$x" ] ; then rm -r "$$x"; fi \
	done

ARB2017: $(ARB2017_OUTS)

Genome879: $(GENOME879_OUTS)

Clusters : $(NIFHCLUST_OUTS)

NCD_cyano: $(NCDCYANO_OUTS)

UCYNA_oligos: $(UCYNA_OLIGOS_OUTS)


## Note the use of the cluster .map file.
$(ANNOT_TAB): $(ARB2017_OUTS) $(GENOME879_OUTS) $(NCDCYANO_OUTS) $(UCYNA_OLIGOS_OUTS) $(NIFHCLUST_OUTS)
	@echo ">>>>> Combining all annotation into one table $(ANNOT_TAB) <<<<<"
	@scripts/makeAnnotationTable.R $(ARB2017_TAB) $(GENOME879_TAB) $(NCDCYANO_TAB) $(UCYNA_OLIGOS_TAB) $(NIFHCLUST_MAP)
	@echo ">>>>> Finished combining annotation <<<<<"
	@echo


$(ARB2017_OUTS) &: $(AUID_FASTA)
	@echo ">>>>> Searching AUIDs against ARB2017 [blastn, $(ARB2017_PCTID) %nid, $(ARB2017_QRYCVG) % coverage] <<<<<"
	@echo "Once completed, see log file in $(ARB2017_DIR)."
	@echo "Main outputs are $(ARB2017_OUTS)"
	@$(BLASTNARB2017) $(AUID_FASTA) $(ARB2017_DIR) $(ARB2017_PCTID) $(ARB2017_QRYCVG) &> $(ARB2017_LOG) ; \
	mv $(ARB2017_LOG) $(ARB2017_DIR)/
	@echo ">>>>> Done searching ARB2017 <<<<<"
	@echo


$(GENOME879_OUTS) &: $(AUID_FASTA)
	@echo ">>>>> Searching AUIDs against Genome879 <<<<<"
	@echo "Once completed, see log file in $(GENOME879_DIR)."
	@echo "Main outputs are $(GENOME879_OUTS)"
	@$(BLASTXGENOME879) $(AUID_FASTA) $(GENOME879_DIR) $(GENOME879_PCTID) $(GENOME879_QRYCVG) &> $(GENOME879_LOG) ; \
	mv $(GENOME879_LOG) $(GENOME879_DIR)/
	@echo ">>>>> Done searching Genome879 <<<<<"
	@echo


$(NIFHCLUST_OUTS) &: $(AUID_FASTA)
	@echo ">>>>> Assigning nifH clusters to the AUIDs (method: Frank et al. 2106) <<<<<"
	@echo "Once completed, see log file in $(NIFHCLUST_DIR) in particular table at end."
	@echo "Main outputs are $(NIFHCLUST_OUTS)"
	@mkdir -p $(NIFHCLUST_DIR)_tmp ; \
	cd $(NIFHCLUST_DIR)_tmp ; \
	$(ASSIGNNIFHCLUSTERS) $(AUID_FASTA) &> $(NIFHCLUST_LOG) ; \
	cd .. ; \
	mv $(NIFHCLUST_DIR)_tmp $(NIFHCLUST_DIR)
	@echo ">>>>> Finished nifH cluster assignments <<<<<"
	@echo


## Make the DB, then search and save as ASN.1. Then create tabular and alignment outputs.
$(NCDCYANO_OUTS) &: $(AUID_FASTA) $(NCDCYANO_FASTA)
	@echo ">>>>> Searching AUIDs against additional reference nifH from cyanos and NCDs <<<<<"
	@mkdir -p $(NCDCYANO_DIR) ; \
	cd $(NCDCYANO_DIR) ; \
	echo "Making database" ; \
	cat $(NCDCYANO_FASTA) | gunzip \
	  | makeblastdb -dbtype nucl -title NCD_cyano_nifH -out NCD_cyano_nifH ; \
	echo "Searching database" ; \
	blastn -db NCD_cyano_nifH -query $(AUID_FASTA) \
	    -perc_identity $(NCDCYANO_PCTID) -qcov_hsp_perc $(NCDCYANO_QRYCVG) -outfmt 11 > asn.tmp ; \
	echo "Reformatting output" ; \
	blast_formatter -archive asn.tmp -outfmt 0 -num_alignments 10 -num_descriptions 100 \
	  | gzip > `basename $(NCDCYANO_ALN)` ; \
	blast_formatter -archive asn.tmp -outfmt 6 > `basename $(NCDCYANO_TAB)` ; \
	echo "Cleaning up" ; \
	rm asn.tmp NCD_cyano_nifH.n??
	@echo ">>>>> Finished searching additional reference nifH <<<<<"
	@echo


# Make the DB, then search and save as ASN.1. Then create tabular and alignment outputs.
$(UCYNA_OLIGOS_OUTS) &: $(AUID_FASTA) $(UCYNA_OLIGOS_FASTA)
	@echo ">>>>> Searching AUIDs against additional reference nifH from UCYN-A oligos <<<<<"
	@mkdir -p $(UCYNA_OLIGOS_DIR) ; \
	cd $(UCYNA_OLIGOS_DIR) ; \
	echo "Making database" ; \
	cat $(UCYNA_OLIGOS_FASTA) | gunzip \
	  | makeblastdb -dbtype nucl -title UCYNA_oligos_nifH -out UCYNA_oligos_nifH ; \
	echo "Searching database" ; \
	blastn -db UCYNA_oligos_nifH -query $(AUID_FASTA) \
	    -perc_identity $(UCYNA_OLIGOS_PCTID) -qcov_hsp_perc $(UCYNA_OLIGOS_QRYCVG) -outfmt 11 > asn.tmp ; \
	echo "Reformatting output" ; \
	blast_formatter -archive asn.tmp -outfmt 0 -num_alignments 10 -num_descriptions 100 \
	  | gzip > `basename $(UCYNA_OLIGOS_ALN)` ; \
	blast_formatter -archive asn.tmp -outfmt 6 > `basename $(UCYNA_OLIGOS_TAB)` ; \
	echo "Cleaning up" ; \
	rm asn.tmp UCYNA_oligos_nifH.n??
	@echo ">>>>> Finished searching additional reference nifH <<<<<"
	@echo

## Consensus taxonomy table
$(CONS_TAXONOMY_TAB): $(ANNOT_TAB)
	@echo ">>>>> Creating consensus taxonomy table $(CONS_TAXONOMY_TAB) <<<<<"
	@python scripts/make_consensus_taxonomy.py --annotation_table $(ANNOT_TAB) --output_table $(CONS_TAXONOMY_TAB)
	@echo ">>>>> Finished creating consensus taxonomy table <<<<<"
	@echo