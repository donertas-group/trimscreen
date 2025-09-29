#!/bin/sh

# Uses preformatted databases from DADA2 (https://benjjneb.github.io/dada2/training.html)
# Note that The DADA2 authors maintain and recommend pre-formatted training sets and species assignment files (hosted on https://benjjneb.github.io/dada2/training.html)
# SILVA 138.2 uses nr99 taxonomy
# The file for taxonomy assignment, identified by containing "train" in the name
gunzip -c *train*gz > assignTaxonomy.fna

# and the file for add species, identified by containing "species" in the name, is renamed
mv *assign*gz addSpecies.fna.gz
