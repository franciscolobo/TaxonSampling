#!/bin/sh

echo "Downloading NCBI Taxonomy"

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

mkdir taxdump

mv new_taxdump.tar.gz taxdump/

tar -xvzf taxdump/new_taxdump.tar.gz

mv *dmp taxdump/

rm new_taxdump.tar.gz

