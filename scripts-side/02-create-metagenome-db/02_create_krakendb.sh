#!/usr/bin/env bash 

# This script will create a custom kraken2 database. It can easily be adapted
# to create extra custom databases. 

# dependencies: GNU parallel, kraken2

tax_from_gtdb=../bin/Metagenomics-Index-Correction/tax_from_gtdb.py
fin_taxonomy=../results/bac120_taxonomy_r95_adapted.tsv
fin_genomes=../data/gtdb_genomes_reps_r95.tar.gz
fin_humantaxonomy=../data/GRCh38_taxonomy.tsv
fin_humangenome=../data/GRCh38_latest_genomic.fna.gz
dout_krakendb=../results/shotgun_16G_GTDB05-bac-rep_GRCh38_lactosubgenera

threads=16

# function to create a krakendb
# args:
# - fin_taxonomy: A tab-separated file with two columns: the name of the genome
#   assembly and the taxonomy (taxa of all ranks, separated by semicolons).
# - din_genomes: A folder with a genome assembly (fna file) for each line in
#   fin_taxonomy; the filename should be structured as 
#   <genomename>_genomic.fna.gz; the tag "GB_" or "RS_" at the beginning of 
#   the genome name can be omitted.
# - dout_db: Output folder to store the kraken2 database in. 
create_krakendb() {
  
  local fin_taxonomy=$1
  local din_genomes=$2
  local dout_db=$3

  # create output folders
  [ -d $dout_db ] || mkdir $dout_db
  [ -d $dout_db/library ] || mkdir $dout_db/library
  [ -d $dout_db/taxonomy ] || mkdir $dout_db/taxonomy 
  
  # STEP 1: convert GTDB db (taxonomy + genomes) to NCBI-formatted db
  # (takes longer than 30 minutes with 16 threads)
  $tax_from_gtdb --gtdb $fin_taxonomy \
    --assemblies $din_genomes \
    --nodes $dout_db/taxonomy/nodes.dmp \
    --names $dout_db/taxonomy/names.dmp \
    --kraken_dir $dout_db/library
  
  # STEP 2: create preliminary seqid/taxid maps
  # (takes ~ 15 minutes with 16 threads)
  find $dout_db/library/ -name *.fa |\
    parallel \
    --jobs $threads \
    --no-notice \
    --verbose \
    kraken2-build --add-to-library {} --db $dout_db \
    2>&1 | tee $dout_db/log_create_taxid_maps.txt
  
  # STEP 3: build the actual kraken2 database 
  # (takes ~ 13 minutes with 16 threads)
  kraken2-build --build --db $dout_db --threads $threads \
    --max-db-size 16000000000 2>&1 | tee $dout_db/log_build_db.txt
    
  # remove large intermediate files
  if [ -f $dout_db/hash.k2d ] ; then
    rm -r $dout_db/library $dout_db/taxonomy
  fi
  
}

# untar bacterial genomes
tar -xzf $fin_genomes --directory $( dirname $fin_genomes )

# add human genome (fna) to bacterial genomes (fnas)
cp $fin_humangenome ${fin_genomes%.tar.gz}

# add human genome taxonomy to GTDB taxonomy file
cp $fin_taxonomy $fin_taxonomy.withhuman
cat $fin_humantaxonomy >> $fin_taxonomy.withhuman

# create krakendb 
create_krakendb $fin_taxonomy.withhuman ${fin_genomes%.tar.gz} $dout_krakendb

# do some cleanup (if db creation successful)
if [ -f $dout_krakendb/hash.k2d ] ; then
  rm -r ${fin_genomes%.tar.gz}
fi
