#!/usr/bin/env bash 

# This script will create various custom 16S rRNA reference databases. They are
# all based on the GTDB taxonomy and 16S sequences, and they are all in the 
# format readable by DADA2.

fin_changes=../results/lactobacillus_subgenera/changes.csv
fin_db_bac_rep=../data/bac120_ssu_reps_r95.tar.gz
fin_db_bac_all=../data/ssu_all_r95.tar.gz
fout_db_bac_rep=../results/SSUrRNA_GTDB05-bac-rep_DADA2.fna
fout_db_bac_all=../results/SSUrRNA_GTDB05-bac-all_DADA2.fna
fout_db_lac_rep=../results/SSUrRNA_GTDB05-lactobacillaceae-rep_DADA2.fna
fout_db_lac_all=../results/SSUrRNA_GTDB05-lactobacillaceae-all_DADA2.fna
fout_db_subgenera_rep=../results/SSUrRNA_GTDB05-lactobacillus-subgenera-rep_DADA2.fna
fout_db_subgenera_all=../results/SSUrRNA_GTDB05-lactobacillus-subgenera-all_DADA2.fna

# create refdb with representative sequences from all species
tar -xzf  $fin_db_bac_rep --directory $( dirname $fin_db_bac_rep )
./reformat_refdb.py ${fin_db_bac_rep%%.tar.gz}.fna $fout_db_bac_rep
rm ${fin_db_bac_rep%%.tar.gz}.fna

# create refdb with all sequences from all species
tar -xzf  $fin_db_bac_all --directory $( dirname $fin_db_bac_all )
./reformat_refdb.py ${fin_db_bac_all%%.tar.gz}.fna $fout_db_bac_all
rm ${fin_db_bac_all%%.tar.gz}.fna

# create refdb with representative sequences from Lactobacillaceae
./subset_refdb.py $fout_db_bac_rep $fout_db_lac_rep 4 Lactobacillaceae

# create refdb with all sequences from Lactobacillaceae
./subset_refdb.py $fout_db_bac_all $fout_db_lac_all 4 Lactobacillaceae

# create refdb with subgenera of Lactobacillus (representative sequences)
./subset_refdb.py $fout_db_lac_rep $fout_db_subgenera_rep.tmp 5 Lactobacillus
./adapt_refdb.py $fout_db_subgenera_rep.tmp $fout_db_subgenera_rep $fin_changes
rm $fout_db_subgenera_rep.tmp

# create refdb with subgenera of Lactobacillus (all sequences)
./subset_refdb.py $fout_db_lac_all $fout_db_subgenera_all.tmp 5 Lactobacillus
./adapt_refdb.py $fout_db_subgenera_all.tmp $fout_db_subgenera_all $fin_changes
rm $fout_db_subgenera_all.tmp