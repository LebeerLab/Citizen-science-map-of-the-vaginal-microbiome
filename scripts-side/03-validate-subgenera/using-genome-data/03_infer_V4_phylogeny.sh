# This script aligns complete V4 sequences and infers a phylogeny. 

# dependencies: MAFFT v7.453, IQTREE v1.6.12

fin_V4=../../results/validation/V4_complete.fna
fout_aln=../../results/validation/V4_complete.aln
dout_tree=../../results/validation/tree

# create output folder 
[ -d $dout_tree ] || mkdir $dout_tree

# align the sequences
mafft $fin_V4 > $fout_aln

# infer phylogeny 
iqtree -s $fout_aln -m GTR+G4 -pre $dout_tree/V4_complete