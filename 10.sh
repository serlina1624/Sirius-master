#!/bin/bash
# 10. Mafft
mkdir 10_mafft
for i in 9_filter/UCE*.fasta
do
cp $i 10_mafft
done
cd 10_mafft 
for i in UCE*.fasta
do 
mafft $i > mafft_$i
done
python ../10.py
