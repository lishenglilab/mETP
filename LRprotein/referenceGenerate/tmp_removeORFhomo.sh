#!/bin/bash
for i in `seq 1 10`
do
{
blastp -query /home/yvzeng/Project/PASM/data/protein/reference/split/ORF/ORF_no_redundancy_${i}.fa -out /home/yvzeng/Project/PASM/data/protein/reference/split/blast/ORF_no_redundancy_${i}.blast -db /home/yvzeng/Project/PASM/data/protein/reference/db/ref.db -outfmt 6 -evalue 1e-5 -num_threads 1
}&
done
wait
