#!/bin/bash
for i in `seq 1 number`
do
{
blastp -query ORF_${i}.fa -out path_${i}.blast -db ref.db -outfmt 6 -evalue 1e-5 -num_threads 1
}&
done
wait
