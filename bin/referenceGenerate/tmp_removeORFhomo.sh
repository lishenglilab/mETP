#!/bin/bash
for i in $(seq 1 10)
do
{
    blastp -query ''/home/yvzeng/Project/mouseDRS/script/nextflow/v1/[:]/work/1a/87ef667b3b4abf836cf989ed6fde36/split/ORF/ORF_no_redundancy'_'${i}'.fa' \
           -out ''/home/yvzeng/Project/mouseDRS/script/nextflow/v1/[:]/work/1a/87ef667b3b4abf836cf989ed6fde36/split/blast/ORF_no_redundancy'_'${i}'.blast' \
           -db '/home/yvzeng/Project/mouseDRS/script/nextflow/v1/[:]/work/1a/87ef667b3b4abf836cf989ed6fde36/db/ref.db' \
           -outfmt 6 -evalue 1e-5 -num_threads 1
}&
done
wait