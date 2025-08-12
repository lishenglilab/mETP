#!/bin/bash
for i in `seq 1 10`
do
{
python /home/yvzeng/Project/mouseDRS/script/nextflow/v1/bin/referenceGenerate/1.ORF_generate.py -i '/home/yvzeng/Project/mouseDRS/script/nextflow/v1/[:]/work/0a/1377c06783946691ea5b4bd107a25b/split/trancript/isoforms'_${i}.fa -o '/home/yvzeng/Project/mouseDRS/script/nextflow/v1/[:]/work/0a/1377c06783946691ea5b4bd107a25b/split/ORF/' -p ORF_${i}
}&
done
wait