#!/bin/bash
for i in `seq 1 10`
do
{
python /home/yvzeng/Project/software/proteomics/script/referenceGenerate/1.ORF_generate.py -i /home/yvzeng/Project/KDP/mouseTAC/data/RNAseq/proteinReference/split/trancript/mouseTACtranscriptExpressedAutosome_${i}.fa -o /home/yvzeng/Project/KDP/mouseTAC/data/RNAseq/proteinReference/split/ORF/ -p ORF_${i}
}&
done
wait