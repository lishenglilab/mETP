Nanopore Direct RNA Sequencing (DRS) Analysis Pipeline
=
Overview
-
This Nextflow pipeline provides a comprehensive analysis solution for Nanopore Direct RNA Sequencing (DRS) data. The workflow performs end-to-end processing including quality control, transcript assembly, quantification, modification detection, alternative splicing analysis, and optional proteomics integration.
Requirements
-
Nextflow

Python：
subprocess
sys
os
argparse
datetime
pyfastx
pysam
re
pandas

R:
getopt
tidyverse
diann
plyr
ComplexHeatmap
ggsci
Rtsne
Hmisc
Biostrings
stringr
rtracklayer

other tools:
tombo
NanoPsu
flair
blast
gffread
suppa2
porechop
samtools
fragpipe(optional)

Installation
git clone https://github.com/yourusername/nanopore-drs-analysis.git
cd nanopore-drs-analysis
