#!/bin/bash

set -e
set -u

sample1="$1"
sample2="$2"
out="$3"

ref="data/atAndSlociColl02.fa"
HOME="/home/CIBIV/arun/"
SW="${HOME}/software"

samtools mpileup -C50 -uf ${ref} ${sample1} | bcftools view -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > data/${sample1}.fq.gz

samtools mpileup -C50 -uf ${ref} ${sample2} | bcftools view -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > data/${sample2}.fq.gz

seqtk mergefa -hq20 ${sample1}.fq.gz ${sample2}.fq.gz \
      | ${SW}/psmc/utils/fq2psmcfa -q30 - > data/${out}.psmcfa

${SW}/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o data/${out}.psmc data/${out}.psmcfa
${SW}/psmc/utils/psmc2history.pl data/${out}.psmc | ${SW}/psmc/utils/history2ms.pl > ms-cmd.sh
${SW}/psmc/utils/psmc_plot.pl ../results/${out} data/${out}.psmc
