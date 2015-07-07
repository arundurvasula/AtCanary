#!/bin/bash -l
#SBATCH -D /home/lv70590/Arun/AtCanary/scripts
#SBATCH -o /home/lv70590/Arun/AtCanary/results/%j-out.txt
#SBATCH -e /home/lv70590/Arun/AtCanary/results/%j-err.txt
#SBATCH -J psmc
#SBATCH -p mem_0064
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.durvasula@gmail.com
set -e
set -u

HOME=/home/lv70590/Arun
SW=${HOME}/software

# run the psmc pipeline on an individual from /home/lv70590/Andrea/data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.bam

${SW}/samtools/samtools mpileup -C50 -uf /home/lv70590/Andrea/runningShore/reference/atAndSlociColl02.fa /home/lv70590/Andrea/data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.bam | ${SW}/bcftools-1.2/bcftools view -c - \
      | ${SW}/bcftools-1.2/vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ../data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.fq.gz

${SW}/psmc/utils/fq2psmcfa -q20 ../data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.fq.gz > ../data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.psmcfa
${SW}/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ../data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.psmc ../data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.psmcfa
${SW}/psmc/utils/psmc2history.pl ../data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.psmc | ${SW}/psmc/utils/history2ms.pl > ms-cmd.sh
${SW}/psmc/utils/psmc_plot.pl ../results/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111 ../data/16509_CTGTAGCC_C2RJUACXX_6_20140111B_20140111.psmc
