set -ev

git clone https://github.com/wdecoster/nanotest.git

NanoComp -h


NanoComp --bam nanotest/alignment.bam nanotest/alignment.bam nanotest/alignment.bam --outdir compare-bams -f pdf -n run1 run2 run3

NanoComp --fastq nanotest/reads.fastq.gz nanotest/reads.fastq.gz --names run1 run2 --plot box
