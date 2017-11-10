set -ev

git clone https://github.com/wdecoster/nanotest.git

NanoComp -h


NanoComp --bam alignment.bam alignment.bam alignment.bam --outdir compare-bams -f pdf

NanoComp --fastq reads.fastq.gz reads.fastq.gz --names run1 run2 --plot box
