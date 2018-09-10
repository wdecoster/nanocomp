set -ev

git clone https://github.com/wdecoster/nanotest.git

NanoComp -h
NanoComp --bam nanotest/alignment.bam nanotest/alignment.bam nanotest/alignment.bam --outdir compare-bams -f pdf -n run1 run2 run3
NanoComp --summary nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt -n A B C
NanoComp --fastq nanotest/reads.fastq.gz nanotest/reads.fastq.gz --names run1 run2 --plot box
NanoComp --fasta nanotest/reads.fa.gz nanotest/reads.fa.gz nanotest/reads.fa.gz nanotest/reads.fa.gz --n run1 run2 run3 run4 --plot violin
NanoComp --summary nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt -n A B --maxlength 20000
