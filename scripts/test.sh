set -ev


if [ -d "nanotest" ]; then
  echo "nanotest already cloned"
else
  git clone https://github.com/wdecoster/nanotest.git
fi

NanoComp -h
NanoComp --bam nanotest/alignment.bam nanotest/alignment.bam nanotest/alignment.bam -f pdf -n run1 run2 run3 -o tests
NanoComp --summary nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt -n A B C -o tests
NanoComp --fastq nanotest/reads.fastq.gz nanotest/reads.fastq.gz --names run1 run2 --plot box -o tests
NanoComp --fastq_rich nanotest/reads.fastq.gz nanotest/reads.fastq.gz --names run1 run2 --plot box -o tests
NanoComp --fasta nanotest/reads.fa.gz nanotest/reads.fa.gz nanotest/reads.fa.gz nanotest/reads.fa.gz --n run1 run2 run3 run4 --plot violin -o tests
NanoComp --summary nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt -n A B --maxlength 20000 -o tests
NanoComp --fastq nanotest/reads.fastq.gz nanotest/reads.fastq.gz --colors red blue -n run1 run2 --format svg --plot ridge --outdir tests
NanoComp --feather nanotest/summary1.feather nanotest/summary2.feather nanotest/summary3.feather -n A B C --outdir tests
NanoComp --feather nanotest/summary1.feather nanotest/summary2.feather nanotest/summary3.feather -n A B C --outdir tests -f pdf png jpeg
