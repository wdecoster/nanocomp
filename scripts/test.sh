#!/bin/bash
set -ev

# Parse command line arguments
RUN_ALL=false
if [ "$1" == "--all" ]; then
    RUN_ALL=true
fi

if [ -d "nanotest" ]; then
  echo "nanotest already cloned"
else
  git clone https://github.com/wdecoster/nanotest.git
fi

# Always run help to verify installation
NanoComp -h

# Fast test (default): only test BAM file
echo "Running quick BAM test..."
NanoComp --bam nanotest/alignment.bam nanotest/alignment.bam nanotest/alignment.bam -n run1 run2 run3 -o tests

# Full test suite (only run with --all flag)
if [ "$RUN_ALL" = true ]; then
    echo "Running full test suite..."
    NanoComp --summary nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt -n A B C -o tests
    NanoComp --fastq nanotest/reads.fastq.gz nanotest/reads.fastq.gz --names run1 run2 --plot box -o tests
    NanoComp --fastq_rich nanotest/reads.fastq.gz nanotest/reads.fastq.gz --names run1 run2 --plot box -o tests
    NanoComp --fasta nanotest/reads.fa.gz nanotest/reads.fa.gz nanotest/reads.fa.gz nanotest/reads.fa.gz --n run1 run2 run3 run4 --plot violin -o tests
    NanoComp --summary nanotest/sequencing_summary.txt nanotest/sequencing_summary.txt -n A B --maxlength 20000 -o tests
    NanoComp --fastq nanotest/reads.fastq.gz nanotest/reads.fastq.gz --colors red blue -n run1 run2 --format svg --plot ridge --outdir tests
    NanoComp --feather nanotest/summary1.feather nanotest/summary2.feather nanotest/summary3.feather -n A B C --outdir tests
    NanoComp --arrow nanotest/summary1.feather nanotest/summary2.feather nanotest/summary3.feather -n A B C --outdir tests
    #NanoComp --feather nanotest/summary1.feather nanotest/summary2.feather nanotest/summary3.feather -n A B C --outdir tests -f pdf png jpeg
    NanoComp --feather nanotest/summary1.feather nanotest/summary2.feather nanotest/summary3.feather -n A B C --color red blue green --outdir tests
    NanoComp --feather nanotest/summary1.feather nanotest/summary2.feather nanotest/summary3.feather --format json -n A B C --outdir tests
else
    echo "Quick test completed. Run with --all flag for full test suite."
fi