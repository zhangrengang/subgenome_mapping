#!/bin/bash
# Trim the reads of every design line with fastp.
#
# usage: bash src/trim_reads.sh [species.design] [readdir=reads] [ncpu=8]
#
# Input : <readdir>/<run>_1.fastq.gz, <readdir>/<run>_2.fastq.gz (single-end: only _1)
# Output: <readdir>/<run>_clean_1.fastq.gz, <run>_clean_2.fastq.gz, <run>.fastp.html
#   (the two routes map/assemble the *_clean_* files)
# Options are the ones used for the project this pipeline comes from:
#   cut the first/last bases while quality <= 3, drop reads with any N,
#   require >= 60 bp after trimming, correct low-quality bases in overlapping pairs.
set -eu

design=${1:-species.design}
readdir=${2:-reads}
ncpu=${3:-8}

while read -r run sp ref dtype
do
	[ -n "${run:-}" ] || continue
	r1=$readdir/${run}_1.fastq.gz
	r2=$readdir/${run}_2.fastq.gz
	o1=$readdir/${run}_clean_1.fastq.gz
	o2=$readdir/${run}_clean_2.fastq.gz
	[ -s "$r1" ] || { echo "ERROR: no $r1 (run src/get_reads.sh first)"; exit 1; }
	if [ -s "$o1" ]; then
		echo "== $run already trimmed, skip"
		continue
	fi
	echo "== $run ($sp): fastp"
	if [ -s "$r2" ]; then
		fastp --in1 "$r1" --in2 "$r2" --out1 "$o1" --out2 "$o2" \
			--html $readdir/${run}.fastp.html \
			--cut_by_quality3 --cut_by_quality5 --n_base_limit 1 \
			--length_required 60 --correction --thread $ncpu
	else
		fastp --in1 "$r1" --out1 "$o1" \
			--html $readdir/${run}.fastp.html \
			--cut_by_quality3 --cut_by_quality5 --n_base_limit 1 \
			--length_required 60 --thread $ncpu
	fi
	echo "   done: $(ls -la $o1 $o2 2>/dev/null | awk '{printf "%.0f MB ", $5/1048576}')"
done < $design

echo
echo "TRIMMED READS READY: $readdir (<run>_clean_{1,2}.fastq.gz)"
