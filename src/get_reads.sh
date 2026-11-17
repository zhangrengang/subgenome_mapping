#!/bin/bash
# Download the reads of every design line from ENA (mirror of SRA), as gzipped fastq.
#
# usage: bash src/get_reads.sh [species.design] [outdir=reads] [ncpu=4]
#
# design.txt col1 is the run accession (SRR/ERR/DRR...); col2/col3 are not used here.
# Only the runs whose files are missing are downloaded, so re-running is cheap.
# Output: <outdir>/<run>_1.fastq.gz and <outdir>/<run>_2.fastq.gz
#   (single-end runs: only <run>_1.fastq.gz)
set -eu

design=${1:-species.design}
outdir=${2:-reads}
ncpu=${3:-4}
mkdir -p $outdir

while read -r run sp ref
do
	[ -n "${run:-}" ] || continue
	if [ -s "$outdir/${run}_1.fastq.gz" ] && [ -s "$outdir/${run}_2.fastq.gz" ]; then
		echo "== $run already downloaded, skip"
		continue
	fi
	echo "== $run ($sp): asking ENA for the fastq links"
	urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv" \
		| awk -F'\t' 'NR>1 && $2!="" {print $2}' | tr ';' '\n')
	[ -n "$urls" ] || { echo "ERROR: ENA returned no fastq for $run"; exit 1; }
	i=0
	for u in $urls
	do
		i=$((i+1))
		case $i in 1) s=1 ;; 2) s=2 ;; *) s=$i ;; esac
		echo "   download $(basename "$u") -> $outdir/${run}_${s}.fastq.gz"
		curl -s -L --retry 5 -C - -o "$outdir/${run}_${s}.fastq.gz" "https://${u}"
	done
	echo "   done: $(ls -la $outdir/${run}_*.fastq.gz | awk '{printf "%.0f MB ", $5/1048576}')"
done < $design

echo
echo "READS READY: $outdir"
