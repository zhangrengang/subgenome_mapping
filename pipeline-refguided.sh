#!/bin/bash
# Subgenome mapping -- reference-guided route.
#
# reads --> competitive mapping to EVERY copy of the reference genome
#       --> targeted assembly of each copy (HybPiper, reusing that mapping)
#       --> per-copy sequences appended as new columns of the subgenome matrix
#       --> subgenome gene trees
#
# usage: bash pipeline-refguided.sh [design.txt] [ncpu] [outgroup]
#
# design.txt: 3 columns, no header, no comments:
#   <sample_id>  <species_id>  <reference_species>
#     sample_id  run accession; reads/<sample_id>_{1,2}.fastq.gz (see src/get_reads.sh)
#     species_id label of this sample; its new columns are <species_id>_<copy>
#     reference  species whose copies are columns of pre-study/merged.alignment.csv
#
# requires on PATH: soi-wgdi, bwa, samtools; HybPiper via HYBPIPER (its bin dir) and
# PYTHON (the interpreter that install needs -- python2 for HybPiper 1.2)
set -eu

design=${1:-design.txt}
ncpu=${2:-8}
outgroup=${3:-Theobroma_cacao}
pre=pre-study
reads=reads
out=out
hyp=${HYBPIPER:-$(dirname "$(command -v reads_first.py || echo /nonexistent/reads_first.py)")}
py=${PYTHON:-python}
[ -s "$hyp/reads_first.py" ] || { echo "ERROR: reads_first.py not found; set HYBPIPER=<HybPiper bin dir>"; exit 1; }
mkdir -p $out

##########################################################################
##### per sample: map, assemble, append, tree ############################
##########################################################################
while read -r sample sp ref
do
	[ -n "${sample:-}" ] || continue
	echo
	echo "########## $sp  (reference: $ref) ##########"

	##################################################################
	##### 1. competitive mapping #####################################
	# every read is assigned to at most one copy of the reference, and the
	# resulting bam only keeps the genes whose copies are all retained; the
	# depth table says whether the copies are covered evenly
	##################################################################
	echo "== 1) compmap: $reads/${sample}_{1,2}.fastq.gz -> $ref copies"
	mkdir -p $out/$sp.$ref
	soi-wgdi compmap $pre/merged.alignment.csv $pre/cds.fa \
		$reads/${sample}_1.fastq.gz $reads/${sample}_2.fastq.gz \
		colnames=$pre/rename.list ref=$ref prefix=$out/$sp.$ref/$sp.$ref ncpu=$ncpu

	##################################################################
	##### 2. targeted assembly of each copy ##########################
	# HybPiper reuses the bam from step 1 (--bwa --no-blast): one spades +
	# exonerate run per reference gene, i.e. one sequence per copy
	##################################################################
	echo "== 2) HybPiper assembly of the retained copies"
	# HybPiper wants cwd = the directory that CONTAINS <prefix>/ and --prefix to be
	# the bare directory name: it looks for <prefix>/<prefix>.bam and
	# <prefix>/<prefix>.target.cds, both written by compmap above.
	# scripts are called through $py, so no chmod +x is needed on the HybPiper install
	( cd $out
	  touch dummy.fastq
	  $py $hyp/reads_first.py -r dummy.fastq -b $sp.$ref/$sp.$ref.target.cds \
		--prefix $sp.$ref --cpu $ncpu --bwa --no-blast
	  $py $hyp/paralog_investigator.py $sp.$ref
	)
	python3 src/retrieve_copies.py $out/$sp.$ref
	echo "   assemblies: $(grep -c '^>' $out/$sp.$ref.pep)"

	##################################################################
	##### 3. append the sample as new columns ########################
	# one column per copy of the reference: <species_id>_<copy>
	##################################################################
	echo "== 3) append_cols: $sp -> new columns"
	soi-wgdi append_cols $pre/merged.alignment.csv $out/$sp.$ref.cds $sp \
		colnames=$pre/rename.list prefix=$out/$sp.$ref-appended

	##################################################################
	##### 4. subgenome tree of this sample ###########################
	##################################################################
	echo "== 4) subgenome tree: $sp x $ref"
	ind=$(awk -F'\t' -v r="$ref" -v s="$sp" -v o="$outgroup" '
		{ n=$2; if (n ~ ("^" r "_") || n ~ ("^" s "_") || n == o) printf "%s%s", (k++ ? "," : ""), $1 }
		' $out/$sp.$ref-appended.rename.list)
	echo "   columns: $ind"
	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' $pre/pep.faa > $out/$sp.$ref-tree.pep.faa
	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' $pre/cds.fa > $out/$sp.$ref-tree.cds.fa
	cat $out/$sp.$ref.pep >> $out/$sp.$ref-tree.pep.faa
	cat $out/$sp.$ref.cds >> $out/$sp.$ref-tree.cds.fa
	bash src/subgenome-tree.sh $out/$sp.$ref-appended.alignment.csv $out/$sp.$ref-appended.rename.list \
		$out/$sp.$ref-tree.pep.faa $out/$sp.$ref-tree.cds.fa \
		$out/$sp.$ref-subgenome "$ind" $outgroup $ncpu
done < $design

##########################################################################
##### one tree per reference: every sample of that reference together ####
##########################################################################
for ref in $(awk '{print $3}' $design | sort -u)
do
	ns=$(awk -v r="$ref" '$3 == r' $design | wc -l)
	[ "$ns" -ge 2 ] || { echo "== combined tree: only $ns sample(s) for $ref, skipped"; continue; }
	echo
	echo "########## combined tree: $ref ##########"
	aln=$pre/merged.alignment.csv
	cols=$pre/rename.list
	k=0
	for sp in $(awk -v r="$ref" '$3 == r {print $2}' $design)
	do
		k=$((k+1))
		if [ $k -eq "$ns" ]; then
			ocsv=$out/$ref-combined.alignment.csv
			ocol=$out/$ref-combined.rename.list
		else
			ocsv=$out/$ref.pass$k.alignment.csv
			ocol=$out/$ref.pass$k.rename.list
		fi
		soi-wgdi append_cols $aln $out/$sp.$ref.cds $sp colnames=$cols prefix=${ocsv%.alignment.csv}
		aln=$ocsv
		cols=$ocol
	done
	rm -f $out/$ref.pass*.alignment.csv $out/$ref.pass*.rename.list
	ind=$(awk -F'\t' -v r="$ref" -v o="$outgroup" '
		{ n=$2; if (n ~ ("^" r "_") || n == o) printf "%s%s", (k++ ? "," : ""), $1 }
		' $cols)
	for sp in $(awk -v r="$ref" '$3 == r {print $2}' $design)
	do
		ind="$ind$(awk -F'\t' -v s="$sp" '$2 ~ ("^" s "_") {printf ",%s", $1}' $cols)"
	done
	echo "   columns: $ind"
	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' $pre/pep.faa > $out/$ref-tree.pep.faa
	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' $pre/cds.fa > $out/$ref-tree.cds.fa
	for sp in $(awk -v r="$ref" '$3 == r {print $2}' $design)
	do
		cat $out/$sp.$ref.pep >> $out/$ref-tree.pep.faa
		cat $out/$sp.$ref.cds >> $out/$ref-tree.cds.fa
	done
	bash src/subgenome-tree.sh $out/$ref-combined.alignment.csv $cols \
		$out/$ref-tree.pep.faa $out/$ref-tree.cds.fa \
		$out/$ref-subgenome "$ind" $outgroup $ncpu
done

echo
echo "ALL DONE"
