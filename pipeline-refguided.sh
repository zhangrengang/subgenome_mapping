#!/bin/bash
# Subgenome mapping -- reference-guided route.
#
# reads --> competitive mapping to EVERY copy of the reference genome
#       --> targeted assembly of each copy (HybPiper, reusing that mapping)
#       --> per-copy sequences appended as new columns of the subgenome matrix
#       --> subgenome gene trees
#
# usage: bash pipeline-refguided.sh [species.design] [ncpu] [outgroup]
#
# species.design: 4 columns, no header, no comments:
#   <run>  <species_id>  <reference_species>  <WGS|HybSeq>
#     run        accession; reads/<run>_clean_{1,2}.fastq.gz (src/get_reads.sh,
#                then src/trim_reads.sh)
#     species_id label of this sample; its new columns are <species_id>_<subgenome>
#     reference  species whose copies are columns of pre-study/merged.alignment.csv
#     dtype      library type, printed with the log line (both types work here)
#
# requires on PATH: soi-wgdi, bwa, samtools; HybPiper via HYBPIPER (its bin dir) and
# PYTHON (the interpreter that install needs -- python2 for HybPiper 1.2)
set -eu

design=${1:-species.design}
ncpu=${2:-8}
outgroup=${3:-Theobroma_cacao}
pre=pre-study
reads=reads
out=out-refguided
hyp=${HYBPIPER:-$(dirname "$(command -v reads_first.py || echo /nonexistent/reads_first.py)")}
py=${PYTHON:-python}
[ -s "$hyp/reads_first.py" ] || { echo "ERROR: reads_first.py not found; set HYBPIPER=<HybPiper bin dir>"; exit 1; }
mkdir -p $out

##########################################################################
##### per sample: map, assemble, append, tree ############################
##########################################################################
while read -r sample sp ref dtype
do
	[ -n "${sample:-}" ] || continue
	echo
	echo "########## $sp  (reference: $ref, $dtype) ##########"

	##################################################################
	##### 1. competitive mapping #####################################
	# every read is assigned to at most one copy of the reference, and the
	# resulting bam only keeps the genes whose copies are all retained; the
	# depth table says whether the copies are covered evenly
	##################################################################
	echo "== 1) compmap: $reads/${sample}_clean_{1,2}.fastq.gz -> $ref copies"
	mkdir -p $out/$sp.$ref
	r1=$reads/${sample}_clean_1.fastq.gz
	r2=$reads/${sample}_clean_2.fastq.gz
	[ -s "$r2" ] || r2=          # single-end: compmap takes one read file
	soi-wgdi compmap $pre/merged.alignment.csv $pre/cds.fa \
		$r1 $r2 \
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

done < $design

##########################################################################
##### one tree per reference: every sample of that reference together ####
##########################################################################
ref=$(awk '{print $3}' $design | sort -u)
[ $(echo "$ref" | wc -l) -eq 1 ] || { echo "ERROR: $(echo "$ref" | wc -l) references in $design; the matrix is one reference, split the design"; exit 1; }
ns=$(awk -v r="$ref" '$3 == r' $design | wc -l)
echo
echo "########## subgenome tree: $ref ($ns sample(s)) ##########"
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
	soi-wgdi baitmap $aln $out/$sp.$ref.cds $sp colnames=$cols prefix=${ocsv%.alignment.csv}
	aln=$ocsv
	cols=$ocol
done
rm -f $out/$ref.pass*.alignment.csv $out/$ref.pass*.rename.list
ind=$(seq -s, 2 $(wc -l < $cols))
echo "   columns: $ind"
awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' $pre/pep.faa > $out/$ref-tree.pep.faa
for sp in $(awk -v r="$ref" '$3 == r {print $2}' $design)
do
	cat $out/$sp.$ref.pep >> $out/$ref-tree.pep.faa
done
bash src/subgenome-tree.sh $out/$ref-combined.alignment.csv $cols \
	$out/$ref-tree.pep.faa $out/$ref-tree.pep.faa \
	$out/$ref-subgenome "$ind" $outgroup 0.5 $ncpu pep

echo
echo "ALL DONE"
