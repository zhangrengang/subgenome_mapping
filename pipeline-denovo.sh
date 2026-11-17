#!/bin/bash
# Subgenome mapping -- de novo route.
#
# reads --> de novo assembly (megahit) --> gene models (MetaEuk)
#       --> orthology against the reference (OrthoFinder)
#       --> every sample gene assigned to a copy of the reference
#       --> subgenome gene trees
#
# usage: bash pipeline-denovo.sh [species.design] [ncpu] [outgroup]
#
# species.design: same file as the reference-guided route, 4 columns, no header:
#   <run>  <species_id>  <reference_species>  <WGS|HybSeq>
#     run        accession; reads/<run>_clean_{1,2}.fastq.gz (src/get_reads.sh,
#                then src/trim_reads.sh)
#     species_id label of this sample; its gene ids become <species_id>|g000001
#     reference  species whose copies are columns of pre-study/merged.alignment.csv
#     dtype      library type; processed the same way, printed with the log line.
#                HybSeq (target capture) reads only cover the captured genes, so the
#                assembly of such a sample is a target assembly, not a genome
#
# requires on PATH: megahit, metaeuk, orthofinder, soi-wgdi
set -eu

design=${1:-species.design}
ncpu=${2:-8}
outgroup=${3:-Theobroma_cacao}
pre=pre-study
reads=reads
out=out-denovo
mkdir -p $out

##########################################################################
##### 1. de novo assembly ###############################################
##########################################################################
while read -r sample sp ref dtype
do
	[ -n "${sample:-}" ] || continue
	[ "$dtype" = "WGS" ] || echo "== note: $sp is $dtype, its assembly covers only the captured genes"
	echo "########## $sp: assembly ##########"
	##################################################################
	##### megahit ####################################################
	##################################################################
	# per-sample directory: assembly and gene models live together
	[ -e $out/$sp ] && rm -rf $out/$sp
	megahit -1 $reads/${sample}_clean_1.fastq.gz -2 $reads/${sample}_clean_2.fastq.gz \
		-o $out/$sp -t $ncpu
	echo "   contigs: $(grep -c '^>' $out/$sp/final.contigs.fa)"
done < $design

##########################################################################
##### 2. gene models: MetaEuk, ids renumbered per sample ################
##########################################################################
while read -r sample sp ref dtype
do
	[ -n "${sample:-}" ] || continue
	echo "########## $sp: gene prediction ##########"
	##################################################################
	##### metaeuk easy-predict #######################################
	# the reference proteins are the evidence, so that the predicted genes
	# are the ones that can be compared with the reference copies
	##################################################################
	awk -v r="$ref" '/^>/{keep=($0 ~ ("^>" r "\\|"))} keep' $pre/pep.faa > $out/$ref.evidence.faa
	metaeuk easy-predict $out/$sp/final.contigs.fa $out/$ref.evidence.faa \
		$out/$sp/$sp.metaeuk $out/$sp/tmp.metaeuk --max-intron 500000 --threads $ncpu
	##################################################################
	##### renumber: protein and CDS share one id series ##############
	##################################################################
	mkdir -p $out/orthofinder-in
	awk -v p="$sp" '/^>/{n++; printf ">%s|g%06d\n", p, n; next} {print}' $out/$sp/$sp.metaeuk.fas      > $out/$sp/$sp.pep
	awk -v p="$sp" '/^>/{n++; printf ">%s|g%06d\n", p, n; next} {print}' $out/$sp/$sp.metaeuk.codon.fas > $out/$sp/$sp.cds
	ln -sf ../$sp/$sp.pep $out/orthofinder-in/$sp.fasta
	echo "   genes: $(grep -c '^>' $out/$sp/$sp.pep)"
done < $design

##########################################################################
##### 3. orthology: reference copies vs sample genes ####################
##########################################################################
echo "########## OrthoFinder ##########"
# one matrix = one reference species; a design mixing references has to be split
nref=$(awk '{print $3}' $design | sort -u | wc -l)
[ "$nref" -eq 1 ] || { echo "ERROR: $nref references in $design but the matrix is one reference; split the design"; exit 1; }
ref=$(awk '{print $3}' $design | sort -u)
mkdir -p $out/orthofinder-in
awk -v r="$ref" '/^>/{keep=($0 ~ ("^>" r "\\|"))} keep' $pre/pep.faa > $out/orthofinder-in/$ref.fasta
awk -v o="$outgroup" '/^>/{keep=($0 ~ ("^>" o "\\|"))} keep' $pre/pep.faa > $out/orthofinder-in/$outgroup.fasta
orthofinder -f $out/orthofinder-in -t $ncpu

##########################################################################
##### 4. assign sample genes to the copies ##############################
##########################################################################
echo "########## orthomap ##########"
##################################################################
##### soi-wgdi orthomap ##########################################
# every sample gene gets at most one cell; one cell holds at most one gene;
# the reference genes keep their rows. ref is inferred from design + results
# + matrix column names (must be unique)
##################################################################
soi-wgdi orthomap $pre/merged.alignment.csv $out/orthofinder-in/OrthoFinder/Results_* \
	design=$design colnames=$pre/rename.list \
	out=$out/$ref-orthomap.alignment.csv prefix=$out/$ref-orthomap

##########################################################################
##### 5. subgenome tree #################################################
##########################################################################
echo "########## subgenome tree: $ref ##########"
aln=$out/$ref-orthomap.alignment.csv
cols=$out/$ref-orthomap.rename.list
# every column except the first (the anchor P_MK) enters the tree:
# reference copies + outgroup + all sample columns
ind=$(seq -s, 2 $(wc -l < $cols))
echo "   columns: $ind"
awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' $pre/pep.faa > $out/$ref-tree.pep.faa
cat $(ls $out/orthofinder-in/*.fasta | grep -vE "/($ref|$outgroup)\.fasta$") >> $out/$ref-tree.pep.faa
bash src/subgenome-tree.sh $aln $cols $out/$ref-tree.pep.faa $out/$ref-tree.pep.faa \
	$out/$ref-subgenome "$ind" $outgroup $ncpu pep

echo
echo "ALL DONE"
