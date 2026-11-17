#!/bin/bash
# Subgenome tree: gene trees -> ASTRAL-hybrid species tree -> ML tree
#
# usage:
#   bash src/subgenome-tree.sh <alignment.csv> <rename.list> <pep.fa> <cds.fa> \
#        <outprefix> <indice> [root] [ncpu] [max_missing] [tree_from]
#
#   alignment.csv : matrix whose columns are subgenome copies (one row per gene family)
#   rename.list   : column-name table of that matrix (1-based index <TAB> name)
#   pep.fa/cds.fa : sequences of every gene id used by the selected columns
#   outprefix     : output prefix (e.g. out/<species>.<ref>-subgenome)
#   indice        : columns to keep, 1-based, e.g. 2,3,5,6 ; empty = all columns
#   root          : outgroup name used by ASTRAL --root and to reroot the ML tree
#   max_missing   : missing data allowed per gene family (default 0.5)
#   tree_from     : pep (protein alignment, default -- same as the source project) or cds
#
# outputs:
#   <outprefix>.genetrees      one Newick per gene family
#   <outprefix>.astral-hyb     ASTRAL-hybrid species tree
#   <outprefix>.alignments     concatenated alignment
#   <outprefix>.alignments.trimal.fasttree   quick ML tree (trimal -gappyout + fasttree -nt -gtr)
#   <outprefix>.alignments.trimal.iqtree.treefile   ML tree, GTR + 1000 ultrafast bootstrap (iqtree2)
#   <outprefix>.tmp/           per-family mafft/trimal alignments and gene trees
set -u

aln=$1
colnames=$2
pep=$3
cds=$4
pre=$5
indice=${6:-}
root=${7:-}
ncpu=${8:-8}
max_missing=${9:-0.5}
tree_from=${10:-pep}

[ -s "$aln" ]      || { echo "ERROR: no alignment file $aln"; exit 1; }
[ -s "$colnames" ] || { echo "ERROR: no colnames file $colnames"; exit 1; }
[ -s "$pep" ]      || { echo "ERROR: no pep file $pep"; exit 1; }
if [ "$tree_from" = "cds" ]; then
	[ -s "$cds" ] || { echo "ERROR: no cds file $cds (tree_from=cds needs codons)"; exit 1; }
else
	[ -s "$cds" ] || cds=$pep      # protein trees do not read the CDS; any file will do
fi

tmpd=${pre}.tmp
mkdir -p "$tmpd"
kargs="colnames=$colnames max_missing=$max_missing ncpu=$ncpu tree_from=$tree_from"
kargs="$kargs outdir=$tmpd genetrees=${pre}.genetrees alignments=${pre}.alignments"
[ -n "$indice" ] && kargs="$kargs indice=$indice"

##########################################################################
##### 1. gene trees (mafft -> trimal -> iqtree, inside soi-wgdi) ##########
##########################################################################
echo "== soi-wgdi genetrees $aln $pep $cds $kargs"
soi-wgdi genetrees "$aln" "$pep" "$cds" $kargs
[ -s "${pre}.genetrees" ] || { echo "ERROR: no ${pre}.genetrees"; exit 1; }
echo "   gene trees: $(grep -c . ${pre}.genetrees)"

##########################################################################
##### 2. species tree from gene trees (ASTRAL-hybrid) ####################
##########################################################################
opts=""
[ -n "$root" ] && opts="--root $root"
echo "== astral-hybrid $opts -i ${pre}.genetrees"
astral-hybrid $opts -i "${pre}.genetrees" -t 8 > "${pre}.astral-hyb" 2> "${pre}.astral-hyb.err" \
	|| { echo "ERROR: astral-hybrid failed, see ${pre}.astral-hyb.err"; exit 1; }
[ -s "${pre}.astral-hyb" ] || { echo "ERROR: astral-hybrid produced nothing"; exit 1; }
grep -o "^Score: .*" "${pre}.astral-hyb.err" | head -1

##########################################################################
##### 3. ML tree on the concatenation (trimal + fasttree) ################
##########################################################################
trimal -in "${pre}.alignments" -gappyout > "${pre}.alignments.trimal" \
	|| { echo "ERROR: trimal failed"; exit 1; }
# protein alignment -> fasttree default model (LG); nucleotide -> -nt -gtr
ft=""
[ "$tree_from" = "cds" ] && ft="-nt -gtr"
if [ -n "$root" ]; then
	fasttree $ft "${pre}.alignments.trimal" | nw_reroot - "$root" > "${pre}.alignments.trimal.fasttree" \
		|| { echo "ERROR: fasttree/nw_reroot failed"; exit 1; }
else
	fasttree $ft "${pre}.alignments.trimal" > "${pre}.alignments.trimal.fasttree" \
		|| { echo "ERROR: fasttree failed"; exit 1; }
fi

##########################################################################
##### 4. ML tree with IQ-TREE (GTR + ultrafast bootstrap 1000) ##########
##########################################################################
o=""
[ -n "$root" ] && o="-o $root"
mset=JTT
[ "$tree_from" = "cds" ] && mset=GTR
iqtree2 -redo -s "${pre}.alignments.trimal" -mset $mset -B 1000 -T $ncpu $o \
	-pre "${pre}.alignments.trimal.iqtree" \
	|| { echo "ERROR: iqtree2 failed"; exit 1; }
[ -s "${pre}.alignments.trimal.iqtree.treefile" ] || { echo "ERROR: no iqtree treefile"; exit 1; }

ls -la "${pre}.genetrees" "${pre}.astral-hyb" "${pre}.alignments" \
	"${pre}.alignments.trimal.fasttree" "${pre}.alignments.trimal.iqtree.treefile"
