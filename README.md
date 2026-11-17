## Table of Contents

   * [What it does](#what-it-does)
   * [Input data](#input-data)
   * [Installation](#installation)
   * [Step 0: get the reads](#step-0-get-the-reads)
   * [Route 1: reference-guided](#route-1-reference-guided)
   * [Route 2: de novo](#route-2-de-novo)
   * [Outputs](#outputs)
   * [The example data](#the-example-data)
   * [FAQ](#faq)

### What it does ###

A polyploid lineage keeps several **subgenomes**: after a whole-genome duplication it carries
two, three or five copies of most genes, one copy per subgenome. For one reference species
those copies have already been assigned, and that assignment is the starting point here: the
file `pre-study/merged.alignment.csv`, described under *Input data*.

The question this repository answers is narrower than "are there several subgenomes": **does
a species known only from sequencing reads share that duplication?** If it does, its genes
should fall into the same subgenomes, i.e. each of its copies should be sister to the
corresponding copy of the reference instead of the sample's copies being sisters to each
other.

So the work is (i) to decide, for every reference gene family, which copy a read-derived
sequence of the sample belongs to, and (ii) to build trees that show whether that assignment
recovers the reference's subgenomes. Two routes do (i) differently:

| route | idea | script |
|---|---|---|
| 1, reference-guided | reads compete for all copies of the reference; a read is assigned to at most one copy, and that copy is then assembled from those reads | `pipeline-refguided.sh` |
| 2, de novo | the sample is assembled and annotated on its own, and its genes are assigned to the copies by orthology | `pipeline-denovo.sh` |

Route 1 keeps the depth information and is the route that tests the duplication. Route 2 does
not use the reference for assembly at all, only afterwards for the assignment; with target
capture (HybSeq) reads its assembly covers only the captured genes, so the result is limited
to them. Both routes write the same matrix format, so samples of both can end up in one tree.

Every step below is a loop over `species.design`, and the two `pipeline-*.sh` scripts are
those loops in order plus one combined tree per reference. Read the steps to understand what
happens; run the scripts when you have many samples.

### Input data ###

One file drives everything, `species.design`: four columns, TAB separated, **no header, no
comments**, one line per sample:

```
SRR34069007	Triplochiton_scleroxylon	Reevesia_pubescens	WGS
ERR7599640	Ungeria_floribunda	Reevesia_pubescens	HybSeq
```

1. **col1, run accession** — locates the reads at ENA/NCBI. If you already have clean reads,
   put them as `reads/<col1>_1.fastq.gz` and `reads/<col1>_2.fastq.gz` and skip Step 0.
2. **col2, species id** — the label of this sample in the output. Matrix columns become
   `<species id>_<subgenome>`; de novo gene ids become `<species id>|g000001`. A species
   sequenced twice needs two ids (`_WGS`, `_1`, `_2`).
3. **col3, reference species** — whose subgenome copies are columns of the matrix. These
   copies are the baits (route 1) and the assignment target (route 2).
4. **col4, library type** — `WGS` or `HybSeq`; only printed with the log line, it tells the
   reader how far the assembly of route 2 can reach.

The reference side is the directory `pre-study/`:

| file | content |
|---|---|
| `merged.alignment.csv` | subgenome-assigned matrix, the output format of `wgdi -at`: one row per gene family, one **column per subgenome copy**, a cell = gene id or empty |
| `rename.list` | column names of that matrix, line *n* = column *n*, name = the second field (the first is just a running number) |
| `pep.faa` / `cds.fa` | proteins / CDS of every gene id that appears in the matrix |
| `all_species_gene.gff` | `contig  gene  start  end  strand` — only needed to split trees by chromosome |

Two conventions are used everywhere:

- gene ids are `<species>|<gene id>` (the `|` separator is what OrthoFinder and the tree tools
  expect);
- a copy column is named `<species>_<subgenome>`, e.g. `Reevesia_pubescens_R` and
  `Reevesia_pubescens_H`; a species without copies (the outgroup, `Theobroma_cacao`) is a
  single column. The labels are free — `_A/_B/_C` works the same way.

### Installation ###

[SOI / OrthoIndex](https://github.com/zhangrengang/SOI#installation) provides `soi-wgdi` and
the environment (mafft, trimal, iqtree, fasttree, newick_utils, orthofinder, blast, diamond…):

```
git clone https://github.com/zhangrengang/SOI.git
cd SOI
mamba env create -f OrthoIndex.yaml
mamba activate OrthoIndex
pip3 install .
```

Not included there, install into the same environment (route 1 and the trees):

```
mamba install -n OrthoIndex -c bioconda bwa samtools spades exonerate
```

`astral-hybrid` comes from [ASTER](https://github.com/chaoszhang/ASTER) ("Accurate Species
Tree EstimatoR"), which is a conda package too:

```
mamba install -n OrthoIndex -c bioconda aster
```

Route 2 additionally needs

```
mamba install -n OrthoIndex -c bioconda megahit metaeuk
```

Route 1 also needs [HybPiper](https://github.com/mossmatters/HybPiper) (it is what calls
spades and exonerate, so install it as its own README says) and `src/retrieve_copies.py`
from this repository (`biopython` is already in the environment). The pipeline calls
HybPiper's scripts through an interpreter, so nothing has to be made executable:

```
export HYBPIPER=/path/to/HybPiper/bin   # the directory holding reads_first.py
export PYTHON=python                    # interpreter that install needs (python2 for HybPiper 1.2)
```

Step 0 uses only `curl`.

### Step 0: get the reads ###

```
bash src/get_reads.sh species.design reads 4      # design, output directory, parallel downloads
```

For each line of the design the script asks ENA (the SRA mirror) for the fastq links of the
run in col1, through the portal API

```
https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR34069007&result=read_run&fields=fastq_ftp&format=tsv
```

and downloads them as `reads/<run>_1.fastq.gz` and `reads/<run>_2.fastq.gz`; single-end runs
give one file. With the example design that is `reads/SRR34069007_{1,2}.fastq.gz` (5.9 GB,
WGS) and `reads/ERR7599640_{1,2}.fastq.gz` (0.3 GB, target capture). Runs already downloaded
are skipped, so you can add lines to the design and re-run. Reads are used as downloaded; if
your own data needs adapter/quality trimming, trim first and keep the same file names.

### Route 1: reference-guided ###

The point of this route is that **no read is allowed to support two copies**. Everything else
follows from that. Two variables are used by all steps below:

```
outgroup=Theobroma_cacao
ncpu=8
```

All loops read four fields. It has to be four: `read` puts the *rest* of the line into the
last variable, so `while read -r run sp ref` would hand `Reevesia_pubescens<TAB>WGS` to `$ref`.
If the fourth column is of no interest, use `while read -r run sp ref _`.

**Step 1 — competitive mapping.** All CDS of the reference species are used as baits, not only
the genes that are complete in the reference: a read coming from one copy has to be able to
compete with the other copies, otherwise a conserved gene would pile up on whichever copy
happens to be in the bait set and the depth comparison would be meaningless. After mapping,
only primary alignments are kept (a read pair lands on one reference gene at most), and only
the genes whose copies are all present are written into the bam:

```
while read -r run sp ref dtype; do
	soi-wgdi compmap pre-study/merged.alignment.csv pre-study/cds.fa \
		reads/${run}_1.fastq.gz reads/${run}_2.fastq.gz \
		colnames=pre-study/rename.list ref=$ref \
		prefix=out/$sp.$ref/$sp.$ref ncpu=$ncpu
done < species.design
```

The prefix decides where the bam of step 2 has to be, so keep the shape
`out/<species>.<reference>/<species>.<reference>`. Outputs in `out/$sp.$ref/`:

- `$sp.$ref.depth.tsv` — mean depth of every reference gene in every subgenome copy;
- `$sp.$ref.depth.pdf` — the same per copy. Read this first: if the copies are covered very
  unevenly, the sample either does not share the duplication or its reads are too shallow for
  this reference. An even profile proves nothing by itself;
- `$sp.$ref.target.cds` — baits of the retained genes, i.e. what step 2 assembles;
- `$sp.$ref.bam` — primary-only alignments of those baits.

**Step 2 — assemble each copy.** HybPiper runs with `--bwa --no-blast`, i.e. it reuses the bam
of step 1 instead of mapping again. Each bait becomes one gene directory with one SPAdes +
exonerate run, so every copy is assembled separately, from the reads that were assigned to it:

```
while read -r run sp ref dtype; do
	( cd out
	  touch dummy.fastq
	  $PYTHON $HYBPIPER/reads_first.py -r dummy.fastq -b $sp.$ref/$sp.$ref.target.cds \
		--prefix $sp.$ref --cpu $ncpu --bwa --no-blast
	  $PYTHON $HYBPIPER/paralog_investigator.py $sp.$ref )
	python3 src/retrieve_copies.py out/$sp.$ref
done < species.design
```

HybPiper must be started from the directory that *contains* `<prefix>/` and gets the bare
directory name as `--prefix` — that is why the loop starts with `cd out`.
`src/retrieve_copies.py` collects the per-gene sequences into `out/$sp.$ref.cds` and `.pep`,
with ids `$sp.$ref-<column>_XG<row>` (the id is what `append_cols` reads back: the part after
the last `-` is the matrix cell, so the column and row are recovered from it). A paralog gets
`-0`, `-1`… appended and is dropped by the next step, since it does not belong to one copy.

**Step 3 — append the sample as new columns.** One column per copy of the reference; a cell is
filled only where the sample has a sequence for that gene family and that copy:

```
while read -r run sp ref dtype; do
	soi-wgdi append_cols pre-study/merged.alignment.csv out/$sp.$ref.cds $sp \
		colnames=pre-study/rename.list prefix=out/$sp.$ref-appended
done < species.design
```

Outputs `out/$sp.$ref-appended.alignment.csv` and `.rename.list`: the reference matrix plus
`<species>_<subgenome>` columns, named after the copies of the reference.

**Step 4 — subgenome tree.** Only the columns of this reference, of this sample and of the
outgroup enter the tree. Three things happen: gene trees per family from the codon-aware CDS
alignment, ASTRAL-hybrid on those gene trees, and an ML tree on the concatenation as a
cross-check. All three are written out below, no wrapper:

```
while read -r run sp ref dtype; do
	ind=$(awk -F'\t' -v r="$ref" -v s="$sp" -v o="$outgroup" \
		'{n=$2; if (n ~ ("^" r "_") || n ~ ("^" s "_") || n == o) printf "%s%s", (k++ ? "," : ""), $1}' \
		out/$sp.$ref-appended.rename.list)

	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' \
		pre-study/pep.faa > out/$sp.$ref-tree.pep.faa
	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' \
		pre-study/cds.fa > out/$sp.$ref-tree.cds.fa
	cat out/$sp.$ref.pep >> out/$sp.$ref-tree.pep.faa
	cat out/$sp.$ref.cds >> out/$sp.$ref-tree.cds.fa

	# 4a. one gene tree per family (mafft -> trimal -> iqtree), plus the concatenation
	soi-wgdi genetrees out/$sp.$ref-appended.alignment.csv out/$sp.$ref-tree.pep.faa \
		out/$sp.$ref-tree.cds.fa colnames=out/$sp.$ref-appended.rename.list \
		indice=$ind max_missing=0.5 ncpu=$ncpu tree_from=cds \
		outdir=out/$sp.$ref-subgenome.tmp \
		genetrees=out/$sp.$ref-subgenome.genetrees alignments=out/$sp.$ref-subgenome.alignments

	# 4b. species tree from the gene trees
	astral-hybrid --root $outgroup -t 8 -i out/$sp.$ref-subgenome.genetrees \
		> out/$sp.$ref-subgenome.astral-hyb

	# 4c. ML tree on the concatenation (gappyout trimming, then fasttree, then reroot)
	trimal -in out/$sp.$ref-subgenome.alignments -gappyout \
		> out/$sp.$ref-subgenome.alignments.trimal
	fasttree -nt -gtr out/$sp.$ref-subgenome.alignments.trimal \
		| nw_reroot - $outgroup > out/$sp.$ref-subgenome.alignments.trimal.fasttree
done < species.design
```

What to look for: if the duplication is shared, the two copies of the sample are **not**
sisters — each sits next to the corresponding reference subgenome. If they are sisters, the
sample's copies do not resolve into the reference subgenomes.

`pipeline-refguided.sh` is these four loops plus one combined tree per reference, with every
sample of that reference in it (`out/<reference>-subgenome.*`); it calls
`src/subgenome-tree.sh`, which is exactly the 4a–4c block above.

### Route 2: de novo ###

Here the sample is assembled and annotated on its own, and only then are its genes compared
with the reference — instead of assembling it on the reference baits as in route 1. Use it when
you want the sample's own gene models (its own introns, its own UTRs) or when the sample is too
divergent for the baits.

The three steps are separable: if you already have an assembly, put it at
`out-denovo/<species>-assembly/final.contigs.fa` and start at step 2 (annotation); if you
already have gene models, write them as `out-denovo/orthofinder-in/<species>.fasta`
(proteins, ids `<species>|g…`) plus `out-denovo/<species>.cds` and start at step 3
(OrthoFinder). A `HybSeq` sample is processed like any other, but its reads cover
only the captured genes, so its assembly — and everything downstream — is limited to them.

**Step 1 — assembly.**

```
while read -r run sp ref dtype; do
	megahit -1 reads/${run}_1.fastq.gz -2 reads/${run}_2.fastq.gz \
		-o out-denovo/$sp-assembly -t $ncpu
done < species.design
```

**Step 2 — gene models.** MetaEuk predicts genes with the reference proteins as evidence, so
the genes found are the ones that can be compared with the reference copies. Protein and CDS
are then renumbered with one shared id series (`<species>|g000001`, `…|g000002`, the same order
in both files); the protein file is written straight into the OrthoFinder input directory
(one file per species, named after the species):

```
while read -r run sp ref dtype; do
	awk -v r="$ref" '/^>/{keep=($0 ~ ("^>" r "\\|"))} keep' pre-study/pep.faa > out-denovo/$ref.evidence.faa
	metaeuk easy-predict out-denovo/$sp-assembly/final.contigs.fa out-denovo/$ref.evidence.faa \
		out-denovo/$sp.metaeuk out-denovo/$sp.metaeuk.tmp --max-intron 500000 --threads $ncpu
	mkdir -p out-denovo/orthofinder-in
	awk -v p="$sp" '/^>/{n++; printf ">%s|g%06d\n", p, n; next} {print}' \
		out-denovo/$sp.metaeuk.fas       > out-denovo/orthofinder-in/$sp.fasta
	awk -v p="$sp" '/^>/{n++; printf ">%s|g%06d\n", p, n; next} {print}' \
		out-denovo/$sp.metaeuk.codon.fas > out-denovo/$sp.cds
done < species.design
```

**Step 3 — orthology.** One OrthoFinder run over the whole input directory — the reference
proteins and the outgroup are added once, the sample files are already there from step 2:

```
for ref in $(awk '{print $3}' species.design | sort -u); do
	awk -v r="$ref" '/^>/{keep=($0 ~ ("^>" r "\\|"))} keep' pre-study/pep.faa > out-denovo/orthofinder-in/$ref.fasta
done
awk -v o="$outgroup" '/^>/{keep=($0 ~ ("^>" o "\\|"))} keep' pre-study/pep.faa > out-denovo/orthofinder-in/$outgroup.fasta
orthofinder -f out-denovo/orthofinder-in -t $ncpu
```

**Step 4 — assign genes to copies.** `soi-wgdi orthomap` turns the orthogroups into matrix
columns: a sample gene goes to at most one cell, a cell holds at most one gene, and the
reference genes keep their rows, so the matrix has the same shape as in route 1:

```
for ref in $(awk '{print $3}' species.design | sort -u); do
	soi-wgdi orthomap pre-study/merged.alignment.csv out-denovo/orthofinder-in/OrthoFinder/Results_* \
		design=species.design colnames=pre-study/rename.list ref=$ref \
		out=out-denovo/$ref-orthomap.alignment.csv prefix=out-denovo/$ref-orthomap
done
```

It also writes `<prefix>.orthomap.tsv`, the audit table: one row per attempted assignment, with
`status` (assigned / gene already used / cell already taken) and `kind` (`unique` = the sample
gene is orthologous to **one** reference copy; `same`/`conflict` = to several, so the copy is
decided by competition, ties randomly). Watch the `unique` fraction: when most links are
`conflict`, the orthogroups do not resolve the copies and the copy labels of that sample mean
nothing — a limit of the data, not of the tree step.

**Step 5 — subgenome tree.** The same three commands as route 1 step 4, on the orthomap
matrix; here the sample genes come from the de novo gene models:

```
for ref in $(awk '{print $3}' species.design | sort -u); do
	ind=$(awk -F'\t' -v r="$ref" -v o="$outgroup" \
		'{n=$2; if (n ~ ("^" r "_") || n == o) printf "%s%s", (k++ ? "," : ""), $1}' \
		out-denovo/$ref-orthomap.rename.list)
	for sp in $(awk -v r="$ref" '$3 == r {print $2}' species.design); do
		ind="$ind$(awk -F'\t' -v s="$sp" '$2 ~ ("^" s "_") {printf ",%s", $1}' \
			out-denovo/$ref-orthomap.rename.list)"
	done

	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' \
		pre-study/pep.faa > out-denovo/$ref-tree.pep.faa
	awk -v r="$ref" -v o="$outgroup" '/^>/{keep=($0 ~ ("^>(" r "|" o ")\\|"))} keep' \
		pre-study/cds.fa > out-denovo/$ref-tree.cds.fa
	for sp in $(awk -v r="$ref" '$3 == r {print $2}' species.design); do
		cat out-denovo/orthofinder-in/$sp.fasta >> out-denovo/$ref-tree.pep.faa
		cat out-denovo/$sp.cds >> out-denovo/$ref-tree.cds.fa
	done

	soi-wgdi genetrees out-denovo/$ref-orthomap.alignment.csv out-denovo/$ref-tree.pep.faa \
		out-denovo/$ref-tree.cds.fa colnames=out-denovo/$ref-orthomap.rename.list \
		indice=$ind max_missing=0.5 ncpu=$ncpu tree_from=cds \
		outdir=out-denovo/$ref-subgenome.tmp \
		genetrees=out-denovo/$ref-subgenome.genetrees alignments=out-denovo/$ref-subgenome.alignments
	astral-hybrid --root $outgroup -t 8 -i out-denovo/$ref-subgenome.genetrees \
		> out-denovo/$ref-subgenome.astral-hyb
	trimal -in out-denovo/$ref-subgenome.alignments -gappyout \
		> out-denovo/$ref-subgenome.alignments.trimal
	fasttree -nt -gtr out-denovo/$ref-subgenome.alignments.trimal \
		| nw_reroot - $outgroup > out-denovo/$ref-subgenome.alignments.trimal.fasttree
done
```

### Outputs ###

The two routes produce files with the same names (`<reference>-subgenome.*`), so they write to
different directories: route 1 to `out/`, route 2 to `out-denovo/`. Mixing them in one
directory would silently overwrite trees.

Route 1, in `out/`:

```
<species>.<reference>.depth.tsv/.depth.pdf       per-copy coverage of the reference genes
<species>.<reference>.cds/.pep                   assembled copies of this sample
<species>.<reference>-appended.*                 matrix with this sample's columns
<species>.<reference>-subgenome.genetrees        one tree per gene family
<species>.<reference>-subgenome.astral-hyb       ASTRAL-hybrid subgenome tree
<species>.<reference>-subgenome.alignments*      concatenation (+ .trimal.fasttree = ML tree)
<reference>-subgenome.*                          the same tree for all samples of a reference
```

Route 2, in `out-denovo/`:

```
<species>-assembly/final.contigs.fa              assembly
<species>.pep/.cds/.gff                          gene models, ids <species>|g000001…
orthofinder-in/OrthoFinder/Results_*/            orthology
<reference>-orthomap.alignment.csv/.rename.list  matrix with the new columns
<reference>-orthomap.orthomap.tsv                assignment audit table
<reference>-subgenome.*                          subgenome tree
```

Trees are plain Newick. For figures use your own viewer (e.g.
[PhyTop](https://github.com/zhangrengang/phytop) or FigTree).

### The example data ###

`species.design` and `pre-study/` are a small real example: 58 gene families of the reference
*Reevesia pubescens* (its two subgenomes `_R` and `_H`), and two samples — *Triplochiton
scleroxylon* (`SRR34069007`, WGS, 5.9 GB) and *Ungeria floribunda* (`ERR7599640`, target
capture, 0.3 GB) — whose full runs are downloaded in step 0:

```
bash src/get_reads.sh species.design reads 4
bash pipeline-refguided.sh species.design 8
bash pipeline-denovo.sh species.design 8
```

Measured with `bash pipeline-refguided.sh` on the downloaded runs:

| | *Ungeria floribunda* |
|---|---|
| subfamily of the sample | Helicteroideae (same as the reference) |
| library | HybSeq (target capture) |
| mean depth per subgenome copy | 4.5x / 4.7x |
| copies assembled (of 114 baits) | 42 |
| gene trees in the subgenome tree | 27 |
| topology | `(Ref_R, Ung_R)` and `(Ref_H, Ung_H)`, pp = 1.00 / 0.9999 |

The sample from the same subfamily resolves cleanly into the two reference subgenomes, which is
the signal this pipeline is built to detect. (*Triplochiton*, from another subfamily, is in the
design and runs the same way; its numbers are not listed here yet.)

Route 2 on the same HybSeq sample, for comparison: 26,921 contigs, 20 gene models, 42 ortholog
pairs of which only 6 are `unique` (36 `conflict`, copy decided randomly), and a tree in which
the two sample copies stay sisters next to one reference subgenome — i.e. the copies are not
resolved. With 58 families and a 2-copy reference most orthogroups contain both copies, so
route 2 is the fallback, not the primary evidence; the `unique` fraction of its audit table
tells you how much of the answer to trust.

### FAQ ###

1. **The design file has no header and no comments.** A comment line would be read as a
   sample, and rows are read by position. Four TAB separated columns, as above.
2. **What are `_R` and `_H`?** The two subgenomes of the reference. A polyploid lineage keeps
   one copy of every gene per subgenome, so each copy is one matrix column named
   `<species>_<subgenome>`; the labels are free (`_A/_B/_C` works the same way), only the
   `<species>_<subgenome>` shape matters.
3. **Few genes assembled (route 1).** Look at `<species>.<reference>.depth.tsv` and
   `.depth.pdf` first: if the copies are covered very unevenly, the sample either does not
   share the duplication or its reads are too shallow for that reference. A target capture
   experiment covers only the captured genes.
4. **Almost no `unique` assignments (route 2).** Several reference copies per orthogroup, so
   the copy is decided by competition and the labels are weak. Add gene families, use a closer
   reference, or use route 1.
5. **Re-running.** Nothing is checkpointed: every step overwrites its outputs. Delete `out/`
   (or `out-denovo/`) for a clean run.
6. **Trees per chromosome.** `soi-wgdi chrsplit` splits the gene trees of
   `<prefix>-subgenome.tmp` by the chromosome of the reference genes in
   `all_species_gene.gff`; each part can be passed to `astral-hybrid` separately.
