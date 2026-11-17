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

A polyploid species carries several **subgenomes**: after a whole-genome duplication it
keeps two, three or five copies of most genes, one copy per subgenome. If another species
shares that duplication, its genes should also fall into the same subgenomes, i.e. each of
its copies should be sister to the corresponding subgenome of a reference species whose
subgenomes are already resolved.

This repository scores exactly that. Starting from the reads of a species, every gene is
assigned to one subgenome copy of a reference species, and gene trees are built to see
whether the assignment makes sense.

There are two routes; they differ in **how the copies are found**:

| route | idea | script |
|---|---|---|
| 1, reference-guided | reads compete for all copies of the reference; a read is assigned to at most one copy, and that copy is then assembled from those reads | `pipeline-refguided.sh` |
| 2, de novo | reads are assembled into contigs, annotated, and the resulting genes are assigned to the copies by orthology | `pipeline-denovo.sh` |

Route 1 keeps the depth information and is the one that tests the subgenomes. Route 2 works
from an assembly, so with target capture (HybSeq) reads the assembly only covers the captured
genes and the result is limited to them — it runs, but it is a target assembly, not a genome.
Both routes write the same matrix format, so samples of both can end up in one tree.

Every step below is a loop over `species.design`, and the two `pipeline-*.sh` scripts are
those loops in order, plus one combined tree per reference. Read the steps to understand what
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
4. **col4, library type** — `WGS` or `HybSeq`. Only printed with the log line; it tells you
   (and the log reader) how far the assembly of route 2 can reach: a HybSeq sample assembles
   only the captured genes.

The reference side is the directory `pre-study/`:

| file | content |
|---|---|
| `merged.alignment.csv` | matrix: one row per gene family, one **column per subgenome copy**; a cell is a gene id or empty |
| `rename.list` | column names, `<index> <TAB> <name>`, 1-based |
| `pep.faa` / `cds.fa` | proteins / CDS of every gene id that appears in the matrix |
| `all_species_gene.gff` | `contig  gene  start  end  strand` — only needed to split trees by chromosome |

Two conventions are used everywhere:

- gene ids are `<species>|<gene id>` (the `|` separator is what OrthoFinder and the tree
  tools expect);
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

`astral-hybrid` is a separate download from [ASTRAL-MP](https://github.com/smirarab/ASTRAL);
put the binary on `PATH`. Route 2 additionally needs

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
`src/retrieve_copies.py` then collects the per-gene sequences into `out/$sp.$ref.cds` and
`.pep`, with ids `$sp.$ref-<column>_XG<row>`; a paralog gets `-0`, `-1`… appended and is
dropped by the next step, since it does not belong to a single copy.

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
outgroup enter the tree. Gene trees are built from the CDS (codon-aware alignment, so indels
do not shift the frame) and combined by ASTRAL-hybrid; an ML tree on the concatenation is
built as a cross-check:

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
	bash src/subgenome-tree.sh out/$sp.$ref-appended.alignment.csv \
		out/$sp.$ref-appended.rename.list \
		out/$sp.$ref-tree.pep.faa out/$sp.$ref-tree.cds.fa \
		out/$sp.$ref-subgenome "$ind" $outgroup $ncpu
done < species.design
```

`src/subgenome-tree.sh` runs `soi-wgdi genetrees` (mafft → trimal → iqtree per gene family),
then `astral-hybrid --root $outgroup`, then `trimal -gappyout` + `fasttree -nt -gtr` +
`nw_reroot` for the ML tree.

What to look for: if the duplication is shared, the two copies of the sample are **not**
sisters — each sits next to the corresponding reference subgenome. If they are sisters, the
sample's copies do not resolve into the reference subgenomes.

`pipeline-refguided.sh` is these four loops plus one combined tree per reference, with every
sample of that reference in it (`out/<reference>-subgenome.*`).

### Route 2: de novo ###

Use this when you want the sample's own gene models, or when only an assembly exists. Every
design line is processed; for a `HybSeq` sample the assembly covers only the captured genes,
so its gene models — and everything downstream — are limited to those genes.

**Step 1 — assembly.**

```
while read -r run sp ref dtype; do
	megahit -1 reads/${run}_1.fastq.gz -2 reads/${run}_2.fastq.gz \
		-o out-denovo/$sp-assembly -t $ncpu --min-contig-len 300
done < species.design
```

**Step 2 — gene models.** MetaEuk predicts genes with the reference proteins as evidence, so
the genes found are the ones that can be compared with the reference copies. Protein and CDS
are then renumbered with one shared id series (`<species>|g000001`, `…|g000002`, the same
order in both files):

```
while read -r run sp ref dtype; do
	awk -v r="$ref" '/^>/{keep=($0 ~ ("^>" r "\\|"))} keep' pre-study/pep.faa > out-denovo/$ref.evidence.faa
	metaeuk easy-predict out-denovo/$sp-assembly/final.contigs.fa out-denovo/$ref.evidence.faa \
		out-denovo/$sp.metaeuk out-denovo/$sp.metaeuk.tmp --max-intron 500000 --threads $ncpu
	awk -v p="$sp" '/^>/{n++; printf ">%s|g%06d\n", p, n; next} {print}' \
		out-denovo/$sp.metaeuk.fas       > out-denovo/$sp.pep
	awk -v p="$sp" '/^>/{n++; printf ">%s|g%06d\n", p, n; next} {print}' \
		out-denovo/$sp.metaeuk.codon.fas > out-denovo/$sp.cds
done < species.design
```

**Step 3 — orthology.** OrthoFinder over the reference proteins, the outgroup and the sample
proteins — one fasta per species, named after the species:

```
mkdir -p out-denovo/orthofinder-in
while read -r run sp ref dtype; do
	awk -v r="$ref" '/^>/{keep=($0 ~ ("^>" r "\\|"))} keep' pre-study/pep.faa > out-denovo/orthofinder-in/$ref.fasta
	cp out-denovo/$sp.pep out-denovo/orthofinder-in/$sp.fasta
done < species.design
awk -v o="$outgroup" '/^>/{keep=($0 ~ ("^>" o "\\|"))} keep' pre-study/pep.faa > out-denovo/orthofinder-in/$outgroup.fasta
orthofinder -f out-denovo/orthofinder-in -t $ncpu
```

**Step 4 — assign genes to copies.** `soi-wgdi orthomap` turns the orthogroups into matrix
columns: a sample gene goes to at most one cell, a cell holds at most one gene, and the
reference genes keep their rows, so the matrix has the same shape as in route 1:

```
while read -r run sp ref dtype; do
	soi-wgdi orthomap pre-study/merged.alignment.csv out-denovo/orthofinder-in/OrthoFinder/Results_* \
		design=species.design colnames=pre-study/rename.list ref=$ref \
		out=out-denovo/$ref-orthomap.alignment.csv prefix=out-denovo/$ref-orthomap
done < species.design
```

It also writes `<prefix>.orthomap.tsv`, the audit table: one row per attempted assignment, with
`status` (assigned / gene already used / cell already taken) and `kind` (`unique` = the sample
gene is orthologous to **one** reference copy; `same`/`conflict` = to several, so the copy is
decided by competition, ties randomly). Watch the `unique` fraction: when most links are
`conflict`, the orthogroups do not resolve the copies and the copy labels of that sample mean
nothing — a limit of the data, not of the tree step.

**Step 5 — subgenome tree.** The same loop as route 1 step 4, reading
`out-denovo/$ref-orthomap.alignment.csv` and building `out-denovo/$ref-subgenome.*`;
`bash pipeline-denovo.sh species.design $ncpu` runs the whole route.

### Outputs ###

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
*Reevesia pubescens* (its two subgenomes `_R` and `_H`), and two samples whose full runs are
downloaded in step 0 — *Triplochiton scleroxylon* (`SRR34069007`, WGS, 5.9 GB) and *Ungeria
floribunda* (`ERR7599640`, target capture, 0.3 GB):

```
bash src/get_reads.sh species.design reads 4
bash pipeline-refguided.sh species.design 8
bash pipeline-denovo.sh species.design 8
```

Measured on the downloaded runs (route 1):

| | *Ungeria floribunda* | *Triplochiton scleroxylon* |
|---|---|---|
| subfamily of the sample | Helicteroideae (same as the reference) | Tilioideae |
| library | HybSeq (target capture) | WGS |
| mean depth per subgenome copy | 4.5x / 4.7x | — |
| copies assembled (of 114 baits) | 42 | — |
| gene trees in the subgenome tree | 27 | — |
| topology | `(Ref_R, Ung_R)` and `(Ref_H, Ung_H)`, pp = 1.00 / 0.9999 | — |

The sample from the same subfamily resolves cleanly into the two reference subgenomes, which
is the signal this pipeline is built to detect.

Route 2 on the same HybSeq sample, for comparison: 26,921 contigs, 20 gene models, 42
ortholog pairs of which only 6 are `unique` (36 `conflict`, copy decided randomly), and a tree
in which the two sample copies stay sisters next to one reference subgenome — i.e. the copies
are not resolved. With 58 families and a 2-copy reference most orthogroups contain both
copies, so route 2 is the fallback, not the primary evidence; the `unique` fraction of its
audit table tells you how much of the answer to trust.

### FAQ ###

1. **The design file has no header and no comments.** A comment line would be read as a
   sample, and rows are read by position. Four TAB separated columns, as above.
2. **What are `_R` and `_H`?** The two subgenomes of the reference. A polyploid lineage keeps
   one copy of every gene per subgenome, so each copy is one matrix column named
   `<species>_<subgenome>`; the labels are free (`_A/_B/_C` works the same way), only the
   `<species>_<subgenome>` shape matters.
3. **Which column is used for what?** `col1` locates the reads (`reads/<col1>_{1,2}.fastq.gz`),
   `col2` labels the output, `col3` selects the reference, `col4` only describes the library
   (it decides nothing; it is read into the log so that a shallow or capture-limited result
   can be recognised).
4. **Few genes assembled (route 1).** Look at `<species>.<reference>.depth.tsv` and
   `.depth.pdf` first: if the copies are covered very unevenly, the sample either does not
   share the duplication or its reads are too shallow for that reference. A target capture
   experiment covers only the captured genes.
5. **Almost no `unique` assignments (route 2).** Several reference copies per orthogroup, so
   the copy is decided by competition and the labels are weak. Add gene families, use a closer
   reference, or use route 1.
6. **Re-running.** Nothing is checkpointed: every step overwrites its outputs. Delete `out/`
   (or `out-denovo/`) for a clean run.
7. **Trees per chromosome.** `soi-wgdi chrsplit` splits the gene trees of
   `<prefix>-subgenome.tmp` by the chromosome of the reference genes in
   `all_species_gene.gff`; each part can be passed to `astral-hybrid` separately.
