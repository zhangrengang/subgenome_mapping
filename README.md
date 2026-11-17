## Table of Contents

   * [What it does](#what-it-does)
   * [Prepare data](#prepare-data)
   * [Installation](#installation)
   * [Route 1: reference-guided (reads)](#route-1-reference-guided-reads)
   * [Route 2: de novo (reads -> assembly -> gene models)](#route-2-de-novo-reads---assembly---gene-models)
   * [Outputs](#outputs)
   * [The example data](#the-example-data)
   * [FAQ](#faq)

### What it does ###

A polyploid species carries several **subgenomes** (its WGD copies). If another species
shares that WGD, its genes should also fall into the same subgenome copies. This pipeline
scores that: starting from reads of a species, it assigns every gene to one copy of a
reference species whose subgenomes are already resolved, then builds gene trees to see
whether each subgenome of the sample is sister to the corresponding subgenome of the
reference.

Two routes, same downstream:

| route | input | how the copies are found | pipeline |
|---|---|---|---|
| 1 reference-guided | reads | reads compete for all copies of the reference, each read being assigned to at most one; the retained copy is then assembled from those reads | `pipeline-refguided.sh` |
| 2 de novo | reads | reads are assembled and annotated first; sample genes are then assigned to the reference copies by orthology | `pipeline-denovo.sh` |

Both write the same matrix format, so samples of both routes can end up in one tree.

### Prepare data ###

Everything is driven by one design file, `design.txt`: three columns, **no header, no
comments**, one line per sample:

```
Triplochiton_scleroxylon	Triplochiton_scleroxylon	Reevesia_pubescens
Ungeria_floribunda	Ungeria_floribunda	Reevesia_pubescens
```

1. **sample id** (`col1`) – locates the reads: `reads/<sample id>_filter_R1.fastq.gz`
   and `reads/<sample id>_filter_R2.fastq.gz` (gzipped, cleaned/trimmed).
2. **species id** (`col2`) – the label of this sample in the output. New matrix columns
   are named `<species id>_<copy>`; de novo gene ids are `<species id>|g000001`.
   Give two datasets of the same species different ids (`_WGS`, `_1`, `_2`).
3. **reference species** (`col3`) – one of the species whose copies are columns of the
   matrix. Its copies are the baits (route 1) or the assignment target (route 2).

The reference side lives in `pre-study/`:

```
pre-study/merged.alignment.csv   matrix, one row per gene family, one COLUMN per
                                 subgenome copy; a cell = a gene id (or empty)
pre-study/rename.list            column names: "<index> <TAB> <name>" (1-based)
pre-study/pep.faa                proteins of every gene id in the matrix
pre-study/cds.fa                 CDS of the same genes
pre-study/all_species_gene.gff   contig/chromosome coordinates (needed only to
                                 split trees by chromosome)
```

Column names follow `<species>_<copy>`, e.g. `Reevesia_pubescens_R` /
`Reevesia_pubescens_H`; a species without copies (the outgroup, e.g.
`Theobroma_cacao`) is a single column. Gene ids are `<species>|<gene id>`.

### Installation ###

[OrthoIndex](https://github.com/zhangrengang/orthoindex#installation) provides `soi-wgdi`
and the whole conda environment:

```
git clone https://github.com/zhangrengang/orthoindex.git
cd orthoindex
conda env create -f OrthoIndex.yaml
conda activate OrthoIndex
python3 setup.py install
```

Route 1 additionally needs [HybPiper](https://github.com/mossmatters/HybPiper) and
`src/retrieve_copies.py` (this repository, needs `biopython` which OrthoIndex already
has). The pipeline calls HybPiper's scripts through an interpreter, so nothing has to be
made executable:

```
export HYBPIPER=/path/to/HybPiper/bin   # directory holding reads_first.py
export PYTHON=python                    # the interpreter that install needs (python2 for HybPiper 1.2)
```

`bwa` and `samtools` are already in OrthoIndex.

Route 2 additionally needs `megahit`, `metaeuk` and `OrthoFinder`.

The tree step uses `astral-hybrid` ([ASTRAL-MP](https://github.com/smirarab/ASTRAL)),
`trimal`, `fasttree` and `newick_utils` (`nw_reroot`); the last three are in OrthoIndex.

### Route 1: reference-guided (reads) ###

Run the whole route for every line of the design:

```
bash pipeline-refguided.sh design.txt 8          # 8 = threads per sample
```

The script does four things per sample; each step can also be run by hand.

**Step 1 – competitive mapping.** The reads are mapped against **all** CDS of the
reference species, not only the genes that are complete in it: a read derived from copy
*A* must be free to compete with copies *B*, *C*… and wins at most one gene. Only primary
alignments are kept (a pair is assigned to a single reference gene), and only the genes
whose copies are all present are written into the bam.

```
soi-wgdi compmap pre-study/merged.alignment.csv pre-study/cds.fa \
    reads/<sample>_filter_R1.fastq.gz reads/<sample>_filter_R2.fastq.gz \
    colnames=pre-study/rename.list ref=<reference> \
    prefix=out/<species>.<reference>/<species>.<reference> ncpu=8
```

Outputs (in `out/<species>.<reference>/`):
`<species>.<reference>.depth.tsv` (mean depth per reference gene per copy),
`.depth.pdf` (depth distribution of the copies — an uneven profile means the copies are
not equally contributed, i.e. probably no shared WGD; an even profile proves nothing by
itself), `.target.cds` (the baits of the retained genes) and `.bam` (primary-only).

**Step 2 – assemble each copy.** HybPiper reuses that bam (`--bwa --no-blast`): every
bait becomes one gene directory, one SPAdes + exonerate run each, i.e. one sequence per
copy. `src/retrieve_copies.py` then collects the per-gene sequences into
`out/<species>.<reference>.cds` and `.pep`, with ids
`<species>.<reference>-<column>_XG<row>` (paralogs get `-0`, `-1`… appended).

```
cd out/<species>.<reference>
reads_first.py -r dummy.fastq -b <species>.<reference>/<species>.<reference>.target.cds \
    --prefix <species>.<reference> --cpu 8 --bwa --no-blast
paralog_investigator.py <species>.<reference>
cd ../.. && python3 src/retrieve_copies.py out/<species>.<reference>
```

**Step 3 – append the sample as new columns.** One column per copy of the reference:

```
soi-wgdi append_cols pre-study/merged.alignment.csv out/<species>.<reference>.cds \
    <species> colnames=pre-study/rename.list prefix=out/<species>.<reference>-appended
```

Outputs: `out/<species>.<reference>-appended.alignment.csv` and `.rename.list` (same
matrix, plus `<species>_<copy>` columns; unmatched cells stay empty).

**Step 4 – subgenome tree.** Only the columns of this reference, of this sample, and the
outgroup go into the tree; gene trees are built from the CDS (codon-aware alignment) and
combined by ASTRAL-hybrid, with an ML tree on the concatenation as a cross-check:

```
bash src/subgenome-tree.sh out/<species>.<reference>-appended.alignment.csv \
    out/<species>.<reference>-appended.rename.list \
    out/<species>.<reference>-tree.pep.faa out/<species>.<reference>-tree.cds.fa \
    out/<species>.<reference>-subgenome 2,3,5,6 Theobroma_cacao 8
```

The last number list is the columns to use (reference copies + sample copies + outgroup).
Shared WGD shows up as each sample copy being sister to the corresponding reference copy.
The script repeats this for all samples of a reference together
(`out/<reference>-subgenome.*`).

### Route 2: de novo (reads -> assembly -> gene models) ###

Run:

```
bash pipeline-denovo.sh design.txt 8
```

**Step 1 – assembly.** `megahit` per sample -> `out/<species>-assembly/final.contigs.fa`.

**Step 2 – gene models.** `metaeuk easy-predict` with the reference proteins as evidence,
then protein and CDS are renumbered with one shared id series (`<species>|g000001`,
`<species>|g000002`…). Outputs `out/<species>.pep`, `.cds`, `.gff`.

**Step 3 – orthology.** `OrthoFinder` on the reference proteins, the outgroup and the
sample proteins (`out/orthofinder-in/`).

**Step 4 – assign genes to copies.** `soi-wgdi orthomap` turns the orthogroups into matrix
columns: every sample gene goes to at most one cell, every cell holds at most one gene,
and the reference genes keep their rows.

```
soi-wgdi orthomap pre-study/merged.alignment.csv out/orthofinder-in/OrthoFinder/Results_* \
    design=design.txt colnames=pre-study/rename.list ref=<reference> \
    out=out/<reference>-orthomap.alignment.csv prefix=out/<reference>-orthomap
```

**Step 5 – subgenome tree.** Same as route 1, on the orthomap matrix.

### Outputs ###

Route 1 writes to `out/`:

```
out/<species>.<reference>.depth.tsv/.depth.pdf   per-copy coverage of the reference genes
out/<species>.<reference>.cds/.pep               assembled copies of this sample
out/<species>.<reference>-appended.*             matrix with this sample's columns
out/<species>.<reference>-subgenome.genetrees    one tree per gene family
out/<species>.<reference>-subgenome.astral-hyb   ASTRAL-hybrid subgenome tree
out/<species>.<reference>-subgenome.alignments*  concatenation (+ .trimal.fasttree = ML tree)
out/<reference>-subgenome.*                      same tree, all samples of one reference
```

Route 2 writes to `out-denovo/` (so the two routes never overwrite each other):

```
out-denovo/<species>-assembly/final.contigs.fa   assembly
out-denovo/<species>.pep/.cds/.gff               gene models (ids <species>|g000001...)
out-denovo/orthofinder-in/OrthoFinder/Results_*/ orthology
out-denovo/<reference>-orthomap.alignment.csv    matrix with the new columns
out-denovo/<reference>-subgenome.*               subgenome tree
```

Trees are plain Newick; for figures use your own viewer (e.g.
[PhyTop](https://github.com/zhangrengang/phytop) or FigTree).

### The example data ###

`pre-study/`, `reads/` and `design.txt` are a small real example, so both routes run in
a few minutes and produce a readable tree:

- reference: *Reevesia pubescens*, the two subgenomes `_R` and `_H`, 58 gene families
  (all of them with both subgenome copies present)
- samples: *Triplochiton scleroxylon* (WGS) and *Ungeria floribunda* (target capture)
- reads: the reads of the contigs that carry the homologous genes of those 58 families
  (all reads of those contigs, nothing sampled away) — 150k read pairs for *Triplochiton*
  (22 MB), 56k for *Ungeria* (7 MB)

```
bash pipeline-refguided.sh design.txt 8
bash pipeline-denovo.sh design.txt 8
```

In this example two subfamilies behave differently, which is what the analysis is for:
*Ungeria* (same subfamily as the reference) puts each of its copies next to the
corresponding reference subgenome, while the two copies of *Triplochiton* (another
subfamily) stay together and attach to one reference subgenome.

### FAQ ###

1. **The design file has no header and no comments.** A comment line would be read as a
   sample, and rows are read by position.
2. **What are `_R` and `_H`?** The two subgenomes of the reference. A polyploid
   lineage keeps several copies of every gene, one per subgenome, so each copy is one
   matrix column named `<species>_<subgenome>`. The labels themselves are free —
   `<species>_A/_B/_C` works the same way; what matters is that every copy column of a
   reference follows `<species>_<subgenome>`.
3. **Where does the sample id have to match a file name?** Only `col1` is used for file
   names (`reads/<col1>_filter_R{1,2}.fastq.gz`); `col2` only labels output.
4. **Few genes assembled (route 1).** Look at `<species>.<reference>.depth.tsv` and
   `.depth.pdf` first: if the copies are covered very unevenly, the sample either has no
   shared WGD or its reads are too shallow/short for that reference. Reads of a target
   capture experiment only cover a few hundred genes.
5. **Re-running.** Every step overwrites `out/`; delete `out/` for a clean run. The
   scripts do not skip finished samples.
6. **One gene family per tree, not per chromosome.** To split the gene trees by
   chromosome, use the GFF and `soi-wgdi chrsplit` on `<prefix>-subgenome.tmp`.
