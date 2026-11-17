#!/usr/bin/env python3
"""Retrieve per-copy CDS/protein sequences from a HybPiper run.

Ported (2to3) from the mature internal script ~/src/HybPiper.py used by the
original pipeline -- same logic, only the external imports were replaced by a
local logger shim, and the unused MultiHybPiper class was dropped.

HybPiper writes every enriched gene into its own directory:
    <prefix>/<gene>/<prefix>/sequences/FNA/<gene>.FNA   (CDS)
    <prefix>/<gene>/<prefix>/sequences/FAA/<gene>.FAA   (protein)
    <prefix>/<gene>/<prefix>/paralogs/<gene>_paralogs.fasta
This script collects them into <prefix>.cds and <prefix>.pep with the ids
"<prefix>-<gene>" (paralogs: "<prefix>-<gene>-<index>"), which is exactly what
`soi-wgdi append_cols` expects.

usage: python3 retrieve_copies.py <prefix dir>
"""
import os
import sys
import glob
import logging
from Bio import SeqIO

logger = logging.getLogger('retrieve')
logging.basicConfig(level=logging.INFO, format='%(message)s')

SEQUENCES = {
    'nucl': '{prefix}/{gene}/{sample}/sequences/FNA/{gene}.FNA',
    'prot': '{prefix}/{gene}/{sample}/sequences/FAA/{gene}.FAA',
    'introns': '{prefix}/{gene}/{sample}/sequences/intron/{gene}_introns.fasta',
    'supercontig': '{prefix}/{gene}/{sample}/sequences/intron/{gene}_supercontig.fasta',
}


class lazyproperty(object):
    """Minimal stand-in for the external lazy_property package."""

    def __init__(self, f):
        self.f = f
        self.__name__ = f.__name__

    def __get__(self, obj, cls=None):
        if obj is None:
            return self
        val = self.f(obj)
        obj.__dict__[self.__name__] = val
        return val

class HybPiper(object):
	def __init__(self, prefix):
		self.prefix = prefix.rstrip('/')
		self.sample = self.prefix.split('/')[-1]
	@lazyproperty
	def genes_with_seqs(self):
		rows = [ line.strip().split() for line in open(
				'{prefix}/genes_with_seqs.txt'.format(prefix=self.prefix))
				]
		bad = [r for r in rows if len(r) != 2]
		if bad:
			logger.warn('{} malformed line(s) in genes_with_seqs.txt skipped'.format(len(bad)))
		return dict([r for r in rows if len(r) == 2])
	@lazyproperty
	def genes_with_paralog_warnings(self):
		paralog_warnings = '{prefix}/genes_with_paralog_warnings.txt'.format(prefix=self.prefix)
		if os.path.exists(paralog_warnings):
			return {line.strip() for line in open(paralog_warnings) }
		else:
			logger.warn('{} do not exist'.format(paralog_warnings))
			return set([])

	@lazyproperty
	def target_file(self):
		db_file = glob.glob('{prefix}/*.amb'.format(prefix=self.prefix))[0]
		return os.path.splitext(db_file)[0]
	
	def get_paralogs(self, gene):
		return '{prefix}/{gene}/{sample}/paralogs/{gene}_paralogs.fasta'.format(prefix=self.prefix, sample=self.sample, gene=gene)

	@lazyproperty
	def prot_sequences(self):
		return glob.glob('{prefix}/*/{sample}/sequences/FAA/*.FAA'.format(prefix=self.prefix, sample=self.sample))
	@lazyproperty
	def nucl_sequences(self):
		return glob.glob('{prefix}/*/{sample}/sequences/FNA/*.FNA'.format(prefix=self.prefix, sample=self.sample))
	@lazyproperty
	def intron_sequences(self):
		return glob.glob('{prefix}/*/{sample}/sequences/intron/*introns.fasta'.format(prefix=self.prefix, sample=self.sample))
	@lazyproperty
	def supercontig_sequences(self):
		return glob.glob('{prefix}/*/{sample}/sequences/intron/*supercontig.fasta'.format(prefix=self.prefix, sample=self.sample))
	def get_sequence(self, gene, seqtype='nucl'):
		seqfile = SEQUENCES[seqtype].format(prefix=self.prefix, sample=self.sample, gene=gene)
		if os.path.exists(seqfile):
			return seqfile
		else:
			return None
	def get_genename(self, path):
		last = False
		for subpath in path.split('/'):
			if last:
				return subpath
			if subpath == self.prefix:
				last = True
	@lazyproperty
	def genes(self):
		return [path.split('/')[-2] for path in \
					glob.glob('{prefix}/*/'.format(prefix=self.prefix))]
	@lazyproperty
	def cdsfile(self):
		return '{}.cds'.format(self.prefix)
	@lazyproperty
	def pepfile(self):
		return '{}.pep'.format(self.prefix)
	def retrieve_sequences(self, ):
		prefix = self.prefix
		logger.info('retrieving CDS and protein sequences from {}'.format(prefix))
		cdsfile = open(self.cdsfile, 'w')
		pepfile = open(self.pepfile, 'w')
		for gene in sorted(self.genes):
			if not gene in self.genes_with_seqs:
				continue
			for seqtype, outfile in zip(['nucl', 'prot'], [cdsfile, pepfile]):
				seqfile = self.get_sequence(gene, seqtype=seqtype)
				rc = SeqIO.read(seqfile, 'fasta')
				rc.id = '{}-{}'.format(self.sample, gene)
				SeqIO.write(rc, outfile, 'fasta')
			if gene in self.genes_with_paralog_warnings: # paralog
				paralog = self.get_paralogs(gene)
				if not os.path.exists(paralog):
					#logger.info('{} do not exists, skipping'.format(paralog))
					continue
				for rc in SeqIO.parse(paralog, 'fasta'):
					index = rc.id.split('.')[-1]
					if index == 'main':
						continue
					rc.id = '{}-{}-{}'.format(self.sample, gene, index)
					SeqIO.write(rc, cdsfile, 'fasta')
					rc.seq = aa = rc.seq.translate()
					SeqIO.write(rc, pepfile, 'fasta')
		cdsfile.close()
		pepfile.close()


def main():
    if len(sys.argv) != 2:
        sys.exit('usage: python3 retrieve_copies.py <prefix dir>')
    HybPiper(sys.argv[1]).retrieve_sequences()


if __name__ == '__main__':
    main()
