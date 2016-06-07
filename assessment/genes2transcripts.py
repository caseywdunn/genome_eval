#!/usr/bin/env python

import sys
import os
import string
from Bio import SeqIO

Usage = """
Takes a fasta file corresponding to genes where
upper case corresponds to exons or trasncribed regions
and lower case or "n" corresponds to non-coding regions,
and removes lower case and n entries from the sequences.

Usage:
	python genes2transcripts.py gene_models.fasta
"""

if len(sys.argv) < 1:
	print Usage
else:
	infile = open(sys.argv[1], "rU")
	filename = os.path.basename(sys.argv[1]).split(".")[0]
	outfile = os.path.join(os.getcwd(), filename + "_transcripts.fna")
	with open(outfile, "w") as f:
		for record in SeqIO.parse(infile, "fasta"):
			interm = str(record.seq)
			result = interm.translate(None, "N")
			result = result.translate(None, string.ascii_lowercase)
			print >>f, '>%s' % record.id
			print >>f, result
