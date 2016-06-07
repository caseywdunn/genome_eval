import os
import argparse
import sys
from Bio import SeqIO

def internal_stops(proteins):
	'''
	Finds if protein has internal stops and prints seq_id
	to a file
	'''
	
	infile = open(proteins, "rU")
	filename = os.path.basename(proteins).split("_")[0]
	outfile=os.path.join(os.getcwd(), filename + "_ids_internal_stops.txt")
	with open(outfile, "w") as f:
		for rec in SeqIO.parse(infile, "fasta"):
			stop = rec.seq.find("*")
			if stop > 0 and stop != len(rec.seq)-1:
				print >>f, rec.id

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes a fasta file \
            of amino acid sequences (proteins), finds all internal stop \
            codons and prints a file with sequence_ids if sequence has \
            internal stop codons')
    parser.add_argument('proteins', help="""protein file in fasta format""")
    args = parser.parse_args()
    internal_stops(args.proteins)
