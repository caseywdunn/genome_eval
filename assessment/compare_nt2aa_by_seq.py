import argparse
import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def make_protein_record(nuc_record, trans_table, to_stop):
	'''
	Returns a new SeqRecord with translated sequences
	'''
	if to_stop:
		return SeqRecord(seq=nuc_record.seq.translate(trans_table, to_stop=True),\
			id = nuc_record.id)
	else: 
		return SeqRecord(seq=nuc_record.seq.translate(trans_table), \
			id = nuc_record.id)

def compare_trans2aa(transcripts, trans_table, proteins, to_stop):
	'''
	Compares translated nucleotide sequences with proteins
	'''
	infile = open(transcripts, "rU")
	filename = os.path.basename(transcripts).split("_")
	sp_name = "_".join(filename[0:2])
	proteins = open(proteins, "rU")
	nucleotides = SeqIO.parse(infile, "fasta")
	translation = (make_protein_record(nuc_rec, trans_table, to_stop) for nuc_rec in nucleotides)
	trans_dict = SeqIO.to_dict(translation)
	aa_dict = SeqIO.to_dict(SeqIO.parse(proteins, "fasta"))
	outfile = os.path.join(os.getcwd(), sp_name + "_mismatch_ids.txt")
	with open(outfile, "w") as f:
		for key, value in trans_dict.iteritems():
			if str(value.seq) not in [str(v.seq) for v in aa_dict.values()]:
				print >>f, '%s' % key
			else:
				continue
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Takes a fasta file of \
		nucleotides, translates it, and compares it to a fasta file \
		of proteins.')
	parser.add_argument('transcripts', help="""transcripts file in fasta \
		format""")
	parser.add_argument('proteins', help="""proteins file in fasta format.""")
	parser.add_argument('--trans_table', type=int, help="""translation table
		[default = 1] (Universal)""", default=1)
	parser.add_argument('--to_stop', help="""translation should NOT include \
	 the stop codon (if present in the transcripts or proteins) \
	 [default = 0]""", default=0)
	args = parser.parse_args()

	print '''
        ====================================================
        +                                                  +    
        +  THERE MAY BE A BioPython WARNING MESSAGE BELOW, +
        +       DO NOT WORRY, THE SCRIPT WILL WORK         +
        +                                                  +
        ====================================================
        '''
	compare_trans2aa(args.transcripts, args.trans_table, args.proteins, args.to_stop)
