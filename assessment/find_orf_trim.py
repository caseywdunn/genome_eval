import argparse
import sys
import os
from Bio import SeqIO

def find_orfs(seq, trans_table, min_protein_length, rev_comp):
    '''
    Finds all ORFs in a given transcript. 
    '''
    all_orfs = []
    seq_len = len(seq)
    if rev_comp:
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                trans = str(nuc[frame:].translate(trans_table))
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if aa_end-aa_start >= min_protein_length:
                        if strand == 1:
                            start = frame+aa_start*3
                            end = min(seq_len,frame+aa_end*3+3)
                        else:
                            start = seq_len-frame-aa_end*3-3
                            end = seq_len-frame-aa_start*3
                        all_orfs.append((len(trans[aa_start:aa_end]), start, end,
                                        strand, trans[aa_start:aa_end]))
                    aa_start = aa_end+1
    else:
        for strand, nuc in [(+1, seq)]:
            for frame in range(3):
                trans = str(nuc[frame:].translate(trans_table))
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if (aa_end-aa_start)+1 >= min_protein_length:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                        all_orfs.append((len(trans[aa_start:aa_end]), start, end,
                                    strand, trans[aa_start:aa_end]))
                    aa_start = aa_end+1
    all_orfs.sort(reverse=True)
    return all_orfs

def find_full_orfs(orfs):
    '''
    Takes a dict of ORFs, looks for longest full ORF (Met-*)
    and trims nucleotide sequence.
    '''
    result =[]
    for i in range(len(orfs)):
        length = orfs[i][0]
        start = orfs[i][1]
        end = orfs[i][2]
        strand = orfs[i][3]
        aa_seq = orfs[i][4]
        new_start = start
        if aa_seq.find("M") != -1:
            new_start = start+aa_seq.find("M")*3
            aa_seq = aa_seq[aa_seq.find("M"):end]
        else:
            new_start = start
            aa_seq = aa_seq[new_start:end]
        result.append((len(aa_seq), new_start, end, strand, aa_seq))
    result.sort(reverse=True)
    return result

#TO DO: This function should filter out full_orf shorter than min_prot_length

def decide_final_orf(seq_name, orfs, min_protein_length):
    '''
    Chooses the final ORF the original transcript will be trimmed to,
    and returns a dict sequence_name : final_ORF
    '''
    final_orf = {}
    if len(orfs) != 0:
        final_orf[seq_name] = orfs[0]
    else:
        final_orf[seq_name] = ()
        print "%s does not have any translation longer than %d aa" % (seq_name, min_protein_length)
    return final_orf

def write_orfs(transcripts, trans_table, min_protein_length, rev_comp, full_longest_orfs):
    '''
    Writes output file of trimmed transcripts
    '''
    infile = open(transcripts, "rU")
    filename = os.path.basename(transcripts).split(".")[0]
    outfile = os.path.join(os.getcwd(), filename + "_trimmed.fna")
    with open(outfile, "w") as f:
        for record in SeqIO.parse(infile, "fasta"):
            orfs = find_orfs(record.seq, trans_table, min_protein_length, rev_comp)
            if full_longest_orfs:
                full_orfs = find_full_orfs(orfs)
                orf_to_trim = decide_final_orf(record.id, full_orfs, min_protein_length)
                for k,v in orf_to_trim.iteritems():
                    if len(v) != 0:
                        print >>f, '>%s' % k
                        print >>f, record.seq[v[1]:v[2]]
            else:
                orf_to_trim = decide_final_orf(record.id, orfs, min_protein_length)
                for k,v in orf_to_trim.iteritems():
                    if len(v) != 0:
                        print >>f, '>%s' % k
                        print >>f, record.seq[v[1]:v[2]]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes a fasta file \
            of nucleotide sequences (transcripts), finds all ORFs, and \
            trims original sequence to the longest full ORF (Met-*)\
            or longest ORF regardless of starting amino acid. Writes \
            a new fasta `_trimmed.fa` file.')
    parser.add_argument('transcripts', 
            help="""transcripts file in fasta format""")
    parser.add_argument('--trans_table', type=int, help="""translation table 
            [default = 1] (i.e., Universal)""", default=1)
    parser.add_argument('--min_prot_length', type=int, help="""minimum 
            protein length [default = 1] (i.e., all proteins)""", default=1)
    parser.add_argument('--rev_comp', type=int, help="""translate sequences
            in 5`-3` AND 3`-5` directions [default = 0] (i.e, only 5`-3`)""", 
            default=0)
    parser.add_argument('--full_longest_orfs', type=int, help="""generate only full
            longest ORFs [default = 1] (i.e., M-*) or generate just longest ORF""", 
            default=1)
    args = parser.parse_args()
    print '''
    ====================================================
    +                                                  +    
    +  THERE MAY BE A BioPython WARNING MESSAGE BELOW, +
    +       DO NOT WORRY, THE SCRIPT WILL WORK         +
    +                                                  +
    ====================================================
    '''
    write_orfs(args.transcripts, args.trans_table, args.min_prot_length, args.rev_comp, args.full_longest_orfs)
