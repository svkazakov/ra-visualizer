
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run(args):
    
    input_file = open(args.input_file, "r")
    out_file   = open(args.output_filename, "w")
    
    print("Reading input contigs from " + args.input_file)
    print("  and writing output to " + args.output_filename + "...")
    
    gap_size = args.gap_size
    print("Using gap size = " + str(gap_size))
    
    
    input_reads = 0
    sequences = []
    for r in SeqIO.parse(input_file, "fasta"):
        input_reads += 1
        sequences.append(str(r.seq))
        
    output_contig = ('N' * gap_size).join(sequences)
    record = SeqRecord(Seq(output_contig), id="Concatenated_scaffold", description='')
    SeqIO.write(record, out_file, "fasta")
    
    input_file.close()
    out_file.close()
    
    print("Done!")
    print("Total input contigs =  " + str(input_reads))
    print("Outputted one resulting scaffold")


def main():
    parser = argparse.ArgumentParser(description='concat_contigs.py -- Python script for concat all contigs to one scaffold '
                                                 'with gaps between contigs (\'N\'s). Input and output to .fasta format.')
    
    parser.add_argument("-i", "--input-file", help="Input contigs .fasta file", required=True)
    parser.add_argument("-o", "--output-filename", help="Output file name (.fasta format)", required=True)

    parser.add_argument("-g", "--gap-size", help="Gap size to insert between 2 contigs (default: 10'000)", type=int, default=10000)

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()

