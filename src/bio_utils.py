from Bio import SeqIO

def calculate_one_contig_size(ref_file):
    input_file = open(ref_file, "r")
    
    resulting_size = 0
    input_reads = 0
    for r in SeqIO.parse(input_file, "fasta"):
        resulting_size += len(r.seq)
        input_reads += 1
    if input_reads > 1:
        raise RuntimeError("Expected only one contig as reference, found " + str(input_reads) + " contigs "
                                "in reference file")
    
    return resulting_size


