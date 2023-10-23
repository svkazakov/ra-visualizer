#!/usr/bin/python3

import pysam
from utils import *


def run(args):
    file = args.input_sam

    print2("Reading alignments from " + file)

    samfile = pysam.AlignmentFile(file, 'r')
    out_f = open(args.output_file, 'w')
    
    print("  and writing statistics to " + args.output_file + "...")
    # out_f.write("no.\tshortFN  \tsignal_len\talbacore_len\tmap_strain\tref_name\tref_start\tref_end\tref_len\t" +\
	# 			"  Interesting\tUseful\n")
    out_f.write("read\tread_len\tmap_strain\tref_name\tref_start\tref_end\tmapping_len\t" +
				"Interesting\tUseful\tclipped_head\tclipped_tail\n")

    in_c = 0
    out_c = 0
    all_reads = 0
    mapped = 0
    unmapped = 0
    suppl_alignments = 0
    
    ref_name_one_value = None
    if len(samfile.references) == 1:
        ref_name_one_value = samfile.references[0]
    print("Found one reference in .sam file = " + str(ref_name_one_value is not None))

    use_all_alignments = args.all
    print("use_all_alignments = " + str(use_all_alignments))


    for read in samfile.fetch(until_eof=True):
        in_c += 1
        read_len = read.query_length if read.is_unmapped else read.infer_read_length()
        # print("Found record: - ", end='')

        if read.is_supplementary:
            suppl_alignments += 1
            # print("Supplementary:")
        else:
            all_reads += 1
            # print("Main alignment:")
        
        if (not read.is_supplementary) or use_all_alignments:
            outS = str(read.query_name) + "\t" + str(read_len) + "\t"
            out_c += 1
            
            if read.is_unmapped:
                unmapped += 1
                # print("Unmapped read:  read.is_reverse = " + str(read.is_reverse) + ", ref_start = " + str(read.reference_end) +\
                #         ", ref_end = " + str(read.reference_end))
                outS += "u\t-\t-\t-\t-\t-\t-\t-\t-"
                # print("Unmapped")
            else:
                # print("Mapped:")
                mapped += 1
                mapping_length = read.reference_end - read.reference_start
                ref_name = read.reference_name
                if ref_name_one_value is not None:
                    assert ref_name == ref_name_one_value
                    ref_name = "chr"
                
                # Write your code to determine if the read is interesting here
                interesting = '.'
                
                # Write your code to determine if the read is useful here
                useful = '.'
                
                cigar_list = read.cigartuples
                clipped_head = 0
                clipped_tail = 0
                first_op = cigar_list[0]
                if (first_op[0] == 4) or (first_op[0] == 5):    # clipping operation
                    # print("Clipped " + str(first_op[1]) + " nucs from head")
                    clipped_head = first_op[1]
                last_op = cigar_list[-1]
                if (len(cigar_list) > 1) and ((last_op[0] == 4) or (last_op[0] == 5)):  # clipping from tail
                    # print("Clipped " + str(last_op[1]) + " nucs from tail")
                    clipped_tail = last_op[1]
                
                outS += ('-' if read.is_reverse else '+') + "\t" + ref_name + "\t" + \
                        str(read.reference_start) + "\t" + str(read.reference_end - 1) + "\t" + \
                        str(mapping_length) + "\t" + \
                        interesting + "\t" + useful + "\t" + \
                        str(clipped_head) + "\t" + str(clipped_tail)
                
            out_f.write(outS + "\n")


    samfile.close()
    out_f.close()
    fl = len(int2str(in_c))
    print2("Done")
    print2("Input  records     =  " + int2str(in_c))
    print2("  suppl_alignments =  " + withP(suppl_alignments, in_c, fl=fl))
    print2()
    print2("All  reads       =    " + withP(all_reads, in_c, fl=fl))
    print2("    mapped       =    " + withP(mapped, all_reads, "of all reads", fl=fl))
    print2("  unmapped       =    " + withP(unmapped, all_reads, "of all reads", fl=fl))
    print2()

    print2("Output records   =    " + withP(out_c, in_c, "of input records", fl=fl))


def main():
    parser = argparse.ArgumentParser(description='sam-gen -- Python script to generate <output-table.stats> '
                                                 'from .sam file')

    parser.add_argument("-i", "--input-sam", help="Input .sam file to process", required=True)
    parser.add_argument("-o", "--output-file", help="Output .stats file name", required=True)
    
    parser.add_argument("--all", action='store_true', help="Output all alignments (with supplementary) (default: False)")

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
