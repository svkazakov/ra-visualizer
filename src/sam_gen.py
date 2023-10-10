#!/usr/bin/python3

import sys, pysam
from utils import *


def run(args):
    file = args.input_sam

    print2("Reading alignments from " + file + "...")

    samfile = pysam.AlignmentFile(file, "r")
    out_f = open(args.output, "w")
    
    print("writting statistics to file " + args.output)
    out_f.write("no.\tshortFN  \tsignal_len\talbacore_len\tmap_strain\tref_name\tref_start\tref_end\tref_len\t" +\
				"  Interesting\tUseful\n")

    in_c = 0
    out_c = 0
    all_reads = 0
    mapped = 0
    unmapped = 0
    suppl_alignments = 0

    progress_start()
    for read in samfile.fetch(until_eof=True):
        in_c += 1
        read_len = read.query_length if read.is_unmapped else read.infer_read_length()

        if read.is_supplementary:
            suppl_alignments += 1
        else:
            all_reads += 1
            if read.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
                out_c += 1
                mapping_len = read.reference_end - read.reference_start
                
                outS = str(readInfo.no) + "\t" + readInfo.file + "   \t" + str(readInfo.signal_len) + "\t" + str(
                    readInfo.albacore_len) + "\t"
                
                if read is None:
                    # print("Read without mapping information!!  shortFN = " + str(readInfo.file))
                    outS += "noInfo\t-\t-\t-\t-\t-\t-"
                elif read.is_unmapped:
                    # print("Unmapped read:  read.is_reverse = " + str(read.is_reverse) + ", ref_start = " + str(read.reference_end) +\
                    #         ", ref_end = " + str(read.reference_end))
                    outS += "u\t-\t-\t-\t-\t-\t-"
                else:  # i.e. mapped
                    mapping_length = read.reference_end - read.reference_start
                    ref_name = read.reference_name
                    if refName_oneValue is not None:
                        assert ref_name == refName_oneValue
                        ref_name = "chr"
                    
                    outS += ('+' if not read.is_reverse else '-') + "\t" + ref_name + "\t" + \
                            str(read.reference_start) + "\t" + str(read.reference_end - 1) + "\t" + str(
                        mapping_length) + "\t" + \
                            readInfo.interesting_by_coords + "\t" + readInfo.useful_by_coords
                
                out_f.write(read.query_name + "\t" + read.reference_name + "\t" + str(read_len) + "\t" + str(
                    mapping_len) + "\n")
        progress_tick(in_c)

    progress_finish()

    samfile.close()
    out_f.close()
    fl = len(int2str(in_c))
    print2("Done")
    print2("Input  records    =  " + int2str(in_c))
    print2("  suppl_alignments =  " + withP(suppl_alignments, in_c, fl=fl))
    print2()
    print2("All  reads      =    " + withP(all_reads, in_c, fl=fl))
    print2("    mapped      =    " + withP(mapped, all_reads, "of all reads", fl=fl))
    print2("  unmapped      =    " + withP(unmapped, all_reads, "of all reads", fl=fl))
    print2()

    print2("Output records    =  " + withP(out_c, in_c, "of input records", fl=fl))


def main():
    parser = argparse.ArgumentParser(description='sam-gen -- Python script to generate <output-table.stats> '
                                                 'from sam mapping')

    parser.add_argument("-i", "--input-sam", help="Input sam file to process")
    parser.add_argument("-o", "--output", help="Output file name <output.stats>")

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
