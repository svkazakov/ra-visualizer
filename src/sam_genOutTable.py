#!/usr/bin/python3

import sys, pysam
from utils import *


if len(sys.argv) < 3:
	print2("Usage:   sam_process.py  <input.sam>  <output.table>")
	print2("")
	exit(1)

file	= sys.argv[1]

print2("Reading alignments from " + file + "...")

inF  = pysam.AlignmentFile(file, "r")
outF = open(sys.argv[2], "w")
outF.write("# read\tref\tread_len\tmapping_len\n")

inC = 0; outC = 0
allReads = 0; mapped = 0; unmapped = 0
supplAlignments = 0

progress_start()
for read in inF.fetch(until_eof=True):
	inC += 1
	read_len = read.query_length if read.is_unmapped else read.infer_read_length()

	if read.is_supplementary:
		supplAlignments += 1
	else:
		allReads += 1
		if read.is_unmapped:
			unmapped += 1
		else:
			mapped += 1
			outC += 1
			mapping_len = read.reference_end - read.reference_start
			outF.write(read.query_name + "\t" + read.reference_name + "\t" + str(read_len) + "\t" + str(mapping_len) + "\n")
	progress_tick(inC)

progress_finish()


inF.close()
outF.close()
fl = len(int2str(inC))
print2("Done")
print2("Input  records    =  " + int2str(inC))
print2("  supplAlignments =  " + withP(supplAlignments, inC, fl=fl))
print2()
print2("All  reads      =    " + withP(allReads, inC,fl=fl))
print2("    mapped      =    " + withP(mapped, allReads, "of all reads",fl=fl))
print2("  unmapped      =    " + withP(unmapped, allReads, "of all reads",fl=fl))
print2()

print2("Output records    =  " + withP(outC, inC, "of input records", fl=fl))



