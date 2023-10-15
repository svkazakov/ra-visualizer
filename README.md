# ra-visualizer
Genome read alignments visualizer

**ra-visualizer** is a Python script for visualizing genome reads and/or contigs alignments to reference.

It can generate the following pictures (in .svg format), with arbitrary zoom and offset, relative to the start 
of the reference:<br/>
Image with reads coloring:<br/>
<img src="pics/pic-useful-reads.jpg" alt="Image with useful reads coloring" width="650"> <br/>
<img src="pics/pic1.png" alt="Image 1" width="650"><br/>
<img src="pics/pic2.jpg" alt="Image 2" width="650"><br/>
<img src="pics/pic4.jpg" alt="Image 4" width="650"><br/>
Only contig mapping output:<br/>
<img src="pics/pic-contigs-mapping.jpg" alt="Image with contigs mapping" width="650"><br/>


Now the project is still a bit raw, but within a few months I will try to improve it.


## Future improvements

* Tidy up the code, add examples of its work.
* Add support for working directly with .sam files.

## Dependencies
* TeX Live <br>
  (`sudo apt install texlive`)
* bwa

Python packages:
* pyx
* pysam

<b>NB.</b> Python package `pysam` is working only for Linux OS, but for Windows users there are several ways
how to make everything work.  One of them is described in the next section.

## For Windows
Install Linux (for example, Ubuntu 22.04.2 LTS) with WSL for Windows (you can use any instruction on Web 
how to do this, for example, see first steps in [this article](https://ruslanmv.com/blog/Python3-in-Windows-with-Ubuntu)).<br>
After that you can install `python3` and all dependencies for it (don't forget to install texlive).<br>
After all, you can run `python3 src/visualize.py` in folder `/mnt/<PATH_TO_PROJECT_DIR>`.

## Run steps
Suppose you have:
1. Reference in `.fasta` format
2. Some assembly (contigs) in `.fasta` format
3. And long reads in `.fastq` format...  

What you should do to visualize it all together? ☺

Run these steps:
1. `bwa index reference.fasta`
2. `bwa mem reference.fasta contigs.fasta > contigs.sam`
3. `bwa mem -x ont2d -A1 -B2 -O2 -E1 -t <number_of_CPUs_to_use> reference.fasta reads.fastq > reads.sam` <br>
  **NB.** This command is used for Oxford Nanopore long reads.  For other sequencing technologies you should use
    different command to align reads to reference.
4. (temporary) Generate `.stats` file from `contigs.sam` file via command: <br>
  `python3 src/gen_stats_table.py -i contigs.sam -o contigs.stats`
5. (temporary) Generate `.stats` file from `reads.sam` file via command: <br>
  `python3 src/gen_stats_table.py -i reads.sam -o reads.stats`
6. Run Visualizer: <br>
  `python3 src/visualize.py --ref-size <reference_size> -c contigs.stats -r reads.stats`
7. See the resulting `output.svg` image and/or change the output options in the previous step


## Example 1
ToDo

## Questions, feedback
You can make issues on the Project page on GitHub or write me directly to svkazakov.me at gmail.com.

## License
The MIT License (MIT)
