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
* TeX Live

Python packages:
* pyx
* pysam

## For Windows
Install Linux (for example, Ubuntu 22.04.2 LTS) with WSL for Windows (you can use any instruction on Web 
how to do this, for example, see first steps in [this article](https://ruslanmv.com/blog/Python3-in-Windows-with-Ubuntu)).

## Questions, feedback
You can make issues on the Project page on GitHub or write me directly to svkazakov.me at gmail.com.

## License
The MIT License (MIT)
