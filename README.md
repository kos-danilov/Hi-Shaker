# Hi-Shaker
Hi-Shaker is the tool for working with Hi-C matrices in *.cool format. It can read, change initial data and save results as images in *.png format or edited *.fasta. Repository include test data - example.fasta (randomly generated sequences based on mat18_100k.cool matrix), sample matrix is avalibale via google drive link: https://drive.google.com/file/d/1XPRoiCD53UmwI8CHcYHkGkh1erT5LA9H/view?usp=sharing.

hs.ipynb file contatins example of usage and it is possible to work with data using jupyter notebook. 

hs.py file is a python3 script for working through command line (*python3 hs.py*). After initialization 7 command are avaliable:
- read *path-to-cool-file*;
- move *num-of-contig where*, after command takes two arguments: contig number in actual order to move and exact new index location for this contig (0-based);
- reverse *num-of-contig*, after command takes one argument: contig number to reverse;
- savepng *name-of-the-file i1 i2 j1 j2*, after command takes five arguments: name of the image, indexes in bins coordinates of region of interest (i1,j1 - top left corner, i2,j2 - bottom right corner);
- savefasta *path-to-original-fasta edited-fasta*, after command takes two arguments: path to the original fasta and name of new fasta;
- help;
- exit.
