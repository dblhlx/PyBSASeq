# PyBSASeq
A Python script for BSA-Seq data analysis

Python 3.6 or above is required to run the script. 

Usage:
> python PyBSASeq.py -i input -o output -p popstrct -f fbsize -s sbsize

Here are the details of the options used in the script:
* input – the name of the input file (the GATK4-generated tsv file)
* output – the name of the output file
* popstrct – population structure; three choices available: F2 for an F2 population, RIL for a population of recombinant inbred lines, or BC for a backcross population
* fbsize – the number of individuals in the first bulk
* sbsize – the number of individuals in the second bulk

The default p-value for identifying ltaSNPs from the SNP data set is 0.01 (alpha), the default p-value for identifying ltaSNPs from the simulation sample is 0.1 (smalpha), and the default size of the simulation sample is 300 (smplesize). These values can be changed using the following options: 
> --alpha p1 --smalpha p2 --smplsize size 
p1 and p2 should be in the range of 0.0 – 1.0, and size can be any positive integer. The greater the p2 value, the higher the threshold and the lower the false positive rate.
