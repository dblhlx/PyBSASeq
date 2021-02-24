#### Usage

`$ python PyBSASeq_WP.py -i input -o output -p popstrct -b fbsize,sbsize`

Here are the details of the options used in the script:
- input – the names of the input files (the GATK4-generated tsv file). For the script PyBSASeq_WP.py, the format is as follows: parents.tsv,bulks.tsv
- output – the name of the output file
- popstrct – population structure; three choices available: F2 for an F2 population, RIL for a population of recombinant inbred lines, or BC for a backcross population
- fbsize – the number of individuals in the first bulk
- sbsize – the number of individuals in the second bulk

The default cutoff p-value for identifying significant SNPs (sSNP) from the SNP dataset is 0.01 (alpha), and the default cutoff p-value for identifying sSNPs from the simulated dataset is 0.1 (smalpha). These values can be changed using the following options:

`-v alpha,smalpha`

alpha and smalpha should be in the range of 0.0 – 1.0, the chosen value should make statistical sense. The greater the smalpha value, the higher the threshold and the lower the false positive rate.

The default size of the sliding window is 2000000 (base pairs) and the incremental step is 10000 (base pairs), and their values can be changed using the following option:

`-s slidingWindowSize,incrementalStep`