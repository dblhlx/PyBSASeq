**Note**: 
1. It is strongly recommended to have [fisher](https://github.com/brentp/fishers_exact_test) installed on your system. It is fast in dealing with large datasets.
2. If you used PyBSASeq in your pulications, please cite: Zhang, J., Panthee, D.R. PyBSASeq: a simple and effective algorithm for bulked segregant analysis with whole-genome sequencing data. BMC Bioinformatics 21, 99 (2020). https://doi.org/10.1186/s12859-020-3435-8

### PyBSASeq
A novel algorithm for BSA-Seq data analysis

Python 3.6 or above is required to run the script.

#### Usage
If only the bulk genomes were sequenced, use the script in the folder /BulksOnly; if the genomes of both the bulks and the parents were sequenced, use the script in the folder /Parents&Bulks.

#### Workflow
1. SNP filtering
2. Perform Fisher's exact test using the AD values of each SNP in both bulks. A SNP would be identified as an sSNP if its p-value is less than alpha. In the meantime, simulated REF/ALT reads of each SNP is obtained via simulation under null hypothesis, and Fisher's exact test is also performed using these simulated AD values; for each SNP, it would be an sSNP if its p-value is less than smalpha. Identification of sSNPs from the simulated dataset is for threshold calculation. A file named "COMPLETE.txt" will be writen to the working directory if Fisher's exact test is successful, and the results of Fisher's exact test are saved in a .csv file. The "COMPLETE.txt" file needs to be deleted in case starting over is desired. 
3. Threshold calculation. The result is saved in the "threshold.txt" file. The "threshold.txt" file needs to be deleted if starting over is desired (e.g, if the size of the sliding window is changed).
4. Plotting.

#### Dataset
The file [snp_final.tsv.bz2](https://github.com/dblhlx/PyBSASeq/blob/master/snp_final.tsv.bz2) contains the [GATK4](https://software.broadinstitute.org/gatk/download/)-generated SNP dataset using the sequencing data from [the work of Yang et al](https://www.ncbi.nlm.nih.gov/pubmed/23935868). The sequence reads were treated as either single-end or paired-end when aligned to the reference genome. Significantly more SNPs were identified by the latter approach; however, the results of BSA-Seq analysis were very similar. Only the tsv file generated by the former approach is included here because of the 25 Mb file size limitation. 

A [small test file](https://github.com/dblhlx/PyBSASeq/blob/master/smallTestFile.tsv) is included here for testing purpose; just issue the command `python PyBSASeq.py` in a terminal to test the Python script.

#### Other methods for BSA-Seq data analysis
BSA-Seq data analysis can be done using either the [SNP index (allele frequency) method](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105) or the [G-statistic method](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255) as well. I implemented both methods in Python for the purpose of comparison: [PySNPIndex](https://github.com/dblhlx/PySNPIndex) and [PyGStatistic](https://github.com/dblhlx/PyGStatistic). The Python implementation of the original [G-statistic method](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255) by Magwene can be found [here](https://bitbucket.org/pmagwene/bsaseq/src/master/) (just found this site today, 6/27/2019), and the R implementation of both methods by Mansfeld can be found [here](https://github.com/bmansfeld/QTLseqr).
