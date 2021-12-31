**Note**: 
It is strongly recommended to have [fisher](https://github.com/brentp/fishers_exact_test) installed on your system. It is fast in dealing with large datasets.

**Citation**
If you used PyBSASeq in your pulications, please cite:
- Zhang, J., Panthee, D.R. PyBSASeq: a simple and effective algorithm for bulked segregant analysis with whole-genome sequencing data. BMC Bioinformatics 21, 99 (2020). https://doi.org/10.1186/s12859-020-3435-8
- Zhang, J., Panthee, D.R. Next-generation sequencing-based bulked segregant analysis without sequencing the parental genomes. G3 Genes|Genomes|Genetics, jkab400, https://doi.org/10.1093/g3journal/jkab400


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
Two SNP datasets for testing purpose, one for the bulks while the other for the parents, are stored in the [Data](https://github.com/dblhlx/PyBSASeq/tree/master/Data) folder. Both .csv files were generated via [GATK4](https://software.broadinstitute.org/gatk/download/) using the sequencing data from [the work of Lahari et al](https://www.ebi.ac.uk/ena/browser/view/PRJEB27629).

#### Other methods for BSA-Seq data analysis
BSA-Seq data analysis can be done using either the [SNP index (allele frequency) method](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105) or the [G-statistic method](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255) as well. I implemented both methods in Python for the purpose of comparison: [PySNPIndex](https://github.com/dblhlx/PySNPIndex) and [PyGStatistic](https://github.com/dblhlx/PyGStatistic). The Python implementation of the original [G-statistic method](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255) by Magwene can be found [here](https://bitbucket.org/pmagwene/bsaseq/src/master/) (just found this site today, 6/27/2019), and the R implementation of both methods by Mansfeld can be found [here](https://github.com/bmansfeld/QTLseqr).
