**Note**: 
It is strongly recommended to have [fisher](https://github.com/brentp/fishers_exact_test) installed on your system. It is fast in dealing with large datasets.

**Citation**:
If you used PyBSASeq in your pulications, please cite:
- Zhang, J., Panthee, D.R. PyBSASeq: a simple and effective algorithm for bulked segregant analysis with whole-genome sequencing data. BMC Bioinformatics 21, 99 (2020). https://doi.org/10.1186/s12859-020-3435-8
- Zhang, J., Panthee, D.R. Next-generation sequencing-based bulked segregant analysis without sequencing the parental genomes. G3 Genes|Genomes|Genetics, jkab400 (2021), https://doi.org/10.1093/g3journal/jkab400


### PyBSASeq
A novel algorithm with high detection power for BSA-Seq data analysis - the significant structural variant method

Python 3.6 or above is required to run the script.

#### Usage

`$ python PyBSASeq.py -i input -o output -p popstrct -b fbsize,sbsize`

Here are the details of the options used in the script:
- `input` – the names of the input files (the GATK4-generated tsv/csv file). If you have variant calling data from both the parents and the bulks, the format is as follows: `parents.tsv,bulks.tsv`. If you have only the variant calling data of the bulks, the format is as follows: `bulks.tsv`. The script `PyBSASeq.py` can handle both situations now. The script and the input files should be in the same folder.
- `output` – the name of the output file. Default is `BSASeq.csv`
- `popstrct` – population structure; three choices available: F2 for an F<sub>2</sub> population, RIL for a population of recombinant inbred lines, or BC for a backcross population
- `fbsize` – the number of individuals in the first bulk
- `sbsize` – the number of individuals in the second bulk

The default cutoff p-value for identifying significant structural variants (sSV) from the SV dataset is 0.01 (`alpha`), and the default cutoff p-value for identifying sSVs from the simulated dataset is 0.01 (`smalpha`). These values can be changed using the following options:

`-v alpha,smalpha`

`alpha` and `smalpha` should be in the range of 0.0 – 1.0, the chosen value should make statistical sense. The greater the `smalpha` value, the higher the threshold and the lower the false positive rate.

The default size of the sliding window is 2000000 (base pairs) and the incremental step is 10000 (base pairs), and their values can be changed using the following option:

`-s slidingWindowSize,incrementalStep`

Four files (`sliding_windows.csv`, `sv_fagz.csv`, `sv_fagz_fep.csv`, and `threshold.txt`) and a folder in the `date_time` format containing `BSASeq.csv` and `BSASeq.pdf` will be generated in the `./Results` folder after succesfully running the script. If gaps between subplots in `PyBSASeq.pdf` are too wide or too narrow, we can rerun the script to fine-tune the gaps by changing the values of `a` and/or `b` using the options below:

`-a True -g a,b,c,d,e,f`

- `a` - the horizontal gap between subplots
- `b` - the vertical gap between subplots
- `c`, `d`, `e`, and `f` - the top, bottom, left, and right margins of the plot, respectively

The default values for `a`, `b`, `c`, `d`, `e`, and `f` are 0.028, 0.056, 0.0264, 0.054, 0.076, 0.002, 0.002, respectively. The script uses the `sliding_windows.csv` and `threshold.txt` files generated previously for plotting, and it is very fast.

If two or more peaks/valleys and all the values in between are beyond the confidence intervals/thresholds, only the highest peak or the lowerest valley will be identified as the peak/valley of this region. We can rerun the script to identify the positions of the other peaks/valleys and test their significance using the option below:

`-e a1,b1,c1,a2,b2,c2,......,an,bn,cn`

- `a` - chromosome id
- `b` - the start point of a chromosomal fragment
- `c` - the end point of a chromosomal fragment

Right now, this option will not work if the chromosome IDs in the reference genome sequences are not digits, with the exception of sex chromosomes; we can use 1000 - 1005 to respectively represent sex chromosomes X, Y, Z, W, U, and V when specify regions on these chromosomes. Multiple regions on the same chromsome can be selected.

#### Workflow
1. Structural variant (SV) filtering
2. Perform Fisher's exact test using the AD values of each SV in both bulks. A SV would be identified as an sSV if its p-value is less than `alpha`. In the meantime, simulated REF/ALT reads of each SV is obtained via simulation under null hypothesis, and Fisher's exact test is also performed using these simulated AD values; for each SV, it would be an sSV if its p-value is less than `smalpha`. Identification of sSVs from the simulated dataset is for threshold calculation.
3. Threshold calculation. The result is saved in the `threshold.txt` file. The `threshold.txt` file needs to be deleted if starting over is desired (e.g, if the size of the sliding window is changed).
4. Plotting.

#### Test datasets
Two SV datasets for testing purpose, one from rice while the other from maize, are stored in the [Data](https://github.com/dblhlx/PyBSASeq/tree/master/Data) folder. The rice .csv files were generated via [GATK4](https://software.broadinstitute.org/gatk/download/) using the sequencing data from [the work of Lahari _et al_](https://www.ebi.ac.uk/ena/browser/view/PRJEB27629) while the maize .csv files were extracted from .csv/.tsv files generated by [Zheng _et al_](https://doi.org/10.1534/g3.120.401192).

#### Bug report

Please report [here](https://github.com/dblhlx/PyBSASeq/issues) if encounter any issue. I can receive issue notifications now.

#### Other methods for BSA-Seq data analysis
BSA-Seq data analysis can be done using either the [SNP index (allele frequency) method](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105) or the [G-statistic method](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255) as well. I implemented both methods in Python for the purpose of comparison: [PySNPIndex](https://github.com/dblhlx/PySNPIndex) and [PyGStatistic](https://github.com/dblhlx/PyGStatistic). The Python implementation of the original [G-statistic method](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255) by Magwene can be found [here](https://bitbucket.org/pmagwene/bsaseq/src/master/) (just found this site today, 6/27/2019), and the R implementation of both methods by Mansfeld can be found [here](https://github.com/bmansfeld/QTLseqr).
