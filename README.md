# ZGA - prokaryotic genome assembly and annotation pipeline

[![version status](https://img.shields.io/pypi/v/zga.svg)](https://pypi.python.org/pypi/zga)
[![Anaconda Cloud](https://anaconda.org/bioconda/zga/badges/installer/conda.svg)](https://anaconda.org/bioconda/zga/)
[![Publication](https://img.shields.io/badge/DOI-published-green.svg)](https://dx.doi.org/10.1101/2021.04.27.441618)

## Main Features

* Wide range of supported reads: Illumina, Oxford Nanopore, PacBio, BGI.
* Short read multi-threaded processing: QC, filtering, trimming, overlapped pairs merging.
* Assemblies from short reads, long reads or hybrid assembly using modern and powerful assemblers: [SPAdes](http://cab.spbu.ru/software/spades/), [Unicycler](https://github.com/rrwick/Unicycler/) or [Flye](https://github.com/fenderglass/Flye).
* Quality control of assembly: completeness and contamination assessment with [CheckM](https://github.com/Ecogenomics/CheckM) as well as PhiX detection.
* Fast annotation of bacterial and archeal genome assemblies with [DFAST](https://github.com/nigyta/dfast_core) .
* No High Performance Computing needed. The pipeline works on laptop or desktop.

## Installation

ZGA is written in Python and tested with Python 3.6 and Python 3.7.

### Install with conda

[![Anaconda latest release](https://anaconda.org/bioconda/zga/badges/latest_release_date.svg)](https://anaconda.org/bioconda/zga/)

The simplest way to install ZGA and all dependencies is **conda**:

1. You need to install conda, e.g. [**miniconda**](https://conda.io/en/latest/miniconda.html). Python 3.7 is preferred.

2. After installation You should add channels - the conda's software sources:  
`conda config --add channels bioconda`  
`conda config --add channels conda-forge`

3. At the end You should install ZGA to an existing active environment (Python 3.6 or 3.7):  
`conda install zga`  
or create a fresh environment and activate it:  
`conda create -n zga zga`  
`conda activate zga`

If You have troubles with bioconda channel try to use my personal channel https://anaconda.org/laxeye/zga `conda install -c laxeye zga`

### Install from PyPI

[![Downloads](https://pepy.tech/badge/zga/month)](https://pypi.python.org/pypi/zga)

Run `pip install zga`. Biopython is the only one dependency installed from PyPI. All other dependencies You should install manually or using **conda** as mentioned above. CheckM is available on **PyPi**, but it's easier to install it using **conda**.

### Get source from Github

You can get ZGA by cloning from the repository with `git clone https://github.com/laxeye/zga.git` or by downloading an archive. After downloading enter the directory `cd zga` and run `python3 setup.py install`.

Don't forget to install dependecies (see bellow).

### Installing dependencies

ZGA uses several software and libraries including:

* [fastp](https://github.com/OpenGene/fastp)
* [BBmap](https://sourceforge.net/projects/bbmap/)
* [NxTrim](https://github.com/sequencing/NxTrim)
* [mash](https://mash.readthedocs.io/en/latest/)
* [SPAdes](http://cab.spbu.ru/software/spades/) (>= 3.12 to support merged paired-end reads, >= 3.5.0 to support Nanopore reads)
* [Unicycler](https://github.com/rrwick/Unicycler/)
* [Flye](https://github.com/fenderglass/Flye) >= 2.6
* [minimap2](https://github.com/lh3/minimap2/)
* [racon](https://github.com/lbcb-sci/racon)
* [CheckM](https://github.com/Ecogenomics/CheckM) >= 1.1.0
* [BioPython](https://biopython.org/)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [DFAST](https://github.com/nigyta/dfast_core)

You may install all dependencies separately using **conda**. It's highly recommended to create a new conda environment:

`conda create -n zga "python>=3.6" fastp "spades>=3.12" unicycler checkm-genome dfast bbmap blast biopython"nxtrim "mash>=2" flye minimap2 racon "samtools>=1.9"`

and activate it

`conda activate zga`

Otherwise you may install dependencies to existing conda environment:

`conda install "python>=3.6" fastp "spades>=3.12" unicycler checkm-genome dfast bbmap blast biopython nxtrim "mash>=2" flye minimap2 racon "samtools>=1.9"`

Of course, it's possible to use *another ways* even compile all tools from source code. In this case you should check if binaries are in your '$PATH' variable.

#### DFAST database download

After installation DFAST downloads some basic databases. It's recommended to download more databases using *dfast_file_downloader.py* command-line script:

Run `dfast_file_downloader.py -h` to see available databases and options.

Default databases may be donloaded with `dfast_file_downloader.py --protein dfast --cdd Cog --hmm TIGR`

Soon auto-download of databases during installation with conda will be available.

If You want to use more databases You need to edit DFAST configuration file and provide it to ZGA with `--dfast-config` option.

### Operating systems requirements

ZGA was tested on Ubuntu 18.04, 19.10 and 20.04. Most probably any modern 64-bit Linux distribuition is suitable.

Your feedback on other OS is welcome!

## Usage

Run `zga -h` to get a help message.

### Pipeleine steps

ZGA includes several steps:

1. Read quality check ('readqc')
2. Read processing ('preprocessing')
3. Genome assembling ('assembling')
4. Genome polishing ('polishing')
5. Genome quality assessment ('check_genome')
6. Genome annotation ('annotation')

You may start from any step and finish at any step providing arguments `-s` or `--first-step` and `-l` or `--last-step` followed by step designation (in brackets in the list above).

E.g. if You like to perform read processing, genome assembling and genome polishing You should run

`zga --first-step preprocessing --last-step polishing ...`

### Input files

ZGA may use unprocessed or processed sequencing reads from different platforms as well as genome assemblies to perform assembly polishing, assembly quality assessment and assembly annotation. FASTQ format gzipped or not is required for sequencing reads. Paired-end reads shoul be provided in separate files, not interleaved. Sequencing reads should be provided as space separated list after corresponding argument:

`-1` or `--pe-1` for forward paired-end reads (Illumina, BGI)  
`-2` or `--pe-2` for reverse paired-end reads  
`-S` or `--single-end` for unpaired short reads  
`--pe-merged` for merged overlapping paired-end reads (if You performed merging earlier)  
`--mp-1` for first mate-pair reads, RF orientation is supposed  
`--mp-2` for second mate-pair reads  
`--pacbio` for PacBio single-end sequencing reads  
`--nanopore` for Oxford Nanopore sequencing reads  

When `bbduk.sh` (short read trimming tool) throws an exception ZGA tries to repair reads with `repair.sh` (from BBMap).

#### Examples

`zga -1 Raw.R1.fq.gz -2 Raw.R2.fq.gz` unprocessed paired-end reads  
`zga -1 Unmerged_1.fq -2 Unmerged_2.fq --pe-merged Merged.fq` reads after processing (overlapping reads merging)  
`zga -1 Lib1.R1.fq.gz Lib2.R1.fq -2 Lib1.R2.fq Lib2.R2.fq` combination of reads from two sequencing libraries  

### Output

ZGA produces up to 4 sub-folders in output folder:

* **readQC** - results of reaq quality control with *fastp*,
* **reads** - processed reads,
* **assembly** - folder produced by genomic assembler,
* **annotation** - annotated genome.

Log-file *zga.log* is available in the output folder.

### Usage examples

Perform all steps: read qc, read trimming and merging, assembly, CheckM assesment with default (bacterial) marker set, DFAST annotation and use 4 CPU threads where possible:

`zga -1 R1.fastq.gz -2 R2.fastq.gz --bbmerge --threads 4 -o my_assembly`

Assemble with SPAdes using paired-end and nanopore reads of archaeal genome (CheckM will use archaeal markers) altering memory limit to 16 GB:

`zga -1 R1.fastq.gz -2 R2.fastq.gz --nanopore MiniION.fastq.gz -a spades --threads 4 --memory-limit 16 --domain archaea -o my_assembly`

*(New in 0.8 development releases)* Short read correction with SPAdes is a computationally expensive step, You may run read-correction with tadpole including 
`--tadpole-correct` option which is much faster and needs less memory.

`zga --tadpole-correct -1 R1.fastq.gz -2 R2.fastq.gz --threads 4 -o my_assembly`

Assemble long reads with Flye skipping long read polishing and perfom short-read polishing with racon:

`zga -1 R1.fastq.gz -2 R2.fastq.gz --nanopore MiniION.fastq.gz -a flye --threads 4 --domain archaea -o my_assembly --flye-short-polish --flye-skip-long-polish`

Assemble from Nanopore reads using unicycler:

`zga -a unicycler --nanopore MiniION.fastq -o nanopore_assembly`

Perform assesment and annotation of genome assembly with e.g. *Pectobacterium* CheckM marker set:

`zga --first-step check_genome -g pectobacterium_sp.fasta --checkm_rank genus --checkm_taxon Pectobacterium -o my_output_dir`

Let CheckM to infer the right marker set:

`zga --first-step check_genome -g my_genome.fa --checkm_mode lineage -o my_output_dir`

## Known issues and limitations

ZGA is in the stage of active development.

Known issues and limitations:

* Unicycler can't use mate-pair reads or multiple libraries of same type.

Don't hesitate to report bugs or features!

## Cite

It's a great pleasure to know, that your software is useful. Please cite ZGA:

Korzhenkov A. 2021. ZGA: a flexible pipeline for read processing, de novo assembly and annotation of prokaryotic genomes. bioRxiv https://doi.org/10.1101/2021.04.27.441618

And of course tools it's using:

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560

Bushnell, B., Rood, J., & Singer, E. (2017). BBMerge–accurate paired shotgun read merging via overlap. PloS one, 12(10). https://doi.org/10.1371/journal.pone.0185056

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., ... & Pyshkin, A. V. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology, 19(5), 455-477. https://dx.doi.org/10.1089%2Fcmb.2012.0021

Wick, R. R., Judd, L. M., Gorrie, C. L., & Holt, K. E. (2017). Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS computational biology, 13(6), e1005595. https://doi.org/10.1371/journal.pcbi.1005595

Vaser, R., Sović, I., Nagarajan, N., & Šikić, M. (2017). Fast and accurate de novo genome assembly from long uncorrected reads. Genome research, 27(5), 737-746. https://genome.cshlp.org/content/27/5/737.full

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. https://dx.doi.org/10.1093/bioinformatics/bty191

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature biotechnology, 37(5), 540-546. https://doi.org/10.1038/s41587-019-0072-8

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055. https://dx.doi.org/10.1101%2Fgr.186072.114

Tanizawa, Y., Fujisawa, T., & Nakamura, Y. (2018). DFAST: a flexible prokaryotic genome annotation pipeline for faster genome publication. Bioinformatics, 34(6), 1037-1039. https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtx713

Camacho, C., Coulouris, G., Avagyan, V. et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics 10, 421. https://doi.org/10.1186/1471-2105-10-421

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423. https://doi.org/10.1093/bioinformatics/btp163

O’Connell, J., et al. (2015) NxTrim: optimized trimming of Illumina mate pair reads. Bioinformatics 31(12), 2035-2037. https://doi.org/10.1093/bioinformatics/btv057

Ondov, B.D., Treangen, T.J., Melsted, P. et al. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x
