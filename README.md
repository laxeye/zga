# ZGA - prokaryotic genome assembly and annotation pipeline


## Installation


### Installing dependencies

ZGA is written in Python and tested with Python 3.6 and Python 3.7. ZGA uses several software and libs including:

* fastqc
* ea-utils
* bbmap (or seqprep)
* spades
* unicycler
* CheckM
* DFast
* BioPython

All of them may be installed using **conda**:

It's highly recommended to create a new conda environment:

`conda create -n newcoolenv python=3.7 fastqc ea-utils spades unicycler checkm-genome seqprep dfast bbmap blast biopython`

and activate it

`conda activate newcoolenv`


Otherwise you may install dependencies to existing conda environment:

`conda install python>=3.6 fastqc ea-utils spades unicycler checkm-genome seqprep dfast bbmap blast biopython`


Of course, it's possible to use *another ways* even compile all tools from source code. In this case you should check if binaries are in your '$PATH' variable.


### Get source from Github

You can get ZGA by cloning from the repository with `git clone https://github.com/laxeye/zga.git` or by downloading an archive.


### Install from PyPi

Run `pip install zga` it will check if You have Biopython and istall it if not. But all other dependencies You should install manually or using **conda**. CheckM is available on **PyPi**, but it's easier to install it using **conda**.


### Operating systems requirements

ZGA was tested on Ubuntu 18.04. Most probably any modern 64-bit Linux distribuition is enough.

Your feedback on other OS is welcome!


## Usage

You should run `zga.py` if You cloned the source code or `zga` otherwise.

Run 'zga.py -h' to get a help message.

Examples:

Perform all steps (read qc, read trimming and merging, assembly, CheckM assesment with default (bacterial) marker set, DFAST annotation) and use 4 threads, where possible:

`zga.py -1 R1.fastq.gz -2 R2.fastq.gz --threads 4 -o my_assembly`

or use SPAdes and provide it with paired-end and nanopore reads of archaeal genome (Checfkm will use archaeal markers)

`zga.py -1 R1.fastq.gz -2 R2.fastq.gz --nanopore MiniION.fastq.gz -a spades --threads 4 --domain archaea -o my_assembly`

or from Nanopore reads only using unicycler

`zga.py --nanopore MiniION.fastq.gz -o nanopore_assembly`

Perform genome assesment and annotation:

With 'Pectobacterium' CheckM marker set: 

`zga.py --step check -g pectobacterium_sp.fasta --checkm_rank genus --checkm_taxon Pectobacterium -o my_output_dir`

Let CheckM to infer the right marker set: 

`zga.py --step check -g my_genome.fa --checkm_mode lineage -o my_output_dir`


## Know issues and limitations

Don't forget: ZGA is in the early testing...

I hope to fix next issues **ASAP**:

* It's not posible to provide multiple read libraries i.e. tow sets of PE reads or two nanopore runs. 
* It's not possible to install all dependencies with Python 3.8 via conda, please use 3.7 or 3.6.
* There is no conda package

Don't hesitate to report bug or feature!


## Cite

It's a great pleasure to know, that your software is useful. Please cite ZAG: 

Korzhenkov A. (2020). ZGA: prokaryotic genome assembly and annotation pipeline.

And of course tools it's using:

Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.

Aronesty, E. (2015). ea-utils: Command-line tools for processing biological sequencing data. 2011. URL https://github. com/ExpressionAnalysis/ea-utils.

Bushnell, B., Rood, J., & Singer, E. (2017). BBMerge–accurate paired shotgun read merging via overlap. PloS one, 12(10).

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., ... & Pyshkin, A. V. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology, 19(5), 455-477.

Wick, R. R., Judd, L. M., Gorrie, C. L., & Holt, K. E. (2017). Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS computational biology, 13(6), e1005595.

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055.

Tanizawa, Y., Fujisawa, T., & Nakamura, Y. (2018). DFAST: a flexible prokaryotic genome annotation pipeline for faster genome publication. Bioinformatics, 34(6), 1037-1039.

Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of molecular biology, 215(3), 403-410.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423.

