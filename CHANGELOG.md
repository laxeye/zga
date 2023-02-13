### ZGA changelog

**0.1.0**

+ Short reads normalization by k-mer count (BBNorm)
+ Processed reads are compressed with pigz (if installed) or gzip, while intermediate still in plain FASTQ
+ Providing custom configuration to DFAST
+ Dockerfile added. Warning: the image takes about 8 GB :(
* Merging of paired-end reads (BBmerge) is now optional
+ BBDuk extra options avalaible with--bbduk-extra
+ BBMerge extra options avalaible with--bbmerge-extra
- Arguments --bbmerge-extend, --bbmerge-extend-kmer and --bbmerge-trim for BBMerge were removed.
+ MEGAHIT assembler was added

**0.0.9**

+ Repairing of paired-end read files with different read count (BBmap)

**0.0.8**

+ Genome assembly stats
+ Read correction with tadpole.sh (from BBmap)
+ Multiple libraries of the same type supported for SPAdes and Flye

**0.0.7**

* FastQC replaced with fastp
* fastq-mcf (from ea-utils) replaced with bbduk.sh (from BBmap)
+ conda package available

**0.0.6**

+ Flye assembler implemented
+ Polishing with short reads using racon

**0.0.5**

+ Deleting provisional read files
- SeqPrep dropped

**0.0.4**

+ Genome size estimation with Mash

**0.0.3**

+ Illumina mate-pair read processing with NxTrim

**0.0.2**

+ Replicon extraction from unicycler assemblies
+ Integrated filterbytile.sh from BBMap

**0.0.1**

* Initial release
