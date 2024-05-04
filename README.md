kir-mapper
=======

Described in: 

Castelli EC et al. ...

Version 1.0 (May 1st, 2024)

Author: Erick C. Castelli (erick.castelli@unesp.br)


## What is kir-mapper?

Kir-mapper is a toolkit for calling SNPs, alleles, and haplotypes for KIR genes from short-read second-generation sequencing (NGS) data. Kir-mapper supports both single-end and paired-end Illumina sequencing data. This toolkit presents methods for:

1. Getting unbiased alignments in the context of the hg38 reference genome.
2. Estimating copy numbers.
3. Calling SNPs and InDels across the KIR genes in the context of the hg38 reference genome.
4. Calling KIR alleles based on copy numbers, with reports listing potential new SNPs.
5. Inferring haplotypes within and among KIR genes.


## Important notes:

Data compatibility: We tested kir-mapper with whole-genome sequencing (WGS), whole-exome sequencing (WES), and with targeted sequencing data. It was tested with short reads from Illumina (WGS and WES). It might work with Ion with some adjustments.

System compatibility: MacOS, Linux, or WSL2/Linux. We have tested it with MacOS 10.15, Ubuntu 22.04 LTS, and Ubuntu 22.04 LTS under WSL2. Other versions might be compatible.

Read depth: Please note that read depth is essential. We recommend coverage of at least 30x for WGS and 50x for WES. 

Read size: You will get better results when dealing with a read size larger than 75 nucleotides and paired-end sequencing data, although Kir-mapper is also compatible with single-end sequencing data.

Sample size: The minimum sample size we tested is 50 samples. Sample size is important to get accurate estimations for copy numbers. 


## Install


### Option 1 (using conda)

Kir-mapper depends on a list of libraries and third-part programs, including samtools, bcftools, freebayes, and others. In addition, it depends on some libraries such as ZLIB and BOOST. To install kir-mapper and all its dependences in the easiest way possible, please use the .yml files from the kir-mapper repository, as follows.

**This tutorial assumes that you are using Miniconda3. You can adapt it as necessary.**

1. First, follow the [Bioconda installation instructions](https://bioconda.github.io/), with a preference for [Miniconda](https://docs.anaconda.com/free/miniconda/)


2. Make sure you added the proper channels:

```
	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda config --set channel_priority true
```

3. Clone the [Kir-mapper gitHub repo](https://github.com/erickcastelli/kir-mapper), e.g. with `git clone https://github.com/erickcastelli/kir-mapper.git`

4. Enter the kir-mapper repository.
```
	cd kir-mapper
```

5. Now, download the last version of the kir-mapper database [here](www.castelli-lab.net/apps/kir-mapper)


6. Unzip the database.
```
	unzip database_file.zip
```

7. Use conda to create an environment for kir-mapper using the kir-mapper.yml from the repository
```
	conda env create -f src/kir-mapper.yml
```

8. Use conda to create an environment to compile kir-mapper using the kir-mapper_compile.yml from the repository.
```
	conda env create -f src/kir-mapper_compile.yml
```

9. From the kir-mapper directory (you are already there), create a new folder named `build` and enter it.
```
	mkdir build && cd build
```

10. Now, activate the kir-mapper_compile environment.
```
	conda activate kir-mapper_compile
```

11. Compile kir-mapper
```
	cmake ../src
	make
```

12. If step 9 does't work, try this (assuming that your conda kir-mapper_compile environment is on /home/USER/miniconda/envs/kir-mapper_compile). Otherwise, skip this step.
```
	BOOST_ROOT=/home/USER/miniconda/envs/kir-mapper_compile ZLIB_ROOT=/home/USER/miniconda/envs/kir-mapper_compile cmake ..
	make
```

13. Copy the kir-mapper binary to the other kir-mapper environment, to /usr/bin, or /usr/local/bin
```
	cp kir-mapper /home/USER/miniconda3/envs/kir-mapper/bin/
```

14. Deactivate the kir-mapper_compile environment:
```
	conda deactivate
```

15. Now, activate the kir-mapper environment.
```
	conda activate kir-mapper
```

16. Run kir-mapper. The setup process usually starts automatically. If it doesn't, you can call it by typing the following:
```
	kir-mapper setup
```

17. Follow the setup steps. Kir-mapper will automatically detect most programs (BWA, samtools, bcftools, freebayes, etc). The only exception is the path for the database (from step 5 and 6) and for PICARD tools. 

18. You are all set. Don't forget to activate the kir-mapper environment before using it.
```
	conda activate kir-mapper
```



### Option 2 (installing everything yourself)

Kir-mapper depends on a list of libraries and third-party programs, as follows:

- boost 1.70 or 1.74
- cmake >= 3.26.4
- make >= 4.3
- zlib >= 1.2.13
- R >= 4.2
- R ggplot2, plotly
- gcc/cpp compiler >= 11
- bwa 0.7.17
- freebayes 1.3.6
- samtools >= 1.18
- bcftools 1.13, or >= 1.18
- STAR 2.7.10b or 2.7.11a
- whatshap 2.2
- shapeit4
- picard-tools and java >= 11
- bgzip and tabix


To compile the program by yourself, assuming that everything above is available, follow these steps:

1. Clone the [Kir-mapper gitHub repo](https://github.com/erickcastelli/kir-mapper), e.g. with `git clone https://github.com/erickcastelli/kir-mapper.git`

2. Enter the kir-mapper repository
```
	cd kir-mapper
```

3. Create a new folder named `build` and enter it
```
	mkdir build && cd build
```

4. Compile the program
```
	cmake ../src/
	make
```

5. Now, download the last version of the kir-mapper database [here](www.castelli-lab.net/apps/kir-mapper)


6. Unzip the database.
```
	unzip database_file.zip
```

7. Run kir-mapper. The setup process usually starts automatically. If it doesn't, you can call it by typing the following:
```
	kir-mapper setup
```

Suppose you want to run kir-mapper without specifying the directory. In that case, you can either add the directory of kir-mapper to the environment variable PATH or create a soft link ("ln -s") of the file "kir-mapper" to a directory in PATH.



## Kir-mapper configuration

Before using kir-mapper, you need to download the lastest version of its database 
 [here](www.castelli-lab.net/apps/kir-mapper). Ignore this step if you already done that during the instalation process.

Once you downloaded it, unzip the database.

```
	unzip database_file.zip
```

If you installed kir-mapper dependences using conda, activate the kir-mapper environment before continue
```
	conda activate kir-mapper
```

Kir-mapper uses a hidden configuration file (.txt) in your home folder containing the path for all necessary programs. If the program does not find this file, it enters the setup mode automatically. You can also call this mode by typing `kir-mapper setup` 

```
	kir-mapper setup
```

Follow the instructions provided to indicate the path of all necessary programs. If the conda kir-mapper environment is activated, kir-mapper will automatically detect most programs. The only exception is the database and the PICARD jar file. 

The setup process will save the configuration file in your home folder. This is an exemple of this file. You can edit it by using `nano ~/.kir-mapper`:

	db=/home/USER/kir-mapper/kir-mapper_db_001.0_KIR/
	samtools=/home/USER/miniconda3/envs/kir-mapper/bin/samtools
	bcftools=/home/USER/miniconda3/envs/kir-mapper/bin/bcftools
	bwa=/home/USER/miniconda3/envs/kir-mapper/bin/bwa
	whatshap=/home/USER/miniconda3/envs/kir-mapper/bin/whatshap
	freebayes=/home/USER/miniconda3/envs/kir-mapper/bin/freebayes
	picard=/home/USER/miniconda3/envs/kir-mapper/bin/picard.jar
	star=/home/USER/miniconda3/envs/kir-mapper/bin/STAR
	shapeit4=/home/USER/miniconda3/envs/shapeit4/bin/shapeit4


## Quick reference for kir-mapper usage

Kir-mapper is a toolkit with methods for aligning, genotyping, and also to infer haplotypes. Please check the full manual for details.

In brief, there are **four main methods**, that should be used in this specific order:
- map
- ncopy
- genotype
- haplotype

When using conda, don't forget to activate the environment `conda activate kir-mapper`

### Aligning reads to the reference genome (map)

	Usage: kir-mapper map [OPTIONS]
	Required:
		-r1  STRING and -r2 STRING: path to paired-end read files
			 or
		-r0  STRING: path to single-end read file
			 or
		-bam STRING: path to BAM file (reads aligned to the hg38 reference with BWA-MEM)
		-sample STRING: name/id for the sample
		-output STRING: full path to the output folder.
		
		
	Optional:
		-threads INT: number of threads
		--exome: indicate that this data is WES
		(check manual for all the optionals)


Example
```	
	# Re-aligning a BAM file
	kir-mapper map -bam original_BAM.bam -sample test -output /home/USER/output 
	# Aligning FASTQ 
	kir-mapper map -r1 R1.fastq.gz -r2 R2.fastq.gz -sample test -output /home/USER/output 
```

When evaluating many samples simultaneosly, run `map` for every sample indicating the same output for all of them.

The outputs from `map` are: BAM with aligned reads to the hg38 reference genome and gene-specific fastq files.


### Estimating copy numbers (function ncopy)

	Usage: kir-mapper ncopy [OPTIONS]
	Required:
		-output STRING: full path to the output folder. The same used by function map.
		
	Optional:
		-threads INT: number of threads
		-reference STRING: the reference with two copies for all samples - KIR3DL3, HLA-G, HLA-E, 5UPKIR
		--exome: indicate that this data is WES
		(check manual for all the optionals)


Example
```	
	kir-mapper ncopy -output /home/USER/output 
```

This function will estimate the number of copies of every KIR gene in all samples. The final outputs are plots (PNG and HTML) with the coverage ratio between the target gene and the reference. Users must evaluate the plots to define the correct thresholds and edit the thresholds.txt file acordingly and run ncopy again to reflect the modifications.




### Calling SNPs, InDels, and alleles (function genotype)

	Usage: kir-mapper genotype [OPTIONS]
	Required:
		-output STRING: full path to the output folder. The same used by function map and ncopy.
		
	Optional:
		-threads INT: number of threads
		--full: call SNPs and Indels also in introns
		(check manual for all the optionals)

Example
```	
	kir-mapper genotype -output /home/USER/output 
```

This function will call SNPs and InDels across all exons from KIR genes, by using freebayes and an internal algorithm to detect and remove unlike genotypes. It also phases the variants using whatshap. After, the program detects which KIR alleles are compatible with the observed variants.

The outputs are VCF files for every gene, and reports with the detected alleles for every samples, listing eventual mismatches.

All the SNPs are reported in the context of the hg38 reference genome. For genes that are not annotated in the primary sequence of chr19 (e.g. KIR2DL5), reads from these  are aligned and reported in an alternative contig.

Sometimes, `kir-mapper genotype` reports ambiguities, i.e., more than one combination of alleles that fit the observed genotypes. The following method may solves ambiguities.


### Calling haplotypes and solving ambuguites (function haplotype)

	Usage: kir-mapper haplotype [OPTIONS]
	Required:
		-output STRING: full path to the output folder. The same used by function map, ncopy, genotype.
		
	Optional:
		-threads INT:   number of threads
		--centromeric:  calling haplotype only on the centromeric genes
		--telomeric:    calling haplotype only on the telomeric genes
		(check manual for all the optionals)

Example
```	
	kir-mapper haplotype -output /home/USER/output --centromeric
```

This function will call full haplotypes (all variants will be phased) by using shapeit4. After, it gererate the predicted sequences for each gene and sample, comparing the sequences with the ones at the IPD-IMGT/KIR database.

The outputs are phased VCFs and reports with the detected alleles for every samples.

All the SNPs are reported in the context of the hg38 reference genome. For genes that are not annotated in the primary sequence of chr19 (e.g. KIR2DL5), reads are aligned and reported in an alternative contig. Please refer to the manual regarding this.



## Other methods 

### Function group

Combines multiple map and ncopy runs in a single output structure.

### Function join

Join variants from all genes in a single VCF file.

### Function select

To extract reads related to KIR genes from FASTQ or BAM files. However, not all reads are KIR. You should use the map function to get gene-specific reads.

### Function setup

Configure (or re-configure) kir-mapper.



## Practical notes

### Custom database, with alleles that are not in the IMGT/HLA database 
For now, it is not possible to add new alleles to kir-mapper. We will update the database regularlly. Please contact the author in case you need to add something.

### Evaluating samples from different ancestry backgrounds 
We do not recommend aplying `kir-mapper ncopy` in samples from different populations. The thresholds are quite different among European, Africans, Asians, and others. In our tests, we applied `kir-mapper map` and `kir-mapper ncopy` in all samples from Europe, or from Africa, etc, separatelly. After, we grouped all samples by using `kir-mapper group` to create a kir-mapper output with all samples before running `kir-mapepr genotype`.


## Support

Create a [GitHub issue](https://github.com/erickcastelli/kir-mapper/issues).

