kir-mapper
=======

Castelli et al., 2024.

Version 1.0 (November, 2024)

Author: Erick C. Castelli (erick.castelli@unesp.br)


## What is kir-mapper?

kir-mapper is a toolkit for calling SNPs, alleles, and haplotypes for KIR genes from short-read second-generation sequencing (NGS) data. kir-mapper supports both single-end and paired-end Illumina sequencing data. It is compatible with Ion Torrent data uppon some adjustments. This toolkit presents methods for:

1. Getting unbiased alignments in the context of the hg38 reference genome.
2. Estimating copy numbers.
3. Calling SNPs and InDels across the KIR genes in the context of the hg38 reference genome.
4. Calling KIR alleles, with reports listing potential new SNPs.
5. Inferring haplotypes within KIR genes, and among KIR genes.

## Summary


[What is kir-mapper](#what-is-kir-mapper)

[Important notes](#important-notes)

[Install](#install)

[-- Installing using Conda/Miniconda](#installing-using-conda)

[-- Installing everything by yourself](#installing-everything-by-yourself)

[kir-mapper configuration](#kir-mapper-configuration)

[Quick reference for kir-mapper usage](#quick-reference-for-kir-mapper-usage)

[-- Aligning reads to the hg38 reference genome - map](#aligning-reads-to-the-hg38-reference-genome---map)

[-- Estimating copy numbers - ncopy](#estimating-copy-numbers---ncopy)

[-- Calling SNPs and alleles - genotype](#calling-snps-and-alleles---genotype)

[-- Calling haplotypes and solving ambiguites - haplotype](#calling-haplotypes-and-solving-ambiguites---haplotype)

[Other methods](#other-methods)

[Practical notes](#practical-notes)

[Support](#support)

[Manual](MANUAL.md)


## Important notes:

Data compatibility: We tested kir-mapper with Illumina short-read data from whole-genome sequencing (WGS), whole-exome sequencing (WES), and targeted sequencing. It might work with Ion Torrent with some adjustments.

System compatibility: MacOS (Intel), Linux, or WSL2/Linux. We have tested it with MacOS 10.15, Ubuntu 22.04 LTS, and Ubuntu 22.04 LTS under WSL2. Other versions might be compatible. For MacOS, we tested only with Intel Macs.

Read depth: Please note that read depth is essential. We recommend coverage of at least 20x for WGS and 50x for WES. 

Read size: You will get much better results when dealing with a read size larger than 100 nucleotides and paired-end sequencing. kir-mapper is also compatible with single-end sequencing data. The pipeline might produced biases results with shorter reads ( < 100).

Sample size: The minimum sample size we tested is 50 samples. The sample size is essential to get accurate estimations for copy numbers. 

Always indicate the full path for any input file or output folder.

[Back to Summary](#summary)


## Install

kir-mapper depends on a list of libraries and third-party programs, including samtools, bcftools, freebayes, and others. In addition, it depends on some libraries such as ZLIB and BOOST. 

You can choose how to install kir-mapper and its dependencies. You can opt for Conda, Docker, or install everything by yourself. We recommend Conda.


### Installing using conda

To install kir-mapper and all its dependencies, use Conda and the kir-mapper.yml file, as follows.

<code style="color : red">This tutorial assumes that you are using Miniconda (version 3). You can adapt it as necessary.</code>

<code style="color : red">USER represents your username. Change it to fit your username.</code>
 


1. If Conda/Bioconda is not installed, follow the [Bioconda installation instructions](https://bioconda.github.io/), with a preference for [Miniconda](https://docs.anaconda.com/free/miniconda/)

2. Make sure you added the proper channels:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority true
```

3. Clone the [kir-mapper GitHub repo](https://github.com/erickcastelli/kir-mapper)
```
git clone https://github.com/erickcastelli/kir-mapper
```

4. Enter the kir-mapper repository.
```
cd kir-mapper
```

5. Now, download the last version of the kir-mapper database and unzip it:
```
wget --no-check-certificate https://www.castelli-lab.net/support/kir-mapper_db_latest.zip
unzip kir-mapper_db_latest.zip
```

6. Download PICARD tools.
```
wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar
```

7. Use conda to create an environment for kir-mapper using the kir-mapper.yml from the repository
```
conda env create -f kir-mapper.yml
```

8. Now, activate the kir-mapper environment.
```
conda activate kir-mapper
```

9. From the kir-mapper directory (you are already there), create a new folder named `build` and enter it.
```
mkdir build && cd build
```

10. Compile kir-mapper from the /build folder. If this doesn't work, try step 11.
```
cmake ../src/
make
```

11. If step 10 doesn't work, try this (assuming that your conda kir-mapper environment is on /home/USER/miniconda/envs/kir-mapper). **Replace USER by your username.**
```
BOOST_ROOT=/home/USER/miniconda3/envs/kir-mapper ZLIB_ROOT=/home/USER/miniconda3/envs/kir-mapper cmake ../src/
make
```

12. Copy the kir-mapper binary to the /usr/bin, or /usr/local/bin, or folder /bin from your kir-mapper environment (e.g.: /home/USER/miniconda3/envs/kir-mapper/bin). Alternativelly, you can run `kir-mapper` from the build folder. **Replace USER by your username.**
```
cp kir-mapper /home/USER/miniconda3/envs/kir-mapper/bin/
```
or
```
cp kir-mapper /usr/local/bin
```


13. Run kir-mapper. The setup process usually starts automatically. If it doesn't, you can call it by typing the following:
```
kir-mapper setup
```

14. Follow the setup steps. kir-mapper will automatically detect most programs (BWA, samtools, bcftools, freebayes, etc). The only exception is the path for the database (from step 5) and for PICARD tools (from step 6). 

15. You are all set. Don't forget to activate the kir-mapper environment before using it.
```
conda activate kir-mapper
```

[Back to Summary](#summary)




<br/><br/>

### Installing everything by yourself

kir-mapper depends on a list of libraries and third-party programs, as follows:

- boost 1.74
- cmake >= 3.26.4
- make >= 4.3
- zlib >= 1.2.13
- R >= 4.2
- R ggplot2, plotly, htmlwidgets, stringr, shiny, dplyr, forcats, pacman
- gcc/cpp compiler >= 11
- bwa 0.7.17
- freebayes 1.3.6
- samtools 1.19.2
- bcftools 1.19
- STAR 2.7.10b or 2.7.11a
- whatshap 2.2
- shapeit4 4.2.2
- picard-tools and java >= 11
- bgzip and tabix


To compile the program by your self, assuming that everything above is available, follow these steps:

1. Clone the [kir-mapper GitHub repo](https://github.com/erickcastelli/kir-mapper)
```
git clone https://github.com/erickcastelli/kir-mapper
```

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

5. Now, download the last version of the kir-mapper database:
```
wget --no-check-certificate https://www.castelli-lab.net/support/kir-mapper_db_latest.zip
```


6. Unzip the database.
```
unzip kir-mapper_db_latest.zip
```

7. Run kir-mapper. The setup process usually starts automatically. If it doesn't, you can call it by typing the following:
```
kir-mapper setup
```

8. Follow the setup steps. kir-mapper will automatically detect most programs if they are available in the system. The only exception is the path for the database (from steps 5 and 6) and for PICARD tools. You can also indicate the path for each binary.


[Back to Summary](#summary)

<br/><br/>

## Kir-mapper configuration

If you followed any installation mode described above, kir-mapper is already configured. If it is not, you can configure it by typing `kir-mapper setup`

Remember, you need a copy of the kir-mapper database to run any analysis. 

```
wget --no-check-certificate https://www.castelli-lab.net/support/kir-mapper_db_latest.zip
unzip kir-mapper_db_latest.zip
```

kir-mapper uses a hidden configuration file (.txt) in your home folder containing the path for all necessary programs. If the program does not find this file, it enters the setup mode automatically. You can also call this mode by typing `kir-mapper setup` 

```
kir-mapper setup
```

Follow the instructions provided to indicate the path of all necessary programs. kir-mapper might find the programs automatically. The only exception is the database and the PICARD jar file. 

The setup process will save the configuration file in your home folder. This is an example of this file. You can edit it by using `nano ~/.kir-mapper`. **Replace USER by your username.**

	db=/home/USER/kir-mapper/kir-mapper_db_latest/
	samtools=/home/USER/miniconda3/envs/kir-mapper/bin/samtools
	bcftools=/home/USER/miniconda3/envs/kir-mapper/bin/bcftools
	bwa=/home/USER/miniconda3/envs/kir-mapper/bin/bwa
	whatshap=/home/USER/miniconda3/envs/kir-mapper/bin/whatshap
	freebayes=/home/USER/miniconda3/envs/kir-mapper/bin/freebayes
	picard=/home/USER/miniconda3/envs/kir-mapper/bin/picard.jar
	star=/home/USER/miniconda3/envs/kir-mapper/bin/STAR
	shapeit4=/home/USER/miniconda3/envs/shapeit4/bin/shapeit4

[Back to Summary](#summary)

<br/><br/>

## Quick reference for kir-mapper usage

For full details on how to use kir-mapper, please check kir-mapper documentation [MANUAL.md](MANUAL.md)

In brief, there are **four main methods**, that should be used in this specific order:
- map
- ncopy
- genotype
- haplotype

Typing `kir-mapper` will display all the functions available.

Typing `kir-mapper map`, for instance, will display all the options for the `map` function.

### Aligning reads to the hg38 reference genome - map

	Usage: kir-mapper map [OPTIONS]
	Required:
		-r1  STRING and -r2 STRING: path to paired-end read files
			 or
		-r0  STRING: path to single-end read file
			 or
		-bam STRING: path to BAM file (reads aligned to the hg38 reference with BWA-MEM)
		
		
	Optional:
		-output STRING: full path to the output folder.
		-sample STRING: name/id for the sample
		-threads INT: number of threads
		--exome: indicate that this data is WES
		(check manual for other options)


Example for a sample tagged as "Test". "Test" will be the name for the sample in all outputs.
```	
# Re-aligning a BAM file
kir-mapper map -bam original_BAM.bam -sample test -output /home/USER/output 

# Aligning FASTQ 
kir-mapper map -r1 R1.fastq.gz -r2 R2.fastq.gz -sample test -output /home/USER/output 

# Aligning FASTQ from exomes 
kir-mapper map -r1 R1.fastq.gz -r2 R2.fastq.gz -sample test -output /home/USER/output --exome 
```

Examples using the sample data provided in /samples
```	
kir-mapper map -r1 HG00096.R1.fastq.gz -r2 HG00096.R2.fastq.gz -sample HG00096 -output /home/USER/output 
kir-mapper map -r1 HG02461.R1.fastq.gz -r2 HG02461.R2.fastq.gz -sample HG02461 -output /home/USER/output 
kir-mapper map -bam HG00403.KIR.bam -sample HG00403 -output /home/USER/output 
kir-mapper map -bam HG01583.KIR.bam -sample HG01583 -output /home/USER/output 
```

When evaluating many samples simultaneously, run `map` for every sample, indicating the same output but different sample names (as indicated in the example above).

The outputs from `map` are BAM files with aligned reads to the hg38 reference genome and gene-specific fastq files. The final BAM is the ".adjusted.bam" when not using PICARD tools or ".adjusted.nodup.bam" when using PICARD tools. You can inspect/explore the BAM files using IGV.




[Back to Summary](#summary)


### Estimating copy numbers - ncopy

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
or
kir-mapper ncopy -output /home/USER/output --exome
```

This function will estimate the number of copies for every KIR gene and sample. The final outputs are plots in PNG and HTML format with the coverage ratio between the target gene and the selected reference. Users must evaluate the plots to define the correct thresholds and edit the thresholds.txt file accordingly. If any threshold is modified, you must run `ncopy` again to reflect the modifications.

To evaluate the thresholds, using a browser, please open the .html files inside folder /home/USER/output/ncopy/plots. Define the thresholds to separate samples with 0, 1, 2, 3, or >3 copies, changing it on the thresholds.txt at /home/USER/output/ncopy

Then, run ncopy again. This will update all the plots and the copy numbers for all samples.
```	
kir-mapper ncopy -output /home/USER/output 
```

Alternatively, you can use the R script named `kir-mapper_plot_app.R` inside the /home/USER/output/ncopy. This script can assist you in defining the thresholds and updating the plots. When using this script, there is no need to run `copy` again in case you change any threshold.

[Back to Summary](#summary)


### Calling SNPs and alleles - genotype

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

The outputs are VCF files for every gene, and reports with the detected alleles for every sample, listing eventual mismatches.

The VCF files are placed inside /home/USER/output/genotype/[GENE_NAME]/vcf
The reports for each sample are placed inside /home/USER/output/genotype/[GENE_NAME]/reports
The summary with all allele calls is placed inside /home/USER/output/genotype/[GENE_NAME]/calls

All the SNPs are reported in the context of the hg38 reference genome. For genes that are not annotated in the primary sequence of chr19 (e.g. KIR2DL5), reads from these  are aligned and reported in an alternative contig.

These are the locations for all genes in the alternative contigs:

- KIR2DL2,	chr19_KI270921v1_alt:53185-67900
- KIR2DL5AB, chr19_KI270921v1_alt:175661-185557
- KIR2DS1, chr19_KI270921v1_alt:204223-218437
- KIR2DS2, chr19_KI270921v1_alt:36890-51500
- KIR2DS3, chr19_KI270921v1_alt:81118-95700
- KIR2DS5, chr19_KI270890v1_alt:36829-52100
- KIR3DP1, chr19_KI270923v1_alt:61981-67693
- KIR3DS1, chr19_KI270921v1_alt:159375-174162



Sometimes, `kir-mapper genotype` reports ambiguities, i.e., more than one combination of alleles that fit the observed genotypes. The following method `kir-mapper haplotype` may solve ambiguities.

[Back to Summary](#summary)


### Calling haplotypes and solving ambiguites - haplotype

	Usage: kir-mapper haplotype [OPTIONS]
	Required:
		-output STRING: full path to the output folder. The same used by function map, ncopy, genotype.
		
	Optional:
		-threads INT:   number of threads
		--centromeric:  calling haplotype only on the centromeric genes
		--telomeric:    calling haplotype only on the telomeric genes
		(check manual for all the optionals)

Example. Attention, this function only works with higher sample sizes.
```	
kir-mapper haplotype -output /home/USER/output --centromeric
```

This function will call full haplotypes (all variants will be phased) using shapeit4. Afterward, it generates the predicted sequences for each gene and sample, comparing them with those in the IPD-IMGT/KIR database.

The outputs are phased VCFs and reports with the detected alleles for every sample.

All the SNPs are reported in the context of the hg38 reference genome. For genes that are not annotated in the primary sequence of chr19 (e.g. KIR2DL5), reads are aligned and reported in an alternative contig. 

[Back to Summary](#summary)

<br/><br/>

## Practical notes

### Custom database, with alleles that are not in the IMGT/HLA database 
For now, it is not possible to add new alleles to kir-mapper. We will update the database regularly. Please contact the author if you need to add something.

### Evaluating samples from different ancestry backgrounds 
We do not recommend applying `kir-mapper ncopy` to samples from different populations. The thresholds are quite different among Europeans, Africans, and Asians, for instance. In our tests, we applied `kir-mapper map` and `kir-mapper copy` to all European or African samples separately. After that, we grouped all samples by using `kir-mapper group` to create a kir-mapper output with all samples before running `kir-mapepr genotype`.


## Support

Create a [GitHub issue](https://github.com/erickcastelli/kir-mapper/issues).

