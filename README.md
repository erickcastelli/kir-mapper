kir-mapper
=======

Version 1.1 (May, 2026), using IPD-KIR version: 2.15

Author: Erick C. Castelli (erick.castelli@unesp.br)

***To use this version of kir-mapper, please download the database again. The old database is not compatible with this version, and the alleles were updated to IPD-KIR version 2.15. Do not forget to run `kir-mapper setup` again, or update the database path in the configuration file.***

Castelli EC et al. kir-mapper: A Toolkit for Killer-Cell Immunoglobulin-Like Receptor (KIR) Genotyping From Short-Read Second-Generation Sequencing Data. HLA 2025 Mar;105(3):e70092. doi: 10.1111/tan.70092.



## What is kir-mapper?

kir-mapper is a toolkit for analyzing KIR genes from short-read second-generation sequencing (NGS) data. The toolkit provides functions for read alignment, variant calling for the hg38 reference genome, phasing, KIR allele calling, and KIR haplotype estimation. kir-mapper supports both single-end and paired-end Illumina sequencing data, in FASTQ or BAM format. It is compatible with Ion Torrent and Nanopore data upon some adjustments. This toolkit presents methods for:

1. Getting unbiased alignments in the context of the hg38 reference genome.
2. Estimating gene copy numbers.
3. Calling SNPs and InDels across the KIR genes in the context of the hg38 reference genome, phasing them.
4. Calling KIR alleles, with reports listing potential new SNPs.
5. Inferring haplotypes within KIR genes, and among KIR genes.

kir-mapper was written in C++, and it is compatible with Linux. It may work on MacOS, but we haven't tested it. kir-mapper is a single program with multiple functions and provides an interface similar to that of samtools, BWA, or Bowtie2.

## Summary


[What is kir-mapper](#what-is-kir-mapper)

[Update history](#update-history)

[Important notes](#important-notes)

[Install](#install)

[-- Installing using Conda/Miniconda](#installing-using-conda)

[-- Installing using Docker](#installing-using-docker)

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

Data compatibility: We tested kir-mapper with Illumina short-read data from whole-genome sequencing (WGS), whole-exome sequencing (WES), and targeted sequencing. We also tested it with Oxford Nanopore whole-genomes (R10.4.1 only). It might work with Ion Torrent with some adjustments.

System compatibility: Linux, or WSL2/Linux. We have tested it with Ubuntu 22.04 LTS, and Ubuntu 22.04 LTS under WSL2. Other versions might be compatible. It may work on MacOS, but we haven't tested it.

Read depth: Please note that read depth is essential. We recommend coverage of at least 30x for WGS and 50x for WES. 

Read size: You will get much better results when dealing with a read size larger than 100 nucleotides and paired-end sequencing. kir-mapper is also compatible with single-end sequencing data. The pipeline may produce biased results with shorter reads ( < 100).

Sample size: The minimum sample size we tested for the ncopy and haplotype functions is 50. The sample size is essential for obtaining accurate estimates of copy numbers and haplotypes. 

Always indicate the full path for any input file or output folder. For example, do not use "~" for your home folder.


[Back to Summary](#summary)




## Update history:

The original release was version 1.0.

### Version 1.01:
- Updates to the filter to introduce missing alleles in regions with low read depth.

### Version 1.1:
- Updates on the descriptions of functions and options.
- Implementation of a downsample system for very high-coverage samples.
- Implementation of a way to skip marking duplicates with Picard.
- Updates in the database data and structure. The database from versions 1.0 and 1.01 is not compatible with this new kir-mapper version.
- Calling genotypes only in exonic regions is now much faster than the previous version.
- Increased sensitivity to SNPs at the edge of exons when genotyping only exonic regions, which is essential to reduce ambiguities and missing alleles.
- Support for Oxford Nanopore R10.4.1 data.
- The manual was updated, with a detailed description of all options and flags.

[Back to Summary](#summary)


## Install

kir-mapper depends on a list of libraries and third-party programs, including samtools, bcftools, freebayes, and others. In addition, it depends on some libraries such as ZLIB and BOOST. 

You can choose how to install kir-mapper and its dependencies. We strongly recommend installing it in a conda environment or using Docker.

### Installing using conda

To install kir-mapper and all its dependencies, use Conda and the kir-mapper.yml file, as follows.

**This tutorial assumes that you are using Miniconda (version 3). You can adapt it as necessary.**

**This tutorial assumes your Miniconda installation is at /home/USER/miniconda3. You can adapt it as necessary.**

**USER is your username. Change it to fit your user name.**


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

5. Now, download the latest version of the kir-mapper database:
```
	wget --no-check-certificate https://www.castelli-lab.net/support/kir-mapper_db_latest.zip
```

6. Unzip the database.
```
	unzip kir-mapper_db_latest.zip
```

7. Download PICARD tools.
```
	wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar
```

8. Use conda to create an environment for shapeit4 using the shapeit4.yml from the repository
```
	conda env create -f shapeit4.yml
```

9. Use conda to create an environment for kir-mapper using the kir-mapper.yml from the repository
```
	conda env create -f kir-mapper.yml
```

10. Now, activate the kir-mapper environment.
```
	conda activate kir-mapper
```

11. From the kir-mapper directory (you are already there), create a new folder named `build` and enter it.
```
	mkdir build && cd build
```

12. Compile kir-mapper from the /build folder
```
	cmake ../src/
	make
```

13. If step 12 doesn't work, try this (assuming that your conda kir-mapper environment is on /home/USER/miniconda/envs/kir-mapper).
```
	BOOST_ROOT=/home/USER/miniconda3/envs/kir-mapper ZLIB_ROOT=/home/USER/miniconda3/envs/kir-mapper cmake ../src/
	make
```

14. Copy the kir-mapper binary to the /usr/bin, or /usr/local/bin, or folder /bin from your kir-mapper environment (e.g.: /home/USER/miniconda3/envs/kir-mapper/bin). Alternatively, you can run `kir-mapper` from the build folder. 
```
	cp kir-mapper /home/USER/miniconda3/envs/kir-mapper/bin/
	    or
	cp kir-mapper /usr/local/bin
```

15. Copy the shapeit4 binary from the shapeit4 environment to the kir-mapper environment. 
```
	cp /home/USER/miniconda3/envs/shapeit4/bin/shapeit4 /home/USER/miniconda3/envs/kir-mapper/bin/
```

16. Run kir-mapper. The setup process usually starts automatically. If it doesn't, you can call it by typing the following:
```
	kir-mapper
	or
	kir-mapper setup
```

17. Follow the setup steps. kir-mapper will automatically detect most programs (BWA, samtools, bcftools, freebayes, etc). The only exceptions are the paths for the database (from steps 5 and 6) and for Picard (from step 7).  

18. You are all set. Don't forget to activate the kir-mapper environment before using it.
```
	conda activate kir-mapper
```

[Back to Summary](#summary)



### Installing using Docker

Here, we assume that Docker is already installed in your system.

1. Clone the [Kir-mapper GitHub repo](https://github.com/erickcastelli/kir-mapper)
```
git clone https://github.com/erickcastelli/kir-mapper
```


2. Enter the kir-mapper repository.
```
cd kir-mapper
```

3. Build the kir-mapper image. This might take a while.
```
docker build -t kir-mapper -f Dockerfile .
or
sudo docker build -t kir-mapper -f Dockerfile .
```

4. Once completed, you may run the kir-mapper docker image like this:
```
docker run -it kir-mapper
or
sudo docker run -it kir-mapper
```

5. Now, type `kir-mapper`. The program should be available, with all dependencies and the database already set.
```
kir-mapper
```

6. By using this method, Docker has already downloaded the database, and the software is configured. There is no need to run `kir-mapper setup`


7. To run Docker with access to the host system, you may use a command like this: `/home/USER` is the path to the folder where the data you want to process is, and `/data` is how you can access it inside the docker image.
```
docker run -it -v /home/USER:/data kir-mapper
```
[Back to Summary](#summary)



### Installing everything by yourself

kir-mapper depends on a list of libraries and third-party programs, as follows:

- boost 1.74
- cmake >= 3.29
- make >= 4.3
- zlib >= 1.2.13
- R >= 4.2
- R ggplot2, plotly, htmlwidgets, stringr, shiny, dplyr, forcats, pacman
- gcc/cpp compiler >= 11
- bwa 0.7.17
- freebayes 1.3.8
- samtools 1.21
- bcftools 1.20
- STAR 2.7.10b or 2.7.11a
- whatshap 2.2
- shapeit4 4.2.2
- picard-tools and java >= 11
- bgzip and tabix


To compile the program by yourself, assuming that everything above is available, follow these steps:

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

5. Now, download the latest version of the kir-mapper database:
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

8. Follow the setup steps. kir-mapper will automatically detect most programs if they are available in the system. The only exceptions are the database path (from steps 5 and 6) and the Picard tools path. You can also indicate the path for each binary.


[Back to Summary](#summary)



## Kir-mapper configuration

If you followed any of the installation modes described above, kir-mapper is already configured. If it is not, you can configure it by typing `kir-mapper setup`

Remember, you need a copy of the kir-mapper database to run any analysis. 

```
	wget --no-check-certificate https://www.castelli-lab.net/support/kir-mapper_db_latest.zip
	unzip kir-mapper_db_latest.zip
```

kir-mapper uses a hidden configuration file (.txt) in your home folder that contains the paths to all necessary programs. If the program cannot find this file, it automatically enters setup mode. You can also call this mode by typing `kir-mapper setup` 

```
	kir-mapper setup
```

Follow the instructions provided to indicate the path of all necessary programs. kir-mapper might find the programs automatically. The only exception is the database and the PICARD jar file. 

The setup process will save the configuration file in your home folder. This is an example of this file. You can edit it by using `nano ~/.kir-mapper`:

	db=/home/USER/kir-mapper/kir-mapper_db_latest/
	samtools=/home/USER/miniconda3/envs/kir-mapper/bin/samtools
	bcftools=/home/USER/miniconda3/envs/kir-mapper/bin/bcftools
	bwa=/home/USER/miniconda3/envs/kir-mapper/bin/bwa
	minimap=/home/USER/miniconda3/envs/kir-mapper/bin/minimap2
	whatshap=/home/USER/miniconda3/envs/kir-mapper/bin/whatshap
	freebayes=/home/USER/miniconda3/envs/kir-mapper/bin/freebayes
	picard=/home/USER/miniconda3/envs/kir-mapper/bin/picard.jar
	star=/home/USER/miniconda3/envs/kir-mapper/bin/STAR
	shapeit4=/home/USER/miniconda3/envs/shapeit4/bin/shapeit4

[Back to Summary](#summary)



## Quick reference for kir-mapper usage

For full details on how to use kir-mapper, please check kir-mapper documentation [MANUAL.md](MANUAL.md)

In brief, there are four main methods that should be used in this specific order:
- map
- ncopy
- genotype
- haplotype

Typing `kir-mapper` will display all the functions available.

Typing `kir-mapper map`, for instance, will display all the options for the `map` function.

### Aligning reads to the hg38 reference genome - map

This function aligns or realigns reads to KIR genes. Type `kir-mapper map` to check all options.

For full details on how to use kir-mapper, please check kir-mapper documentation [MANUAL.md](MANUAL.md)

***kir-mapper supports both BAM files and FASTQ files as input. However, it is much faster if the input file is a BAM. If your input is FASTQ files, the best practice is to align the reads to the hg38 reference genome using BWA-MEM, producing a sorted BAM file and its BAI index. Use the same hg38 reference used by the 1000 Genomes Project. Then, use this BAM file as input to kir-mapper. Attention: Always provide the full paths for the input files and the output folder.***


Examples on how to run `kir-mapper map` for a sample tagged as "Test". "Test" will be the sample name in all outputs.
```	
# Re-aligning a BAM file
kir-mapper map -bam original_BAM.bam -sample test -output /home/USER/output 

# Re-aligning ONT BAM file
kir-mapper map -bam original_BAM.bam -sample test -output /home/USER/output --nanopore 

# Aligning FASTQ 
kir-mapper map -r1 R1.fastq.gz -r2 R2.fastq.gz -sample test -output /home/USER/output 

# Aligning FASTQ from exomes 
kir-mapper map -r1 R1.fastq.gz -r2 R2.fastq.gz -sample test -output /home/USER/output --exome

# Re-aligning a BAM file from an exome 
kir-mapper map -bam original_BAM.bam -sample test -output /home/USER/output --exome
```

When evaluating many samples simultaneously, run `kir-mapper map` for each sample, indicating a different sample name (-sample) and the same output folder for all of them. You can run multiple instances of kir-mapper with the same output folder, but the sample names must be different..

The outputs from `map` are BAM files with aligned reads to the hg38 reference genome and gene-specific fastq files. The final BAM is the ".adjusted.bam" when not using Picard or ".adjusted.nodup.bam" when using Picard. 

You may inspect the BAM files using [IGV](https://igv.org/). In IGV, change the genome for "Human (hg38 1kg/GATK)" and open the ".adjusted.bam" or ".adjusted.nodup.bam" file. For KIR genes annotated at the main chr19 sequence (e.g., KIR3DL3), type the gene name to locate it. For genes that are not annotated at the main chr19 chromosome, their locations in alternative contigs are as follows:

- KIR2DL2,	chr19_KI270921v1_alt:53185-67900
- KIR2DL5AB, chr19_KI270921v1_alt:175661-185557
- KIR2DS1, chr19_KI270921v1_alt:204223-218437
- KIR2DS2, chr19_KI270921v1_alt:36890-51500
- KIR2DS3, chr19_KI270921v1_alt:81118-95700
- KIR2DS5, chr19_KI270890v1_alt:36829-52100
- KIR3DP1, chr19_KI270923v1_alt:61981-67693
- KIR3DS1, chr19_KI270921v1_alt:159375-174162


WARNING FOR ONT DATA:

When dealing with Oxford Nanopore data, kir-mapper extracts the fragment corresponding to a KIR gene from large reads and renames it. Therefore, the final BAM alignment does not represent the original FASTQ data. The read size is maintained as large as possible to facilitate the definition of internal haplotypes. Still, they do not contribute to the definition of haplotypes between genes.


***When visualizing kir-mapper BAM files in IGV, keep in mind that two classes of reads are marked to be hidden from genotyping tools. Reads mapping to more than one location are flagged as secondary alignments. Picard marks PCR/optical duplicates. Both classes are excluded from genotyping tools such as FreeBayes, even though they are still part of the true alignment. Because IGV does not hide these reads by default, the displayed coverage will appear higher than the actual sequencing depth considered in downstream analysis. Therefore, to give the user a true sense of the alignment and how the genotyping tools will handle it, turn off secondary and duplicated reads in IGV.***



***After that, if you are experiencing very low read depth with high-coverage sequencing data, try turning off mark duplicates with Picard (--skip-markdup). Sequencing libraries enriched by PCR are usually not compatible with Picard's MarkDuplicates.***



[Back to Summary](#summary)


### Estimating copy numbers - ncopy

This function estimates the number of copies for each KIR gene and sample. The `map` function must be applied before `ncopy` to each sample. 

For full details on how to use kir-mapper, please check kir-mapper documentation [MANUAL.md](MANUAL.md)

**Attention: kir-mapper tries to estimate copy numbers for LILR genes (LILRB1, LILRB2, etc). However, this is still in beta mode and the results should not be considered.**

Examples of how to run `kir-mapper ncopy`
```	
	kir-mapper ncopy -output /home/USER/output 
	   or
	kir-mapper ncopy -output /home/USER/output --exome
	   or
	kir-mapper ncopy -output /home/USER/output --nanopore
```

This function estimates the number of copies for each KIR gene and sample. The final outputs are plots in PNG and HTML format with the coverage ratio between the target gene and the selected reference. Users must evaluate the plots to determine the correct thresholds and edit thresholds.txt accordingly. If any threshold is modified, you must rerun `ncopy` to reflect the changes.

To evaluate the thresholds, please open the .html files in/home/USER/output/ncopy/plots in a browser. Define the thresholds to separate samples with 0, 1, 2, 3, or >3 copies, changing it in the thresholds.txt at /home/USER/output/ncopy

Then, run ncopy again. This will update all the plots and the copy numbers for all samples.

Alternatively, you can use the R script named `kir-mapper_plot_app.R` inside the /home/USER/output/ncopy. This script can help you define thresholds and update the plots. When using this script, there is no need to run `ncopy` again in case you change any threshold.

[Back to Summary](#summary)


### Calling SNPs and alleles - genotype

This function will call SNPs and InDels across KIR genes, phase them using WhatsHap, and evaluate how these SNPs and haplotypes align with known KIR alleles. 

For full details on how to use kir-mapper, please check kir-mapper documentation [MANUAL.md](MANUAL.md)

Note that you must indicate the same output folder used in the previous step (map and ncopy).


Examples of how to run `kir-mapper ncopy`
```	
	kir-mapper genotype -output /home/USER/output 
	   or
	kir-mapper genotype -output /home/USER/output --exome
       or
	kir-mapper genotype -output /home/USER/output --nanopore
```


The outputs from `genotype` are placed in a " genotype " folder inside the output folder. By default, variants are called only in exons and placed under the folder " genotype/cds ". When using `--full`, variants are placed under the folder " genotype/full ".

The outputs are VCF files for each gene and reports listing the detected alleles for each sample, including any mismatches.

The VCF files are placed inside output/genotype/cds/[GENE_NAME]/vcf

The reports for each sample are placed inside output/genotype/cds/[GENE_NAME]/reports

The summary with all allele calls is placed inside output/genotype/cds/[GENE_NAME]/calls

All the SNPs are reported in the context of the hg38 reference genome. For genes not annotated in the primary sequence of chr19 (e.g., KIR2DL5), reads from these genes are aligned and reported in an alternative contig.

These are the locations for all genes in the alternative contigs:

- KIR2DL2,	chr19_KI270921v1_alt:53185-67900
- KIR2DL5AB, chr19_KI270921v1_alt:175661-185557
- KIR2DS1, chr19_KI270921v1_alt:204223-218437
- KIR2DS2, chr19_KI270921v1_alt:36890-51500
- KIR2DS3, chr19_KI270921v1_alt:81118-95700
- KIR2DS5, chr19_KI270890v1_alt:36829-52100
- KIR3DP1, chr19_KI270923v1_alt:61981-67693
- KIR3DS1, chr19_KI270921v1_alt:159375-174162


Sometimes, `kir-mapper genotype` reports ambiguities, i.e., more than one combination of alleles that fit the observed genotypes. This was observed above for sample Test2. The following method, `kir-mapper haplotype,` may solve ambiguities. If there are too many allele combinations, kir-mapper will indicate "*unresolved".

***Tip: You can inspect the kir-mapper VCF and BAM files simultaneously in IGV by setting hg38 as the reference.***


[Back to Summary](#summary)


### Calling haplotypes and solving ambiguites - haplotype

This method will use shapeit4 to call haplotypes within KIR genes and among KIR genes (all SNPs and InDels will be phased). This function will not work properly with fewer than 100 samples, and it is not available for fewer than 20 samples. It generates predicted sequences for each gene and sample, then compares them with those in the IPD-KIR database.

For full details on how to use kir-mapper, please check kir-mapper documentation [MANUAL.md](MANUAL.md)

Example
```	
	kir-mapper haplotype -output /home/USER/output
	   or
	kir-mapper haplotype -output /home/USER/output --centromeric
	   or
	kir-mapper haplotype -output /home/USEr/output -target KIR2DL1,KIR2DL2
```

The outputs are phased VCFs and reports with the detected alleles for each sample and how these alleles are arranged across chr19.

**Attention: the phased VCF uses dummy/fake positions, keeping the expected order of each gene and SNP. Do not use this VCF or consider these positions.**

After running `kir-mapper haplotype`, the user must compare the results from the `haplotype` function (h1 and h2) with those from the `genotype` function (call, ratio, miss). Usually, h1 and h2 will indicate alleles also identified by the `genotype` function.


[Back to Summary](#summary)



## Practical notes

### Custom database, with alleles that are not in the IPD-KIR database 
For now, it is not possible to add new alleles to kir-mapper. We will update the database regularly. Please get in touch with the author if you need to add something.

### Evaluating samples from different ancestry backgrounds 
We do not recommend applying `kir-mapper ncopy` to samples from different populations simultaneously. The thresholds vary widely across populations. In our tests, we applied `kir-mapper map` and `kir-mapper ncopy` in all samples from Europe, Africa, etc., separately. Afterward, we grouped all samples using `kir-mapper group` to create a kir-mapper output with all samples before running `kir-mapper genotype`.

### Always check a few BAM files before processing large datasets
Always analyze a few samples and check the BAM files in IGV before continuing with large datasets. To do that, turn off secondary and duplicate reads on IGV. You might need to adjust the map function, such as turning off Picard markduplicates (--skip-markdup).


## Support

Create a [GitHub issue](https://github.com/erickcastelli/kir-mapper/issues).

