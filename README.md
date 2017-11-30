# MANUAL deconvoNorm

2017-03-14

Marianne S. Felix

marianne.sabourin-felix@hotmail.com
-----------------------------------

marianne.sabourin-felix.1@ulaval.ca
-----------------------------------

Paper : A Deconvolution Protocol for ChIP-Seq Reveals Analogous Enhancer Structures on the Mouse and Human Ribosomal RNA Genes

PubMed : https://www.ncbi.nlm.nih.gov/pubmed/29158335

  1. INTRODUCTION
  2. WARNINGS
  3. SOFTWARE DEPENDENCIES
  4. FILES REQUIREMENTS
  5. HOW TO USE IT
      1. REQUIRED ARGUMENTS
      2. OPTIONAL ARGUMENTS
      3. EXAMPLES
  6. PRECONDITIONS
  7. OUTPUT FILES
  8. HOW IT WORKS

## INTRODUCTION

Deconvolution of sequencing data (ChIP-seq or DNase-seq) in ultra-deep sequencing context.

## WARNINGS

This script was only tested on Ubuntu 14.04 LTS and Ubuntu 16.04 LTS with Python 2.7.10 and 2.7.11.
You can redistribute and/or modify it under the terms of the GNU General Public License v3 as published by the Free Software Foundation. This software is provided "as is" without warranty of any kind, either express or implied, regarding the software, its merchantability, or its fitness for any particular purpose. See the GNU General Public License for more details.

## SOFTWARE DEPENDENCIES

* Samtools
* Bedtools

## FILES REQUIREMENTS

The input files must be in BAM format.

## HOW TO USE IT

For help type :
```
python deconvoNorm.py -h
```
or
```
python deconvoNorm.py --help
```

###### REQUIRED ARGUMENTS

```
-i INPUT, --input INPUT                 Input DNA file (BAM format)
-f IPFILE, --ipfile IPFILE              Immunoprecipitation file (BAM format) (mutually exclusive with -d)
-d IPDIRECTORY, --ipdirectory INPUT     Directory with many immunoprecipitation files (mutually exclusive with -f)
-c CHRNAME, --chrname CHRNAME           Chromosome of interest (to normalize on)
-o OUTPUT, --outputname OUTPUT          Output folder name (if folder already exists, it will be overwritten)
```

###### OPTIONAL ARGUMENTS

```
-l FRAGMENTLENGTH, --fragmentlength FRAGMENTLENGTH      Sequenced fragment length (Default = 100)
-w WINDOWSIZE, --windowsize WINDOWSIZE                  Smoothing window size (Default = 25)
-t THRESHOLD, --threshold THRESHOLD                     Coverage threshold (Default = 10)
-k, --keepfiles                                         Keep intermediates files (Default = False)
-r, --norpm                                             Don't adjust in Read Per Million (Default = False)
```

FRAGMENTLENGTH : The read will be extended toward their 3' end until reaching the FRAGMENTLENGTH.

WINDOWSIZE : The coverage will be averaged base by base with a sliding window of WINDOWSIZE bp.

THRESHOLD : If there is less than THRESHOLD reads, the coverage at that position will be replaced by zero.

###### EXAMPLES

* Example 1 (One file) :

```
python deconvoNorm.py -i input.bam -f ip.bam -c MmrDNA -l 150 -w 75 -t 50 -o outputFolderName1
```
This command will normalize ip.bam with a FRAGMENTLENGTH of 150, a WINDOWSIZE of 75 and a THRESHOLD of 50.

* Example 2 (Many files) :

```
python deconvoNorm.py -i input.bam -d ipDirectory -c HsrDNA -o outputFolderName2 -k -r
```
This command will normalize the ip.bam in the ipDirectory with the defaults parameters (-l 100, -w 25 -t 10), will keep the intermediates files and won't correct with the RPM ratio.

## PRECONDITIONS

* The input files and directory must exist
* The input file must be in BAM format
* The input directory must contain BAM file(s)
* The FRAGMENTLENGTH (-l) must be an integer
* The WINDOWSIZE (-w) must be an odd number
* The THRESHOLD (-t) must be a float (e.g. 21.5, 10.00, 50)

## OUTPUT FILES

###### FINAL FILE

```
ip-chrName-lFRAGMENTLENGTH-cov-wWINDOWSIZE-tTHRESHOLD-rpm_norm.bedgraph
```

###### INTERMEDIATES FILES

```
ip-chrName.bam
ip-chrName.bed
ip-chrName-lFRAGMENTLENGTH.bed
ip-chrName-lFRAGMENTLENGTH-cov.bed
ip-chrName-lFRAGMENTLENGTH-cov-wWINDOWSIZE-tTHRESHOLD.bed
ip-chrName-lFRAGMENTLENGTH-cov-wWINDOWSIZE-tTHRESHOLD-rpm.bed
ip-chrName-lFRAGMENTLENGTH-cov-wWINDOWSIZE-tTHRESHOLD-rpm_norm.bed
```

## HOW IT WORKS

This script will extract the chromosome of interest (CHRNAME) from the BAM file then convert it into BED format. Each read will be extended to FRAGMENTLENGTH toward their 3' end without exceeding the chromosomes boundaries. The coverage will be extracted (how many reads is found at each position of the chromosome). A smoothing (averaging) will be performed with a sliding window of WINDOWSIZE bp. Each position where the coverage is less than THRESHOLD reads will be replaced by zero. The coverage will be multiplied by the RPM ratio. The IP coverage will be divided base by base by the input coverage.
