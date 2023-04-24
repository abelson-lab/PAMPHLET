# PAMPHLET

Here we present a manual detailing how to use PAMPHLET to generate recurrently mutated genomics coordinates for targeted
sequencing panels based on prior knowledge, such as from mutation data obtained from the Catalogue of Somatic Mutations
in Cancer (COSMIC)
and cBio Cancer Genomics Portal (CBioPortal).

## Quick Start

A user can execute the code of PAMPHLET from the terminal. There is no need for installation. Copy the code files to
your folder of choice.

## Dependencies

This pipeline requires the following software and packages:

| Program                         | Packages                                     |
|---------------------------------|----------------------------------------------|
| Python (https://www.python.org) | argparse, matplotlib.pyplot, csv, re, typing, panda |        

## Substitution and Indels

### To download the COSMIC mutation file

1. Go to COSMIC (https://cancer.sanger.ac.uk/cosmic) and log in. Registration might be required
2. On the download page (https://cancer.sanger.ac.uk/cosmic/download) download COSMIC mutation data that includes both 
targeted and genome-wide screens. On the website it is named 'CosmicMutantExport.tsv.gz'

### COSMIC mutation file description

The COSMIC mutation file is a table of mutations associated with information that is presented in columns such as
mutation CDS (the nucleotide mutation), mutation AA (the amino acid mutation), mutation description (amino acid level
mutation type), genomics coordinate, cancer type, gene name, tumour id. Each entry in the file is linked to a single
tumour ID and one genomic coordinate.

### To download the gene track file

1. Go to the UCSC Table browser (https://genome.ucsc.edu/cgi-bin/hgTables)
2. Select the following options; 'Genes and Gene predictions' in group, 'HGNC' in track and 'hgnc' in table.

### Gene track file description

The gene track list the genomic coordinate for the coding region and transcription region of each gene.

### Run with example data
```
python3 PAMPHLET.py -t sub_indel -m <cosmic mutation file path> -r <refseq gene file path> -d default -a m 
```
To run the example, either use the 'Download Whole file' option or use the 'Download Filtered File' option, but also
was fill in 'haematopoietic_and_lymphoid_tissue' for the 'Filter by cancer' option.  

### Tool Description

#### Step 0: Read file and Specify User Preferences
PAMPHLET takes a cosmic mutation file and an UCSC gene track as input to output range that targets substitutions and
small indels. Users can choose only to target mutations that originated from any cancer. Separately, the PAMPHLET
provides an option to cover the entire coding region of any gene in addition to the ranges covering substitutions and
small indels, which is called "cover entire gene."

#### Step 1: Cancer Type
First, PAMPHLET prompts the user to specify what cancer type the user wishes to target. 

#### Step 2: Filtering
Second, it applies user-defined filters to remove synonymous mutations, non-exonic mutations, and large indels.

#### Step 3: Merge Entries
Thirdly, PAMPHLET merges entries that share the same genomic location and counts the number of tumours in which each 
mutated location is found.

#### Step 4: Rank Mutations
Fourth, it ranks each unique mutation by the number of tumours the mutations are reported in, with a decreasing order.
At this point, PAMPHLET produces a table ranking the mutations by the number of tumours, where the most recurrent
mutation is at the top, along with their locations, amino acid changes, corresponding tumour sets, and corresponding
genes. See the intermediate output below.
![img.png](./new%20intermediate%20output.png)

#### Step 5: 
In the fifth step, PAMPHLET begins with the most recurrent unique mutation. It searches for nearby unique mutations to
find the optimal range based on the sum of all mutation recurrences in the ranges. Nearby unique mutation means that
unique mutation is close enough so that the same targeting window can cover that unique mutation and the most recurrent
unique mutation. 

<p style="text-align: center;"> **Targeting window size** is a user input on how many bases one can a probe from the user's sequencing method
of choice reliably sequence. </p> 

PAMPHLET only considers nearby recurrent mutations defined by user input. Are two mutations occurring
at the same coordinate necessary to be considered recurrent, or three? 

To maximize the number of mutations covered by
each range, PAMPHLET does not consider, in subsequent ranges, mutations that are already covered by preceding ranges.
PAMPHLET repeats the same process for the second-most un-removed recurrent mutation, the third, and so on. This process
is repeated until it reaches the cumulative contribution threshold, defined as a user-specified percentage of tumours
covered. Tumours without any recurrent mutations are omitted. The 'merge other' parameter is whether the user chooses to
cover more than one mutation per range, so if this is turned off, PAMPHLET does not search for nearby mutations. It
simply uses one range per mutation.

The final output includes a list of optimal ranges and the mutation they cover along with their corresponding tumours
and corresponding gene, and a subset for ranges that covers indels. Users can use the final output to find ranges that
cover X% of tumours based on the cumulative contribution. See the final output below.
![img.png](./final%20output%20(all).png)
![img.png](./final%20output%20indel.png)

### Running the code

Go to the directory where you put the code files and type into the terminal.

```
python3 PAMPHLET.py -t sub_indel -m <cosmic mutation file path> -r <refseq gene file path>
```

## CNV

### To download the COSMIC mutation file

1. Go to COSMIC (https://cancer.sanger.ac.uk/cosmic) and log in. Registration might be required
2. Download the file for "Copy Number Variants."

### COSMIC CNV File description

The important columns of this file are nearly identical to the mutation file, with the exception that the genomic
coordinates refer to the CNV in which the gene is involved rather than the gene itself.

### To download the cbioportal CNV file

1. Go to cBioPortal (https://www.cbioportal.org)
2. Query the cancer studies in the tissue filter the user would like to target
3. Filter based on user choice of cancer type and others, and download CNV genes

### cBioPortal CNV File description

In the cBioPortal CNV file, the number of samples with a CNV column indicates the number of samples in which CNV was
found for each gene out of the samples analyzed for CNVs. Since only some samples were analyzed for every gene, this
number may differ for some genes. The file also contains the gene name.

### To download the dbSNP file

1. Go to the UCSC Table browser (https://genome.ucsc.edu/cgi-bin/hgTables)
2. Select the following options; group (Variation), track (common SNPs(151)), table (snp151Common)

### dbSNP file description

The dbSNP file lists single nucleotide polymorphisms (SNP). Each SNP also includes the chromosome, start and end
positions, class, allele frequency count, and major and minor allelic frequencies. Class indicates the simple nucleotide
polymorphism type, such as single nucleotide variation. Allele frequency count indicates the number of alleles.

### Tool description

The user can use the COSMIC CNV file, the cBioPortal CNV file, or provide a gene list. For a bi-allele heterozygous SNP
in the control or a healthy sample, the allelic frequency of any SNP site is expected to be 0.5 since half of all the
bases from all cells is the major allele, and the other half is the minor allele. Allelic frequency is defined as the
number of minor alleles divided by the total number. However, in a tumour sample, if there is a copy number change
affecting the heterozygous site, the allelic frequency shifts away from 0.5. We cannot use homozygous SNPs because their
allelic frequency is not 0.5. It would be either 1.0 or 0.0. The difference between point mutations and CNVs is that
shifts in allelic frequency will occur at more sites in CNV, and these shifts will occur in sites that are concentrated
or close to each other.

PAMPHLET's CNV workflow starts by first sorting the genes by the number of tumours in which they were found to be
involved. It then calculates the required population minor allelic frequency (required MAF) by dividing the number of
heterozygous SNP per gene per individual that the user desires (heterozygous SNP number, a user input) by the number of
probes supplied per gene per individual (number of probes, another user input). The heterozygous SNP number is the
expected number all individuals have on average with the heterozygous genotype. To illustrate, suppose only 100 probes
are supplied per individual per gene, and assume one probe targets at least 1 SNP site. If the user wishes to obtain
informative output in 5 of those SNP sites for every individual, in that case, the MAF needs to be around 5%.

Next, for the top X CNV gene, a user-defined input, find SNPs around the MAF. After SNPs are found, they undergo a
similar process to shift and delete. However, this time it is the number of SNPs instead of the sum of mutation
recurrence. If all SNPs within the required MAF range are covered, and there are still probes left, then they are placed
in the part of the gene where a range has not been placed.

The final output is a list of range that targets single nucleotide polymorphism (SNP)
located in the transcription region of genes recurrently involved in copy number variation (CNV)
in cancers of the user's choice within the required MAF range

## Running the code

Go to the directory where you put the code files and type into the terminal.
(if user is using cbioportal for CNV)

```
python3 PAMPHLET.py -t CNV -s cbioportal -c <cbioportal CNV file path> -r <refseq gene file path> -p <common snp filename>
```

(if user is using cosmic for CNV)

```
python3 PAMPHLET.py -t CNV -s cosmic -c <cosmic CNV file path>  -r <refseq gene file path> -p <common snp file path>
```

