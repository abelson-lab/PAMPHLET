# PAMPHLET

Here we present a manual with examples detailing how to use PAMPHLET to generate recurrently mutated genomics
coordinates for targeted sequencing panels based on prior knowledge, such as from mutation data obtained from the
Catalogue of Somatic Mutations in Cancer (COSMIC) and cBio Cancer Genomics Portal (CBioPortal).

## Quick Start

A user can execute the code of PAMPHLET from the terminal. There is no need for installation. Copy the code files to
your folder of choice.

## Dependencies

This pipeline requires the following software and packages:

| Program                         | Packages                                            |
|---------------------------------|-----------------------------------------------------|
| Python (https://www.python.org) | argparse, matplotlib.pyplot, csv, re, typing, panda |        

## Substitution, Insertions, and Deletions

### To download the COSMIC mutation file

1. Go to COSMIC (https://cancer.sanger.ac.uk/cosmic) and log in. Registration might be required
2. On the download page (https://cancer.sanger.ac.uk/cosmic/download), download COSMIC mutation data that includes both
   targeted and genome-wide screens. On the website, it is named 'CosmicMutantExport.tsv.gz'
3. To run the example, use either the 'Download Whole file' option or the 'Download Filtered File' option, and fill in '
   haematopoietic_and_lymphoid_tissue' for the 'Filter by cancer' option.

### COSMIC mutation file description

The COSMIC mutation file is a table of mutations associated with information that is presented in columns such as
mutation CDS (the nucleotide mutation), mutation description (amino acid level mutation type), genomics coordinate,
cancer type, gene name, tumour id. Each entry in the file is linked to a single tumour ID and one genomic coordinate.

### To download the gene track file

1. Go to the UCSC Table browser (https://genome.ucsc.edu/cgi-bin/hgTables)
2. Select the following options; 'Genes and Gene predictions' in group, 'HGNC' in track and 'hgnc' in table.

### Gene track file description

The gene track list the genomic coordinate for the coding region and transcription region of each gene.

### Run with example data

```
python3 PAMPHLET.py -t sub_indel -m <cosmic mutation file path> -r <refseq gene file path> -d default -a m 
```

### Tool Description

#### Step 0: Read file and Specify User Preferences

PAMPHLET takes a COSMIC mutation file and an UCSC gene track as input to output range that targets substitutions and
small insertions and deletions (small indels). Users can choose only to target mutations that originated from any
cancer. Separately, the PAMPHLET provides an option to cover the entire coding region of any gene in addition to the
ranges covering substitutions and small indels, which is called "cover entire gene."

<p style="text-align: center;"> <strong>Cover Entire Gene:</strong> Whether to cover the entire coding region of any gene.</p> 

#### Step 1: Cancer Type

First, PAMPHLET prompts the user to specify what cancer type the user wishes to target. For running with example data,
the default cancer type will be myeloid cancers.

For the example, the following is displayed in the terminal.

```
searching the mutation file to choose cancer
['haematopoietic_and_lymphoid_tissue']
['haematopoietic_neoplasm']
['myelodysplastic-myeloproliferative_neoplasm', 'acute_myeloid_leukaemia',
'acute_leukaemic_transformation_of_chronic_myelomonocytic_leukaemia', 'chronic_myeloid_leukaemia',
'myeloproliferative_neoplasm', 'myelodysplastic_transformation_of_essential_thrombocythaemia',
'myelodysplastic-myeloproliferative_disease-unclassifiable', 'acute_myeloid_leukaemia_associated_with_mastocytosis',
'acute_myeloid_leukaemia_associated_with_MDS', 'myeloid_neoplasm_unspecified_therapy_related',
'myelodysplastic-myeloproliferative_neoplasm-unclassifiable', 'primary_myelofibrosis',
'blast_phase_chronic_myeloid_leukaemia', 'plasma_cell_myeloma',
'chronic_myelomonocytic_leukaemia_therapy_related', 'myelodysplastic_syndrome',
'juvenile_myelomonocytic_leukaemia', 'acute_leukaemic_transformation_of_myeloproliferative_neoplasm',
'acute_myeloid_leukaemia_therapy_related', 'myelofibrosis', 'acute_leukaemic_transformation_of_primary_myelofibrosis',
'myelodysplastic_syndrome_therapy_related', 'chronic_myelomonocytic_leukaemia',
'myeloproliferative_neoplasm_unclassifiable', 'acute_myeloid_leukaemia_myelodysplastic_syndrome_therapy_related_NOS']
```

#### Step 2: Filtering

Second, it applies user-defined filters to remove synonymous mutations through mutation description, non-exonic
mutations through mutation CDS, and large indels through genomics coordinates. 

<p style="text-align: center;"> <strong>Indel Filtering Threshold:</strong> Any indels larger than this number of bases
will be filtered out.</p> 

<p style="text-align: center;"> <strong>Include Non-Exonic Mutation:</strong> Whether the user wants to include mutations
not in exons.</p> 

For the example, in the terminal, it lists the number of mutations and tumours before and after each filter.

```
total number of mutations before intronic filter 167483
total number of tumours before intronic filter 63842
total number of mutations after intronic before indel filter 47951
total number of tumours after intronic before indel filter 63761
total number of mutations after indel filter before recurrent 47137
total number of tumours after indel filter before recurrent 34643
```

#### Step 3: Merge Entries

Thirdly, PAMPHLET merges entries that share the same genomic location and counts the number of tumours in which each
mutated location is found. Non-recurrent mutations are excluded.

<p style="text-align: center;"><strong>Recurrent Definition:</strong> Are two mutations occurring
at the same coordinate necessary for mutation at the coordinate to be considered recurrent, or three?</p> 

For the example, recurrent mutations are displayed in a similar in step 2.

```
total number of recurrent mutations after indel filter 11214
total number of tumours with only recurrent mutations after indel filter 33078
```

#### Step 4: Rank Mutations

Fourth, it ranks each unique mutation by the number of tumours the mutations are reported in, with a decreasing order.
At this point, PAMPHLET produces a table ranking the mutations by the number of tumours, where the most recurrent
mutation is at the top, along with their locations, amino acid changes, corresponding tumour sets, and corresponding
genes. See the example intermediate output below.
![img.png](./new%20intermediate%20output.png)

#### Step 5: Shift and Delete

In the fifth step, PAMPHLET begins with the most recurrent unique mutation. It searches for nearby unique mutations to
find the optimal range based on the sum of all mutation recurrences in the ranges. Nearby unique mutation means that
unique mutation is close enough so that the same targeting window can cover that unique mutation and the most recurrent
unique mutation. PAMPHLET only considers nearby recurrent mutations defined by user input.

<p style="text-align: center;"> <strong>Targeting Window Size:</strong> The number of bases one probe from the user's sequencing method
of choice can reliably sequence. </p> 

To maximize the number of mutations covered by each range, PAMPHLET does not consider, in subsequent ranges, mutations
that are already covered by preceding ranges. PAMPHLET repeats the same process for the second-most un-removed recurrent
mutation, the third, and so on. This process is repeated until it reaches the cumulative contribution threshold.

<p style="text-align: center;"> <strong>Cumulative Contribution Threshold:</strong> The percentage of tumours
covered. Tumours without any recurrent mutations are omitted.</p> 

If the 'merge other' option is turned off, PAMPHLET does not search for nearby mutations. It simply uses one range per
mutation.
<p style="text-align: center;"> <strong>Merge Other</strong> is whether the user chooses to
cover more than one mutation per range </p> 

The final output includes a list of optimal ranges, the mutation each range cover, their corresponding tumours and
corresponding gene, and a subset for ranges that covers indels. Users can use the final output to find ranges that cover
X% of tumours based on the cumulative contribution. See the example final output below.
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
2. On the download page (https://cancer.sanger.ac.uk/cosmic/download), download the file for "Copy Number Variants.". On
   the website, it is named 'CosmicCompleteCNA.tsv.gz'
3. To run the example, use either the 'Download Whole file' option or the 'Download Filtered File' option, and fill in '
   haematopoietic_and_lymphoid_tissue' for the 'Filter by cancer' option.

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

### Run with example data

The example only works with COSMIC file.

```
python3 PAMPHLET.py -t CNV -s cosmic -m <cosmic mutation file path> -r <refseq gene file path> -p <common snp file path> 
-d default -a m
```

### Tool description

The user can use the COSMIC CNV file, the cBioPortal CNV file, or provide a gene list. It targets bi-allele heterozygous
SNP because tho allelic frequency of those SNP sites shifts away from 0.5.

1. PAMPHLET's CNV workflow starts by first sorting the genes by the number of tumours in which they were found to be
   involved. These are the CNV genes. See the example ranking below.
   ![img.png](./CNV%20genes%20ranked.png)
2. It then calculates Required MAF by dividing Heterozygous SNP Number by Number of Probes

<p style="text-align: center;"> <strong>Required MAF:</strong> The required population minor allelic frequency </p> 
<p style="text-align: center;"> <strong>Heterozygous SNP Number:</strong> The number of SNP all individuals have on average
in each targeted gene with the heterozygous genotype that the user expects </p> 
<p style="text-align: center;"> <strong>Number of Probes Number:</strong> the number of probes supplied per gene per 
individual </p> 

3. To illustrate, suppose only 100 probes are supplied per individual per gene, and assume one probe targets at least 1
   SNP site. If the user wishes to obtain informative output in 5 of those SNP sites for every individual, in that case,
   the MAF needs to be around 5%.
4. Next, for the top CNV genes, find SNPs around the MAF. It excludes all simple nucleotide polymorphism, except SNP,
through the class column in dbSNP file. As well, it excludes non bi-allele SNP through allele frequency count.

<p style="text-align: center;"> <strong>Top CNV Genes:</strong> For example, the top 10 CNV genes by number of tumour 
or samples it is found in.</p>

5. After SNPs are found, they undergo a similar process to shift and delete. However, this time it is the number of SNPs
   instead of the sum of mutation recurrence. This step gives the snp ranges. 
6. If all SNPs within the required MAF range are covered, and there are still probes left, then they are placed in the
   part of the gene where a range has yet to be placed. This step gives the read cov ranges. 

The final output is a list of range that targets the SNP of CNV genes. First and last snp is the first and last snp of
each range. Read coverage ranges (starts) refer to the starting coordinate of each read cov range. See the example final
output below
![img.png](./CNV%20workflow%20final%20output.png)

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

(if user is providing their own list for CNV)

```
python3 PAMPHLET.py -t CNV -s user -c <cosmic CNV file path>  -r <refseq gene file path> -p <common snp file path>
```

