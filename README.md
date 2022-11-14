# Int1-epicPCR


This script runs all the steps for the sequence processing and filtering of class 1 integron-16S epicPCR products sequenced with Oxford Nanopore Technologies. The expected epicPCR products are generated from a fusion PCR that physically links class 1 integron gene cassette arrays amplified from single cells, isolated in aqueous phase droplets, to the 16S rRNA V4 region from the same cells. 

As input, this pipeline takes a raw uncompressed fastq file (**WARNING: DO NOT trim adapters during basecalling**), and outputs a single fasta file containing a set of full-length, primer-oriented epicPCR amplicon consensus sequences that contain a complete *attI1* sequence, epciPCR bridging primer, and the 16S rRNA marker gene.

## epicPCR primers
### Expected epicPCR bridging primer used in the initial fusion PCR step:

**R519-*qacE***: 5'-GWATTACCGCGGCKGCTGGCGGAAGGGGCAAGCTTAGT-3' # Bridging primer that fuses class 1 integron cassette array amplicons with 16S amplicons

### Expected nested primer pair used in the final nested PCR step:

***intI1*_nested**: 5'-CGAACGCAGCGGTGGTAA-3' # Nested forward primer that anneals to the class 1 integron-integrase gene (*intI1*).

**AP27_short**: 5'-GCTCTTCCGATCTGGACTACHVGGGTWTCTAAT-3' # Nested reverse primer that anneals to the 16S rRNA gene (R806 annealing site).

***N.B.** If custom primers have been used, the script must be modified with their respective sequences in lines 57-58 (nested primer pair) and line 127 for bridging primer. Also, the length of the bridging primer blast alignment should be modified from '34' in line 133 to suit the custom bridging primer (recommended ~90% of primer length)*


## Dependencies
A system-wide installation (e.g., via conda) is required for the following programs:  
*N.B. The versions listed have been tested using the pipeline, but newer and/or older version will most likely also work*

[NanoFilt](https://github.com/wdecoster/nanofilt) v2.8.0  
[Pychopper](https://github.com/epi2me-labs/pychopper) v2.7.0  
[isONclust](https://github.com/ksahlin/isONclust) v0.0.6.1  
[isONcorrect](https://github.com/ksahlin/isONcorrect) v0.0.8  
[Spoa](https://github.com/rvaser/spoa) v4.0.7  
[BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) v35  
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.2.31  
[SAMtools](http://www.htslib.org/) v1.12  
[Metaxa2](https://microbiology.se/software/metaxa2/) v2.2.3  


## Usage

```
usage: ./epicPCR.sh -r [optional arguments]

Mandatory arguments:
    -r : Nanopore reads as a single uncompressed fastq file
    
Optional arguments:
    -t : number of CPUs | default: 1
    -o : output directory | default: current directory
```
