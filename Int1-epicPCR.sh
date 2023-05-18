#!/bin/bash

usage() { 
    printf "\nInt1-epicPCR version 2.0\n\nusage: ./Int1-epicPCR.sh -r [optional arguments]\n\nMandatory arguments:\n-r : Nanopore reads as a single file\n\nOptional arguments:\n-t : number of CPUs | default: 1\n-o : output directory | default: current directory\n\nDependencies:\nNanoFilt\nPychopper\nisONclust\nisONcorrect\nSpoa\nBBMap\nBLAST\nSAMtools\nMetaxa2\n"
    exit 0
}

[ $# -eq 0 ] && usage

while getopts r:t:o: flag; do
    case "${flag}" in
        r) reads=${OPTARG};;
        t) threads=${OPTARG};;
        o) outdir=${OPTARG};;
        *) usage exit 0;;
    esac
done

if [ -z "$threads" ]; then
    threads="1"
fi

if [ -z "$outdir" ]; then
    outdir="."
fi

if [ -z "$reads" ]; then
    printf "\nMissing: -r <Nanopore reads as a single file>\n" >&2
    usage
    exit 0
fi


# Check for required software tools
if ! command -v NanoFilt &> /dev/null
then
    echo "ERROR: NanoFilt is not installed. Please install it and try again."
    exit 1
fi

if ! command -v cdna_classifier.py &> /dev/null
then
    echo "ERROR: Pychopper is not installed. Please install it and try again."
    exit 1
fi


if ! command -v isONclust &> /dev/null
then
    echo "ERROR: isONclust is not installed. Please install it and try again."
    exit 1
fi

if ! command -v run_isoncorrect &> /dev/null
then
    echo "ERROR: isONcorrect is not installed. Please install it and try again."
    exit 1
fi

if ! command -v spoa &> /dev/null
then
    echo "ERROR: Spoa is not installed. Please install it and try again."
    exit 1
fi

if ! command -v dedupe.sh &> /dev/null
then
    echo "ERROR: BBMap is not installed. Please install it and try again."
    exit 1
fi

if ! command -v blastn &> /dev/null
then
    echo "ERROR: BLAST is not installed. Please install it and try again."
    exit 1
fi

if ! command -v samtools &> /dev/null
then
    echo "ERROR: SAMtools is not installed. Please install it and try again."
    exit 1
fi

if ! command -v metaxa2 &> /dev/null
then
    echo "ERROR: Metaxa2 is not installed. Please install it and try again."
    exit 1
fi


# Ensure that the script will exit immediately if any of the commands in the pipeline fail
set -e

# Create pipeline output directory
mkdir -p $outdir

# check if the input file is compressed, and uncompressed it with gunzip if necessary
if [[ $reads == *.gz ]]; then
    gunzip -c $reads > "${reads%.gz}"
    reads="${reads%.gz}"
fi

#Extract sample name from fastq file
name=$(basename "$reads" | rev | cut -d'.' -f2- | rev);

echo '#################################################'
echo 'Starting NanoFilt'
echo '#################################################'

# Create NanoFilt output directory
mkdir -p $outdir/NanoFilt

# Run NanoFilt
cat $reads | NanoFilt -q 7 -l 650 > $outdir/NanoFilt/$name.fastq

# Check the exit status of NanoFilt and exit if it failed
if [ $? -ne 0 ]; then
    echo "Error: NanoFilt failed"
    exit 1
fi

echo '#################################################'
echo 'Finished NanoFilt'
echo '#################################################'
echo '#################################################'
echo 'Starting Pychopper'
echo '#################################################'

# Create Pychopper output directory
mkdir -p $outdir/Pychopper

# If primer_config.txt file does not exist, create it with default primer Ids
for f in $outdir/Pychopper/primer_config.txt; do
    if [[ ! -e "$f" ]]; then
        printf "+:intI1N,-16SNT|-:16SNT,-intI1N\n" > $outdir/Pychopper/primer_config.txt
    fi
done

#If primers.fas file does not exist, create it with default primer sequences
for f in $outdir/Pychopper/primers.fas; do
    if [[ ! -e "$f" ]]; then
        printf ">intI1N\nCGAACGCAGCGGTGGTAA\n>16SNT\nGCTCTTCCGATCTGGACTAC\n" > $outdir/Pychopper/primers.fas
    fi
done

# Run Pychopper
cdna_classifier.py -m edlib -b $outdir/Pychopper/primers.fas \
-c $outdir/Pychopper/primer_config.txt -t $threads \
-p $outdir/NanoFilt/$name.fastq $outdir/Pychopper/$name.pychopper.fastq

# Check the exit status of Pychopper and exit if it failed
if [ $? -ne 0 ]; then
    echo "Error: Pychopper failed"
    exit 1
fi

# Move cdna_classifier_reports
for f in cdna_classifier_report.*; do 
  mv $f $outdir/Pychopper/$name.$f
done

echo '#################################################'
echo 'Finished Pychopper'
echo '#################################################'
echo '#################################################'
echo 'Starting isONclust'
echo '#################################################'

# Create output directory fro isONclust and isONcorrect
mkdir -p $outdir/isONcorrect/$name

# Run isONclust
isONclust --t $threads --ont --fastq $outdir/Pychopper/$name.pychopper.fastq \
  --outfolder $outdir/isONcorrect/$name/clustering

# Write clusters to fastq files
isONclust write_fastq --N 1 \
  --clusters $outdir/isONcorrect/$name/clustering/final_clusters.tsv \
  --fastq $outdir/Pychopper/$name.pychopper.fastq \
  --outfolder $outdir/isONcorrect/$name/clustering/fastq_files

# Check the exit status of isONclust and exit if it failed
if [ $? -ne 0 ]; then
  echo "Error: isONclust failed"
  exit 1
fi

echo '#################################################'
echo 'Finished isONclust'
echo '#################################################'
echo '#################################################'
echo 'Starting isONcorrect'
echo '#################################################'

# Run isONcorrect
run_isoncorrect --t $threads \
  --fastq_folder $outdir/isONcorrect/$name/clustering/fastq_files \
  --outfolder $outdir/isONcorrect/$name/correction/

# Check the exit status of isONcorrect and exit if it failed
if [ $? -ne 0 ]; then
  echo "Error: isONcorrect failed"
  exit 1
fi

echo '#################################################'
echo 'Finished isONcorrect'
echo '#################################################'
echo '#################################################'
echo 'Extracting consensus sequences using Spoa'
echo '#################################################'

# Create Spoa output directory, including temporary working directory
mkdir -p $outdir/isONcorrect/$name/spoa
mkdir -p $outdir/isONcorrect/$name/spoa/tmp

# Generate consensus sequence for each cluster using Spoa
for i in $outdir/isONcorrect/$name/correction/*; do
    number=$(basename "$i")
    spoa -r 0 $outdir/isONcorrect/$name/correction/$number/corrected_reads.fastq > $outdir/isONcorrect/$name/spoa/tmp/$name.$number.cluster_consensus.fasta
done

# Check the exit status of Spoa and exit if it failed
if [ $? -ne 0 ]; then
    echo "Error: Spoa failed"
    exit 1
fi

# Rename headers of each sequence using its filename
for f in $outdir/isONcorrect/$name/spoa/tmp/*.cluster_consensus.fasta; do
    id=$(basename "$f" .cluster_consensus.fasta)
    sed 's/>.*/>'$id'.cluster_consensus/' $f > $outdir/isONcorrect/$name/spoa/$id.fasta
done

# Remove temp files
rm -r $outdir/isONcorrect/$name/spoa/tmp

echo '#################################################'
echo 'Finished Spoa consensus'
echo '#################################################'
echo '#################################################'
echo 'Combining all consensuses and removing redundancies (including reverse complement redundancies)'
echo '#################################################'

# Combined all consensus sequences and remove any sequence redundancies suing BBMap's dedupe.sh
cat $outdir/isONcorrect/$name/spoa/$name.*.fasta | dedupe.sh in=stdin out=$outdir/isONcorrect/$name.fasta

# Check the exit status of dedupe.sh and exit if it failed
if [ $? -ne 0 ]; then
    echo "Error: dedupe.sh failed"
    exit 1
fi

echo '#################################################'
echo 'BLAST sequences for bridging primer'
echo '#################################################'

# Create directory for Blast database
mkdir -p $outdir/Bridge_primer/Blast_db

# If Bridge_primer.fasta file does not exist, create it with default sequence
for f in $outdir/Bridge_primer/Bridge_primer.fasta; do
    if [[ ! -e "$f" ]]; then
        printf ">Bridge_primer\nGWATTACCGCGGCKGCTGGCGGAAGGGGCAAGCTTAGT\n" > $outdir/Bridge_primer/Bridge_primer.fasta
    fi
done

# If Blast database does not exist, create it
for f in $outdir/Bridge_primer/Blast_db/*.nhr; do
    if [[ ! -e "$f" ]]; then
        makeblastdb -in $outdir/Bridge_primer/Bridge_primer.fasta -dbtype nucl -out $outdir/Bridge_primer/Blast_db/bridge -title "bridge"
    fi
done

# Run Blastn with output format 6 and check exit status

blastn -query $outdir/isONcorrect/$name.fasta -task blastn-short -perc_identity 90 \
-db $outdir/Bridge_primer/Blast_db/bridge -out $outdir/Bridge_primer/$name.blast_bridge.out \
-outfmt '6 qseqid qstart qend sstart send length pident evalue mismatch gaps' -num_threads $threads

# Check the exit status of Blastn and exit if it failed
if [ $? -ne 0 ]; then
    echo "Error: Blastn failed"
    exit 1
fi

# Filter Blast output and extract sequences using samtools

awk '$6>34' $outdir/Bridge_primer/$name.blast_bridge.out | cut -f1 | sort | uniq -c | awk '$1==1' | awk '{print $2}' | xargs -n 1 samtools faidx $outdir/isONcorrect/$name.fasta > $outdir/Bridge_primer/$name.filt_seqs.fasta

echo '#################################################'
echo 'BLAST sequences for end of intI1 up to and including attI1'
echo '#################################################'

# Create directory for Blast database
mkdir -p "$outdir/attI1/Blast_db"

# If attI1.fasta file does not exist, create it with default sequence
for f in "$outdir/attI1/attI1.fasta"; do
  if [[ ! -e "$f" ]]; then
    printf ">attI1\nCGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAG\n" > "$outdir/attI1/attI1.fasta"
  fi
done

# If Blast database does not exist, create it
for f in "$outdir/attI1/Blast_db/"*.nhr; do
  if [[ ! -e "$f" ]]; then
    makeblastdb -in "$outdir/attI1/attI1.fasta" -dbtype nucl -out "$outdir/attI1/Blast_db/attI1" -title "attI1"
  fi
done

# Run Blastn with output format 6 and check exit status
blastn -query "$outdir/Bridge_primer/$name.filt_seqs.fasta" -task blastn \
-db "$outdir/attI1/Blast_db/attI1" -out "$outdir/attI1/$name.blast_attI1.out" \
-outfmt '6 qseqid qstart qend sstart send length pident evalue mismatch gaps' -num_threads $threads

if [ $? -ne 0 ]; then
  echo "Error: Blastn failed"
  exit 1
fi

# Filter Blast output and extract sequences using samtools
awk '$6>149' "$outdir/attI1/$name.blast_attI1.out" | cut -f1 | sort | uniq -c | awk '$1==1' | awk '{print $2}' | xargs -n 1 samtools faidx "$outdir/Bridge_primer/$name.filt_seqs.fasta" > "$outdir/attI1/$name.filt_seqs.fasta"

echo '#################################################'
echo 'Screen sequences for 16S gene using Metaxa2'
echo '#################################################'

# Create directories for Metaxa2 output and final filtered sequences
mkdir -p "$outdir/metaxa2"
mkdir -p "$outdir/Final_filtered"

# Run Metaxa2 and check exit status
metaxa2 -i "$outdir/attI1/$name.filt_seqs.fasta" -o "$outdir/metaxa2/$name" --cpu $threads

if [ $? -ne 0 ]; then
  echo "Error: Metaxa2 failed"
  exit 1
fi

#Extract  final filtered sequences
cut -f1 "$outdir/metaxa2/$name.taxonomy.txt" | xargs -n 1 samtools faidx "$outdir/attI1/$name.filt_seqs.fasta" > "$outdir/Final_filtered/$name.fasta"

echo '#################################################'
echo 'Finished Metaxa2'
echo '#################################################'

echo '#################################################'
echo 'epicPCR.sh has Finished'
echo '#################################################'

