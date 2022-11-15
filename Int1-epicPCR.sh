#!/bin/bash

usage() { printf "\nusage: ./Int1-epicPCR.sh -r [optional arguments]\n\nMandatory arguments:\n-r : Nanopore reads as a single fastq file\n\nOptional arguments:\n-t : number of CPUs | default: 1\n-o : output directory | default: current directory\n\nDependencies:\nNanoFilt\nPychopper\nisONclust\nisONcorrect\nSpoa\nBBTools\nBLAST\nSAMtools\nMetaxa2\n\n"; exit 0;}
[ $# -eq 0 ] && usage

while getopts r:t:o: flag
do
case "${flag}"
in
r) reads=${OPTARG};;
t) threads=${OPTARG};;
o) outdir=${OPTARG};;
*) usage exit 0;;
esac
done

if [ -z "$threads" ]
then
threads="1"
fi

if [ -z "$outdir" ]
then
outdir="."
fi

if [ -z "$reads" ]
then
printf "\nMissing: -r <Nanopore reads as fastq>\n" >&2
usage
exit 0
fi

mkdir -p $outdir

#Extract sample name from fastq file
name=$(basename "$reads" | rev | cut -d'.' -f2- | rev);


echo '#################################################'
echo 'Starting NanoFilt'
echo '#################################################'

mkdir -p $outdir/NanoFilt

cat $reads | NanoFilt -q 7 -l 650 > $outdir/NanoFilt/$name.fastq

echo '#################################################'
echo 'Finished NanoFilt'
echo '#################################################'
echo '#################################################'
echo 'Starting Pychopper'
echo '#################################################'

mkdir -p $outdir/Pychopper

for f in $outdir/Pychopper/primer_config.txt; do if [[ ! -e "$f" ]]; then printf "+:intI1N,-16SNT|-:16SNT,-intI1N\n" > $outdir/Pychopper/primer_config.txt; fi; done
for f in $outdir/Pychopper/primers.fas; do if [[ ! -e "$f" ]]; then printf ">intI1N\nCGAACGCAGCGGTGGTAA\n>16SNT\nTGCTCTTCCGATCT\n" > $outdir/Pychopper/primers.fas; fi; done

cdna_classifier.py -m edlib -b $outdir/Pychopper/primers.fas -c $outdir/Pychopper/primer_config.txt -t $threads -p $outdir/NanoFilt/$name.fastq $outdir/Pychopper/$name.pychopper.fastq

echo '#################################################'
echo 'Finished Pychopper'
echo '#################################################'
echo '#################################################'
echo 'Starting isONclust'
echo '#################################################'


mkdir -p $outdir/isONcorrect/$name

isONclust --t $threads --ont --fastq $outdir/Pychopper/$name.pychopper.fastq --outfolder $outdir/isONcorrect/$name/clustering

isONclust write_fastq --N 1 --clusters $outdir/isONcorrect/$name/clustering/final_clusters.tsv --fastq $outdir/Pychopper/$name.pychopper.fastq --outfolder $outdir/isONcorrect/$name/clustering/fastq_files

echo '#################################################'
echo 'Finished isONclust'
echo '#################################################'
echo '#################################################'
echo 'Starting isONcorrect'
echo '#################################################'

run_isoncorrect --t $threads --fastq_folder $outdir/isONcorrect/$name/clustering/fastq_files --outfolder $outdir/isONcorrect/$name/correction/

echo '#################################################'
echo 'Finished isONcorrect'
echo '#################################################'
echo '#################################################'
echo 'Extracting consensus sequences using Spoa'
echo '#################################################'

mkdir -p $outdir/isONcorrect/$name/spoa

mkdir -p $outdir/isONcorrect/$name/spoa/tmp

for i in $outdir/isONcorrect/$name/correction/*
do
	number=$(basename "$i")
	spoa -r 0 $outdir/isONcorrect/$name/correction/$number/corrected_reads.fastq > $outdir/isONcorrect/$name/spoa/tmp/$name.$number.cluster_consensus.fasta
done

# Rename headers of each sequence using its filename

for f in $outdir/isONcorrect/$name/spoa/tmp/*.cluster_consensus.fasta;
do
	id=$(basename "$f" .cluster_consensus.fasta)
	sed 's/>.*/>'$id'.cluster_consensus/' $f > $outdir/isONcorrect/$name/spoa/$id.fasta
done

rm -r $outdir/isONcorrect/$name/spoa/tmp

echo '#################################################'
echo 'Finished spoa consensus'
echo '#################################################'
echo '#################################################'
echo 'Combining all consensuses and removing redundancies (including reverse complement redundancies)'
echo '#################################################'

cat $outdir/isONcorrect/$name/spoa/$name.*.fasta | dedupe.sh in=stdin out=$outdir/isONcorrect/$name.fasta

echo '#################################################'
echo 'BLAST sequences for bridging primer'
echo '#################################################'

mkdir -p $outdir/Bridge_primer/Blast_db

for f in $outdir/Bridge_primer/Bridge_primer.fasta; do if [[ ! -e "$f" ]]; then printf ">Bridge_primer\nGWATTACCGCGGCKGCTGGCGGAAGGGGCAAGCTTAGT\n" > $outdir/Bridge_primer/Bridge_primer.fasta; fi; done

for f in $outdir/Bridge_primer/Blast_db/*.nhr; do if [[ ! -e "$f" ]]; then makeblastdb -in $outdir/Bridge_primer/Bridge_primer.fasta -dbtype nucl -out $outdir/Bridge_primer/Blast_db/bridge -title "bridge"; fi; done

blastn -query $outdir/isONcorrect/$name.fasta -task blastn-short -perc_identity 90 -db $outdir/Bridge_primer/Blast_db/bridge -out $outdir/Bridge_primer/$name.blast_bridge.out -outfmt '6 qseqid qstart qend sstart send length pident evalue mismatch gaps' -num_threads $threads

awk '$6>34' $outdir/Bridge_primer/$name.blast_bridge.out | cut -f1 | sort | uniq -c | awk '$1==1' | awk '{print $2}' | xargs -n 1 samtools faidx $outdir/isONcorrect/$name.fasta > $outdir/Bridge_primer/$name.filt_seqs.fasta

echo '#################################################'
echo 'BLAST sequences for end of intI1 up to and including attI1'
echo '#################################################'

mkdir -p $outdir/attI1/Blast_db

for f in $outdir/attI1/attI1.fasta; do if [[ ! -e "$f" ]]; then printf ">attI1\nCGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAG\n" > $outdir/attI1/attI1.fasta; fi; done

for f in $outdir/attI1/Blast_db/*.nhr; do if [[ ! -e "$f" ]]; then makeblastdb -in $outdir/attI1/attI1.fasta -dbtype nucl -out $outdir/attI1/Blast_db/attI1 -title "attI1"; fi; done

blastn -query $outdir/Bridge_primer/$name.filt_seqs.fasta -task blastn -db $outdir/attI1/Blast_db/attI1 -out $outdir/attI1/$name.blast_attI1.out -outfmt '6 qseqid qstart qend sstart send length pident evalue mismatch gaps' -num_threads $threads

awk '$6>149' $outdir/attI1/$name.blast_attI1.out | cut -f1 | sort | uniq -c | awk '$1==1' | awk '{print $2}' | xargs -n 1 samtools faidx $outdir/Bridge_primer/$name.filt_seqs.fasta > $outdir/attI1/$name.filt_seqs.fasta


echo '#################################################'
echo 'Screen sequences for 16S gene using Metaxa2'
echo '#################################################'

mkdir -p $outdir/metaxa2
mkdir -p $outdir/Final_filtered

metaxa2 -i $outdir/attI1/$name.filt_seqs.fasta -o $outdir/metaxa2/$name --cpu $threads; cut -f1 $outdir/metaxa2/$name.taxonomy.txt | xargs -n 1 samtools faidx $outdir/attI1/$name.filt_seqs.fasta > $outdir/Final_filtered/$name.fasta

echo '#################################################'
echo 'Finished Metaxa2'
echo '#################################################'

echo '#################################################'
printf "Int1-epicPCR.sh has Finished\n\nFinal output saved as $outdir/Final_filtered/$name.fasta"
printf "#################################################\n"



