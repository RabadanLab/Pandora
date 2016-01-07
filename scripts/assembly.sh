#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=12G,time=12::
#$ -pe smp 8

# This script performs assembly on the reads leftover after host removal

noclean=${1}	# no clean boolean

echo "------------------------------------------------------------------"
echo ASSEMBLY START [[ `date` ]]

mkdir -p assembly_trinity assembly

Trinity --seqType fq --max_memory 50G --CPU 8 --normalize_reads --output assembly_trinity --left host_separation/unmapped_1.fastq.gz --right host_separation/unmapped_2.fastq.gz
sed -e 's/\(^>.*$\)/#\1#/' assembly_trinity/Trinity.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | awk '/^>/{print ">contig_" ++i; next}{print}' > assembly/contigs_trinity.fasta

# get number of contigs
NUM_CONTIGS=$(( `cat assembly/contigs_trinity.fasta | wc -l` / 2 ))

# compute simple distribution
cat assembly/contigs_trinity.fasta | paste - - | awk '{print length($2)}' | sort -nr | ${d}/scripts/tablecount | awk -v tot=${NUM_CONTIGS} 'BEGIN{x=0}{x+=$2; print $1"\t"$2"\t"x"/"tot"\t"int(100*x/tot)"%"}' > assembly/contigs.distrib.txt

echo clean up
if [ ${noclean} -eq 0 ]; then
	rm -r assembly_trinity
fi

echo ASSEMBLY END [[ `date` ]]
echo "------------------------------------------------------------------"