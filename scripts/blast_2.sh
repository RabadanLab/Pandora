#!/bin/bash -x

# grab arguments
blast_type=$1
prefix=$2
db=$3
threads=$4

echo [start]
echo [pwd] `pwd`
echo [date] `date -u`

# for input: prefix + NUM.fasta
# for output: prefix + NUM.result
# divide threads by two for each blast query

fmt="6 qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle"

echo [database] ${db}

for i in {1..4}; do
	
	echo ${i}

	input=${prefix}_${i}.fasta
	output=${prefix}_${i}.result
	
	echo "input "${input}
	echo "output "${output}
	
	{ ${blast_type} -outfmt "${fmt}" -query ${input} -db ${db} -num_threads $((${threads}/4)) > ${output} & }
	pids[${i}]=$!
done


# wait for all pids
echo [wait-start] 

for pid in ${pids[*]}; do
    wait $pid
done

echo [wait-end]
echo [finish]
echo [date] `date -u`
