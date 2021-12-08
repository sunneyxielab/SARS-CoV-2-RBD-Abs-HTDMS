#!/bin/bash

IN=./datasets/fastq
OUT=./datasets/processed
table=./reference/codon_variant_table_clean.csv
wt=./reference/wildtype_sequence.fasta

scripts=./scripts
tmp=./tmp

mkdir -p $tmp

ls $IN | grep _lib |  while read batch
do
    if [ ! -f $IN/$batch/${batch}.fastq.gz ];then
        echo "Empty $batch"
        rm $IN/$batch -d
        continue
    fi
    if [ -f $OUT/$batch/${batch}_variant_counts.csv ];then
        echo "Skipped $batch"
        continue
    fi

    echo "#!/bin/bash" > $tmp/_${batch}.tmp.sh
    echo "python $scripts/01-align.py -i $IN/$batch/${batch}.fastq.gz -o $OUT -t $table -w $wt" >> $tmp/_${batch}.tmp.sh
    echo $batch

    sbatch --partition=compute_new -c 1 --mem=2g -o $tmp/_${batch}.o.txt -e $tmp/_${batch}.e.txt $tmp/_${batch}.tmp.sh
done
