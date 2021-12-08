#!/bin/bash

IN=../00-rawdata
OUT=../01-cleandata
scripts=../../../scripts

mkdir -p ./tmp
ls $IN | grep .fastq.gz | grep -v ref | while read batch
do

f=$batch
bn=$(echo $f | cut -d'_' -f 1 | python $scripts/parse_name_new.py)

if [ -d $OUT/$bn ];then
    echo "Skipped $bn"
    continue
fi

cat << EOF > ./tmp/_00_clean_${batch}.sh
#!/bin/bash

if [ -f $IN/$f ];then
    mkdir -p $OUT/$bn
fi

#zcat $IN/$f | python $scripts/cut_reads.py 26 | gzip -c > $OUT/$bn/${bn}.fastq.gz
cp $IN/$f $OUT/$bn/${bn}.fastq.gz
#fastqc $OUT/$bn/${bn}.fastq.gz -o $OUT/$bn
EOF
sbatch -c 1 --mem=2g --partition=compute_new -o ./tmp/_00_clean_${batch}.o.txt -e ./tmp/_00_clean_${batch}.e.txt ./tmp/_00_clean_${batch}.sh
done

