#!/bin/bash

IN=../02-processed

echo "sample,valid,lowQ,invalid,unparsable,invalidOrUnparsableRatio,detected_barcodes" > ./barcode_stat.csv

ls $IN | while read batch
do
fates=$IN/$batch/${batch}_fates.csv
vc=$IN/$batch/${batch}_variant_counts.csv

if [ ! -f $fates ];then
    continue
fi

    valid=$(cat $fates | grep -v invalid | grep valid | cut -d ',' -f 3)
    lowQ=$(cat $fates | grep low | cut -d ',' -f 3)
    invalid=$(cat $fates | grep invalid | cut -d ',' -f 3)
    unparsable=$(cat $fates | grep unpars | cut -d ',' -f 3)
    ratio=0`echo "scale=3;(${invalid}+${unparsable})/(${invalid}+${unparsable}+${valid}+${lowQ})" | bc`
    detected_barcodes=$(cat $vc | awk 'BEGIN{FS=",";sum=0}{if($5>0)sum+=1}END{print sum}')

echo "$batch,$valid,$lowQ,$invalid,$unparsable,$ratio,$detected_barcodes" >> ./barcode_stat.csv
done

