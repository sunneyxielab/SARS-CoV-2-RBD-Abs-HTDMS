#!/bin/bash

reference=../../../reference
scripts=../../../scripts

IN=../02-processed
OUT=../03-tasks

table=$reference/JDBloom_libs/codon_variant_table_clean.csv
wt=$reference/wildtype_sequence.fasta
singlemut=$reference/single_mut_effects.csv
expr=$reference/expression_meanFs.csv
bind=$reference/binding_Kds.csv

tmp=./tmp
siteav=naive

cat ./tasks.csv | grep _lib[12m] | while read task
do
escape=$(echo $task | cut -d , -f 1)
ref=$(echo $task | cut -d , -f 2)
flag=$(echo $task | cut -d , -f 3)
echo $task


# !!! For RD libraries --exprmin and --bindmin arguments are ignored !!! #
cat << EOF > $tmp/_02_${escape}_${ref}.tmp.sh
#!/bin/bash
python $scripts/02-call_escape_score.py \
       -r $IN/$ref/${ref}_variant_counts.csv \
       -e $IN/$escape/${escape}_variant_counts.csv \
       -d $flag -exp MACS \
       -w $wt -t $table \
       -o $OUT \
       -p SARS-CoV-2 -S 330 \
       --mutbindexpr=$singlemut --variant_bind=$bind --variant_expr=$expr \
       --exprmin=-1.0 --bindmin=-2.35 \
       --site_average=$siteav \
       --score=ratio
EOF

sbatch --partition=compute_new -c 1 --mem=2g -o $tmp/_02_${escape}_${ref}.o.txt -e $tmp/_02_${escape}_${ref}.e.txt $tmp/_02_${escape}_${ref}.tmp.sh
done
