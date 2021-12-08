#!/bin/bash

echo "condition,library,mutation,n_single_mut_measurements,n_any_mut_measurements,epistasis_model_score,raw_single_mut_score,site,wildtype,mut_escape,protein_chain,single_mut_escape,protein_site,label_site,site_total_escape,site_mean_escape" > ./summary_df.csv

ls ../03-tasks |  while read batch
do
cat ../03-tasks/$batch/effects_df.csv | grep -v protein_chain >> ./summary_df.csv
done

cat ./summary_df.csv | cut -d',' -f 1,8,15,16 | uniq > ./summary_df_site.csv
