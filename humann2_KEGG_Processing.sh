for i in `cat LIST.txt`
do
        echo ${i}_genefamilies.tsv
        humann2_regroup_table --input ${i}_genefamilies.tsv -c /home/tarini/KEGG_MAPPINGS/map_ko_uniref90.txt.gz --output ${i}_KO.tsv >> KO_file.log 2>&1
done

for i in `cat LIST.txt`
do
        humann2 --input ${i}_KO.tsv --output KEGG_output_directory --pathways-database /mnt/workspace/tarini/TravellerHealthMET_Shotgun-111190079/curated_KEGG_mapping --input-format genetable
done

humann2_join_tables --input KEGG_Pathabundance/ --output TM_KEGG_Pathways.txt
grep -v "UNMAPPED\|UNINTEGRATED" TM_KEGG_Pathways.txt > TM_KEGG_Pathways_Filtered.txt
humann2_split_stratified_table --input TM_KEGG_Pathways_1.txt --output .


