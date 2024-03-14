echo "********Files with issues***********"

echo "********test_touch_EMPTY.csv***********"

python3 scripts/query_cmap.py ../GatherMetadata/test_touch_EMPTY.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83        

echo "********test_EMPTY.csvv***********"

python3 scripts/query_cmap.py ../GatherMetadata/test_EMPTY.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83                                             
echo "********test_just_headers.csv***********"

python3 scripts/query_cmap.py ../GatherMetadata/test_just_headers.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83                                            

echo "********test_just_headers_incorrect.csv***********"

python3 scripts/query_cmap.py ../GatherMetadata/test_just_headers_incorrect.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83   

echo "********test_incorrect_value.csv.csv***********" 

python3 scripts/query_cmap.py ../GatherMetadata/test_incorrect_value.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83

echo "********test_incorrect_value_ALL.csv.csv***********"

python3 scripts/query_cmap.py ../GatherMetadata/test_incorrect_value_all.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83

echo "******** invalid API key ***********"
python3 scripts/query_cmap.py ../GatherMetadata/test_correct.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc73 





echo "********Files that should work***********"
echo "********test_correct_extra_columns***********"

python3 scripts/query_cmap.py ../GatherMetadata/test_correct_extra_columns.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83

python3 scripts/query_cmap.py ../GatherMetadata/test_correct.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83    

python3 scripts/query_cmap.py ../GatherMetadata/metadata.cmap.csv data/CMAP_metadata.csv.gz  4f3123a0-7fa3-11ec-a96e-6968a44bdc83


