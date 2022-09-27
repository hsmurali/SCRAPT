#! /usr/bin/sh


for test_name in small_test.fna test_600k.fna ; do 
    for threads in 1 2 3 4 ; do
    
        #run_cmd="${cmd} --no-k-mer-filter "
        eval time ./dnaclust --no-k-mer-filter -t ${threads} ../dnaclust/${test_name} > results/${test_name}.nk.${threads}.clusters
        awk '{print NF}' results/${test_name}.nk.${threads}.clusters > results/${test_name}.nk.${threads}.counts 

        eval time ./dnaclust -t ${threads} ../dnaclust/${test_name} > results/${test_name}.k.${threads}.clusters
        awk '{print NF}' results/${test_name}.k.${threads}.clusters > results/${test_name}.k.${threads}.counts 

    done
done
