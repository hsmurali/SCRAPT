from os import listdir, mkdir
from os.path import isdir, isfile
import pandas as pd
import numpy as np
from io import BufferedReader
import sys

scrapt_sim = sys.argv[1]
alpha = sys.argv[2]
dada2_sim = sys.argv[3]

#data_dir = '/fs/cbcb-lab/mpop/projects/SCRAPT/'
#counts_path = data_dir + 'Datasets/Tara_Oceans_Polar/Reads_410_BP/Counts.dict'
#scrapt_path = data_dir + 'Experiments/Tara_Oceans_Polar/Adaptive_With_Modeshifting/sim_'+scrapt_sim+'/alpha_'+alpha+'/'
#dada2_centroid = data_dir + 'Experiments/Tara_Oceans_Polar/DADA2_Benchmarks/dada2_centroids.fna' 
#dada2_results = data_dir+'Experiments/Tara_Oceans_Polar/DADA2_Benchmarks/tara_ocean_larger_dts_dada2_centroids_counts.txt'
#dada2_scrapt = data_dir+'Experiments/Tara_Oceans_Polar/DADA2_Benchmarks/SCRAPT_DADA2/sim_'+scrapt_sim+'_alpha_'+alpha+'_'+dada2_sim
#outdir = data_dir+'Experiments/Tara_Oceans_Polar/DADA2_Benchmarks/Merged_DADA2_SCRAPT/SCRAPT_'+scrapt_sim+'_Alpha_'+alpha+'_DADA2_'+dada2_sim

data_dir = '/fs/cbcb-lab/mpop/projects/SCRAPT/'
counts_path = data_dir+'Datasets/Earth_Microbiome/Counts.soil.dict'
scrapt_path = data_dir+'Experiments/Spatil_Soil/Adaptive_With_Modeshifting/sim_'+scrapt_sim+'/alpha_'+alpha+'/'
dada2_centroid = data_dir+'Experiments/Spatil_Soil/Dada2Benchmark_8/scrapt_vs_dada2/dada2_soil_join.fna' 
dada2_results = data_dir+'Experiments/Spatil_Soil/Dada2Benchmark_8/soil_joining'
dada2_scrapt = data_dir+'Experiments/Spatil_Soil/Dada2Benchmark_8/DADA2_SCRAPT/sim_'+scrapt_sim+'_alpha_'+alpha+'_'+dada2_sim
outdir = data_dir+'Experiments/Spatil_Soil/Dada2Benchmark_8/Merged_DADA2_SCRAPT/SCRAPT_'+scrapt_sim+'_Alpha_'+alpha+'_DADA2_'+dada2_sim

#data_dir = '/fs/cbcb-lab/mpop/projects/SCRAPT/'
#counts_path = data_dir+'Datasets/Lupus-Microbiome-Published/Counts.dict'
#scrapt_path = data_dir+'Experiments/Lupus-Microbiome-Published/Adaptive_With_Modeshifting/sim_0.98/alpha_0.1/'
#dada2_centroid = data_dir+'Experiments/Lupus-Microbiome-Published/DADA2_Benchmarks/by_sample_dada2_vs_scrapt/dada2vsScrapt_ap/by_sample_luo.fna'
#dada2_results = data_dir+'Experiments/Lupus-Microbiome-Published/DADA2_Benchmarks/by_sample_dada2_vs_scrapt/dada2vsScrapt_ap/luo_by_sample.csv'
#dada2_scrapt = data_dir+'Experiments/Lupus-Microbiome-Published/DADA2_Benchmarks/DADA2_SCRAPT/sim_'+scrapt_sim+'_alpha_'+alpha+'_'+dada2_sim
#outdir = data_dir+'Experiments/Lupus-Microbiome-Published/DADA2_Benchmarks/Merged_DADA2_SCRAPT/SCRAPT_'+scrapt_sim+'_Alpha_'+alpha+'_DADA2_'+dada2_sim

with open(counts_path,'r') as f:
    S = f.read().replace('\n','')
counts_dict = eval(S)

def Parse_DNACLUST(filepath, cluster):
    op = {}
    counts = []
    centroids = []
    lines = open(filepath,'r').readlines()
    err_ctr = 0
    with open(filepath) as fileobject:
        for l in fileobject:
            seqs = l.rstrip().split()
            c = seqs[0]
            for s in list(set(seqs)):
                try:
                    x = op[s]
                    err_ctr += 1
                    continue
                except KeyError:
                    pass
                op[s] = cluster
                centroids.append(c)
                try:
                    counts.append(int(counts_dict[s]['Counts']))
                except KeyError:
                    counts.append(1)
            cluster += 1
    df_scrapt_clusters = pd.DataFrame(data = {'Sequence_ID':list(op.keys()), 
                                              'Cluster_ID':list(op.values()),
                                              'Centroid':centroids,
                                              'SCRAPT_Counts':counts})
    return df_scrapt_clusters, cluster

def Load_SCRAPT_Results(data_path):
    cluster = 0
    df_scrapt = pd.DataFrame()
    for s in sorted(listdir(data_path)):
        if s.startswith('Iteration'):
            print(s)
            filepath = data_path+s+'/dnaclust_mode_bait'
            df_scrapt_iter, cluster = Parse_DNACLUST(filepath, cluster)
            df_scrapt = df_scrapt.append(df_scrapt_iter, ignore_index = True)
    print(len(df_scrapt))
    df_scrapt = df_scrapt.groupby('Centroid').aggregate({'Cluster_ID':'count','SCRAPT_Counts':sum})
    df_scrapt = df_scrapt.reset_index().rename(columns = {'Centroid':'SCRAPT_Centroid'})
    df_scrapt = df_scrapt.set_index('SCRAPT_Centroid')
    return df_scrapt

def Load_DADA2_Results(seq_path, results_path):
    buf = []
    op = []
    
    with open(seq_path) as fileobject:
        for l in fileobject:
            buf.append(l)
    for i in range(0, len(buf), 2):
        seq_id = buf[i].replace("\n","").replace(">","")
        seq = buf[i+1].replace("\n","")
        op.append({'DADA2_Centroid':seq_id, 'Seq':seq})
    df_centroids = pd.DataFrame(op)
    df_centroids = df_centroids.set_index('Seq')
    df_cluster_counts = pd.read_csv(results_path, sep = ";", skiprows=1, names=['Seq', 'DADA2_Counts'], 
                                    index_col = 'Seq')
    
    df_cluster_counts = df_cluster_counts.join(df_centroids)
    return df_cluster_counts.set_index('DADA2_Centroid')[['DADA2_Counts']]

def Load_DADA2_Map(dada2_scrapt, dada2_seed):
    lines = open(dada2_scrapt).readlines()
    op = []
    for l in lines:
        seqs = l.rstrip().split()
        d = {}
        centroid = seqs[0]
        for s in seqs[1:]:
            if dada2_seed:
                op.append({'DADA2_Centroid':centroid, 'SCRAPT_Centroid':s})
            else:
                op.append({'DADA2_Centroid':s, 'SCRAPT_Centroid':centroid})
    df_map = pd.DataFrame(op)
    return df_map

df_SCRAPT = Load_SCRAPT_Results(scrapt_path)
df_DADA2 = Load_DADA2_Results(dada2_centroid, dada2_results)
df_map = Load_DADA2_Map(dada2_scrapt, False)
df_DADA2_Map = pd.merge(df_DADA2, df_map, on = 'DADA2_Centroid', how = 'outer')
counter = len(df_DADA2_Map[df_DADA2_Map['SCRAPT_Centroid'].isna()])
df_DADA2_Map.loc[df_DADA2_Map['SCRAPT_Centroid'].isna(), 'SCRAPT_Centroid'] = np.arange(-1*counter, 0)
df_DADA2_Map_sum = df_DADA2_Map.groupby('SCRAPT_Centroid').sum()[['DADA2_Counts']]
df_DADA2_Map_sum = df_DADA2_Map_sum.join(df_SCRAPT, how = 'outer')
df_DADA2_Map_sum = df_DADA2_Map_sum.fillna(-1)
df_DADA2_Map_sum.to_csv(outdir+'_sum.txt', sep = "\t")


df_DADA2_Map_max = df_DADA2_Map.groupby('SCRAPT_Centroid').max()[['DADA2_Counts']]
df_DADA2_Map_max = df_DADA2_Map_max.join(df_SCRAPT, how = 'outer')
df_DADA2_Map_max = df_DADA2_Map_max.fillna(-1)
df_DADA2_Map_max.to_csv(outdir+'_max.txt', sep = "\t")
