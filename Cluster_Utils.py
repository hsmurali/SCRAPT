import sys
import subprocess
from os import popen, listdir, remove, mkdir 
from os.path import isfile, split, realpath, isdir
import numpy as np 
import pandas as pd
import random
import time
import sys

def DE_DUPLICATE_SEQUENCES(filepath, out_path):
	counts = out_path+"Counts.txt"
	out = out_path+"deduplicated.seqs.fna"
	command = "seqkit rmdup -s -D "+counts+" "+filepath+" > "+out
	subprocess.Popen(command,shell=True).wait()
	if isfile(out) and isfile(counts):
		return True
	return False

def SAMPLE(filepath, frac_seq, output):
	head, tail = split(filepath)
	command = 'seqkit sample -p '+str(frac_seq)+' -s ' + str(random.randint(0,100)) + ' ' +filepath +' > '+ output + '/Sample.fasta'
	subprocess.Popen(command,shell=True).wait()
	if isfile(output + '/Sample.fasta'):
		return True
	return False

def WRITE_COUNTS_DICT(filepath):
	head, tail = split(filepath)
	lines = open(filepath).readlines()
	d = {}
	for l in lines:
		a = l.split('\t')
		counts = a[0]
		seq_ids = a[1].split(', ')
		d[seq_ids[0]] = {'Counts':counts, 'Dup_List':list(seq_ids[1:])}

	f = open(head+'/Counts.dict','w')
	f.write(str(d))
	f.close()
	#remove(filepath)

def Pick_New_Alpha(curr_counts, prev_counts, prev_alpha, delta = 0.008):
        if (curr_counts < prev_counts):
                return prev_alpha + (1-(curr_counts/prev_counts))*delta
        return prev_alpha

def Parse_DNACLUST_output(dnaclust_op, counts):
        lines = open(dnaclust_op,'r').readlines()
        singletons = []
        nonsingletons = []
        cluster_id = 0
        out_list = []
        for l in lines:
                reads = l.rstrip().split()
                centroid = reads[0]
                density = len(reads)
                mode= 0
                mode_seq  = ''
                for r in reads:
                        try:
                                ctr = int(counts[r]['Counts'])
                        except KeyError:
                                ctr = 1
                        if ctr > mode:
                                mode = ctr
                                mode_seq = r
                d = {'Cluster_Id':cluster_id, 'Centroid':centroid,'Density':density,'Mode':mode, 'Mode_Seq':mode_seq}
                cluster_id += 1
                if density == 1:
                        singletons += reads
                elif density > 1:
                        nonsingletons += reads
                out_list.append(d)
        df = pd.DataFrame(out_list)
        return df, singletons, nonsingletons

def Cluster_Sequences(sampled_seq_path, dist_cutoff, out_path, counts):
        dnaclust_call = 'dnaclust/dnaclust_linux_release3/dnaclust ' + (sampled_seq_path) + ' -s ' + str(dist_cutoff) + ' > ' + out_path + '/cluster_step'
        start_time = time.time()
        subprocess.call(dnaclust_call, shell=True)
        clust_summary, singletons, nonsingletons = Parse_DNACLUST_output(out_path+'/cluster_step', counts)
        result = " %s minutes " % str(round((time.time() - start_time)/60, 2))
        return result, clust_summary

def Bait_Sequences(center_list, unclustered_seq_path, dist_cutoff, counts, fnames, num_threads):
        ###fnames[0]: pattern name, fnames[1]: center_names fnames[2]: unclust fnames[3]: bait_op
        start_time = time.time()
        f = open(fnames[0],'w')
        lines = [c+'\n' for c in center_list]
        f.writelines(lines)
        f.close()
        flag = Extract_Sequences(unclustered_seq_path, fnames[0], 0, fnames[1:])
        if(flag):
                dnaclust_bait = 'dnaclust/dnaclust_linux_release3/dnaclust ' + fnames[2] + ' -s ' + dist_cutoff + ' -p ' + fnames[1] + ' -r  -t '+ num_threads + ' > ' + fnames[3]
                subprocess.call(dnaclust_bait, shell=True)
                clust_summary, singletons, nonsingletons = Parse_DNACLUST_output(fnames[3], counts)
                result = " %s minutes " % str(round((time.time() - start_time)/60, 2))
                remove(fnames[0])
                remove(fnames[1])
                remove(fnames[2])

                return result, clust_summary, singletons, nonsingletons
        else:
                print('Failed to extract sequences... check commands')
                sys.exit(1)

def Extract_Sequences(seq_file, pattern_file, filter_type, names):
        #both = 0, found = 1, not_found = 2
        command_extract_found = "seqkit grep --pattern-file "+pattern_file +" "+seq_file
        command_extract_notfound =  "seqkit grep -v --pattern-file "+pattern_file +" "+seq_file

        if filter_type == 0:
                if len(names) == 1:
                        return False
                subprocess.call(command_extract_found + " > "+names[0],shell=True)
                subprocess.call(command_extract_notfound+" > "+names[1],shell=True)
                return True
        elif filter_type == 1:
                if len(names) == 0:
                        return False
                subprocess.call(command_extract_found + " > "+names[0],shell=True)
                return True

        elif filter_type == 2:
                if len(names) == 0:
                        return False
                subprocess.call(command_extract_notfound+" > "+names[0],shell=True)
                return True


def Get_Summary(df_clust_bait, df_clust_modeshift, min_cluster_size):
        d = {}
        d['Number of clusters(Baiting on dnaclust centers)'] =  len(df_clust_bait.loc[df_clust_bait['Density'] > 1])
        d['Number of clusters(After shifting the centers)'] = len(df_clust_modeshift.loc[df_clust_modeshift['Density'] > 1])
        d['Clustered sequences(Baiting on dnaclust centers)'] = df_clust_bait.loc[df_clust_bait['Density'] > 1,'Density'].sum()
        d['Clustered sequences(After shifting the centers)'] = df_clust_modeshift.loc[df_clust_modeshift['Density'] > 1,'Density'].sum()
        d['Densest cluster(After shifting the centers)'] = df_clust_bait['Density'].max()
        d['Densest cluster(Baiting on dnaclust centers)'] = df_clust_modeshift['Density'].max()
        d['Average number of sequences per cluster(Baiting on dnaclust centers)'] = df_clust_bait.loc[df_clust_bait['Density'] > 1, 'Density'].mean()
        d['Average number of sequences per cluster(After shifting the centers)']  = df_clust_modeshift.loc[df_clust_modeshift['Density']>1, 'Density'].mean()
        d['Median number of sequences per cluster(Baiting on dnaclust centers)'] = df_clust_bait.loc[df_clust_bait['Density'] > 1, 'Density'].median()
        d['Median number of sequences per cluster(After shifting the centers)']  = df_clust_modeshift.loc[df_clust_modeshift['Density']>1, 'Density'].median()
        d['Num_Clusters_Above_Min_Cluster_Size'] = len(df_clust_modeshift.loc[df_clust_modeshift['Density'] >= min_cluster_size])
      
        return d

def SCRAPT_Iteration(unclustered_seqs, sampling_rate, out_path, iter_id, min_cluster_size, counts, d_cutoff, num_threads):
        #####SAMPLING
        iteration_time = time.time()
        mkdir(out_path)
        f = SAMPLE(unclustered_seqs, sampling_rate, out_path)
        if f == False:
                print("Failed to sample sequences. Check")
                sys.exit(0)
        sampled_seq = out_path+'/Sample.fasta'
        
        #####CLUSTERING
        cluster_time, cluster_summary = Cluster_Sequences(sampled_seq, str(d_cutoff), out_path, counts)
        
        #####BAITING
        bait_centroids = cluster_summary.loc[((cluster_summary['Density'] > 1) | (cluster_summary['Mode'] > 1)), 'Centroid'].tolist()
        op_filenames = [out_path+'/bait_centroids.txt', out_path+'/bait_centroids.fna', out_path+'/bait_subjects.fna', out_path+'/dnaclust_bait']
        (bait_time, bait_summary, 
         bait_singletons, bait_nonsingletons) = Bait_Sequences(bait_centroids, unclustered_seqs, str(d_cutoff), counts, op_filenames, num_threads)
        
        #####MODE SHIFTING
        mode_centroids = bait_summary.loc[(bait_summary['Density'] > 1), 'Mode_Seq'].tolist()
        op_filenames = [out_path+'/bait_mode.txt', out_path+'/bait_mode.fna', out_path+'/bait_mode_subjects.fna', out_path+'/dnaclust_mode_bait']
        (mode_shift_time, mode_shift_summary, 
         mode_shift_singletons, mode_shift_nonsingletons) = Bait_Sequences(mode_centroids, unclustered_seqs, str(d_cutoff), counts, op_filenames, num_threads)
        mode_shift_summary.to_csv(out_path+'/Cluster_Summary.txt', sep = '\t')
        
        #####CREATE UNCLUSTERED SEQUENCES
        fp = open(out_path+'/unclustered_seqs','w')
        lines = [c + '\n' for c in mode_shift_nonsingletons]
        fp.writelines(lines)
        fp.close()
        unclust_extract = Extract_Sequences(unclustered_seqs, out_path+'/unclustered_seqs', 2, [out_path+'/unclustered_seqs_'+str(iter_id)])
        if(unclust_extract == False):
                print('Failed to extract unclustered sequences')
                sys.exit(0)   
        iter_bm = " %s minutes " % str(round((time.time() - iteration_time)/60, 2))
        remove(sampled_seq)
        if iter_id > 1:
                remove(unclustered_seqs)

        ######SUMMARIZE CLUSTERS
        unclustered_seqs = out_path+'/unclustered_seqs_'+str(iter_id)
        d = Get_Summary(bait_summary, mode_shift_summary, min_cluster_size)
        d['Sampling Rate'] = sampling_rate
        d['Cluster Singletons'] = len(cluster_summary[cluster_summary['Density'] == 1])
        d['Singletons(Baiting on dnaclust centers)'] = len(bait_singletons)
        d['Singletons(After shifting the centers)'] = len(mode_shift_singletons)
        d['Time(Cluster)'] = cluster_time
        d['Time(Bait)'] = bait_time
        d['Time(Mode shift)'] =  mode_shift_time
        d['Time(Total)'] = iter_bm
        d['Iteration'] = iter_id

        print(iter_bm)
        return d, unclustered_seqs