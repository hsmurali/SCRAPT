import sys
import subprocess
from os import popen, listdir, remove, mkdir 
from os.path import isfile, split, realpath, isdir,dirname
import numpy as np 
import pandas as pd
import random
import time
import sys
from typing import Tuple
import multiprocessing
from shutil import rmtree

def DE_DUPLICATE_SEQUENCES(filepath, out_path)->bool:
	counts = out_path+"Counts.txt"
	out = out_path+"deduplicated.seqs.fna"
	command = "seqkit rmdup -s -D "+counts+" "+filepath+" > "+out
	subprocess.Popen(command,shell=True).wait()
	if isfile(out) and isfile(counts):
		return True
	return False

def SAMPLE(filepath, frac_seq, output)->bool:
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
		d[seq_ids[0]] = {'Counts':counts}#, 'Dup_List':list(seq_ids[1:])}

	f = open(head+'/Counts.dict','w')
	f.write(str(d))
	f.close()
	#remove(filepath)

def Pick_New_Alpha(curr_counts, prev_counts, prev_alpha, delta = 0.008)->float:
        if (curr_counts < prev_counts):
                return prev_alpha + (1-(curr_counts/prev_counts))*delta
        return prev_alpha

def Parse_DNACLUST_output(dnaclust_op, counts)->Tuple[pd.DataFrame, list, list]:
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

def Cluster_Sequences(sampled_seq_path, dist_cutoff, out_path, counts, num_threads)->Tuple[str, pd.DataFrame]:
        head, tail = split(dirname(realpath(__file__)))
        prog_path = head+'/dnaclust/bin/dnaclust -i '
        dnaclust_call = prog_path + (sampled_seq_path) + ' -s ' + str(dist_cutoff) + ' -t '+num_threads+' > ' + out_path + '/cluster_step'
        start_time = time.time()
        subprocess.call(dnaclust_call, shell=True)
        clust_summary, singletons, nonsingletons = Parse_DNACLUST_output(out_path+'/cluster_step', counts)
        result = " %s minutes " % str(round((time.time() - start_time)/60, 2))
        return result, clust_summary

def Split_Sequences(unclustered_seq_path, split_size, out_dir):
        mkdir(out_dir)
        command = "seqkit split "+ unclustered_seq_path +" -s " +str(split_size)+ " -O "+out_dir
        subprocess.Popen(command,shell=True).wait()

def Job(cmd):
        subprocess.call(cmd, shell=True) 

def Bait_Sequences_Split_Merge(center_list, unclustered_seq_dir, unclustered_seq_path, dist_cutoff, counts, fnames, num_threads, dna_clust_out_dir):
        start_time = time.time()
        f = open(fnames[0],'w')
        lines = [c+'\n' for c in center_list]
        f.writelines(lines)
        f.close()

        flag = Extract_Sequences(unclustered_seq_path, fnames[0], 1, fnames[1:])
        mkdir(dna_clust_out_dir)

        bait_samples = listdir(unclustered_seq_dir)
        head, tail = split(dirname(realpath(__file__)))
        prog_path = head+'/dnaclust/bin/dnaclust -i '
        
        commands = []

        for i in range(len(bait_samples)):
                if bait_samples[i].endswith(".fna"):
                        seq_path = unclustered_seq_dir + bait_samples[i]
                        dna_clust_op_path = dna_clust_out_dir+'dnaclust.'+str(i)+'.out'
                        dnaclust_bait = prog_path + seq_path + ' -s ' + dist_cutoff + ' -p ' + fnames[1] + ' -r  -t '+ num_threads + ' > ' + dna_clust_op_path
                        print(dnaclust_bait)
                        commands.append(dnaclust_bait)

        pool = multiprocessing.Pool(int(num_threads))
        result = pool.map(func=Job, iterable=commands)
        pool.close()
        pool.join()
        
        d = {}
        for i in listdir(dna_clust_out_dir):
                lines = open(dna_clust_out_dir + i,'r').readlines()
                for l in lines:
                        l = l.rstrip().split()
                        if len(l) == 1:
                                try:
                                        d[l[0]] = d[l[0]] + []
                                except KeyError:
                                        d[l[0]] = []
                        else:
                                try:
                                        d[l[0]] = d[l[0]] + l[1:]
                                except KeyError:
                                        d[l[0]] = l[1:]
        o = open(fnames[2],'w')
        for k in d:
                app = " ".join(d[k])
                o.write(k + " " + app+"\n")
        o.close()

        clust_summary, singletons, nonsingletons = Parse_DNACLUST_output(fnames[2], counts)
        result = " %s minutes " % str(round((time.time() - start_time)/60, 2))
        rmtree(dna_clust_out_dir)

        return result, clust_summary, singletons, nonsingletons

def Bait_Sequences(center_list, unclustered_seq_path, dist_cutoff, counts, fnames, num_threads)->Tuple[str, pd.DataFrame, list, list]:
        ###fnames[0]: pattern name, fnames[1]: center_names fnames[2]: unclust fnames[3]: bait_op
        start_time = time.time()
        f = open(fnames[0],'w')
        lines = [c+'\n' for c in center_list]
        f.writelines(lines)
        f.close()
        flag = Extract_Sequences(unclustered_seq_path, fnames[0], 0, fnames[1:])
        if(flag):
                head, tail = split(dirname(realpath(__file__)))
                prog_path = head+'/dnaclust/bin/dnaclust -i '
                dnaclust_bait = prog_path + fnames[2] + ' -s ' + dist_cutoff + ' -p ' + fnames[1] + ' -r  -t '+ num_threads + ' > ' + fnames[3]
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

def Extract_Sequences(seq_file, pattern_file, filter_type, names)->bool:
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


def Get_Summary(df_clust_bait, df_clust_modeshift, min_cluster_size)->dict:
        d = {}
        d['Number of clusters(Baiting on dnaclust centers)'] =  len(df_clust_bait.loc[df_clust_bait['Density'] > 1])
        d['Clustered sequences(Baiting on dnaclust centers)'] = df_clust_bait.loc[df_clust_bait['Density'] > 1,'Density'].sum()
        d['Densest cluster(Baiting on dnaclust centers)'] = df_clust_bait['Density'].max()
        d['Average number of sequences per cluster(Baiting on dnaclust centers)'] = df_clust_bait.loc[df_clust_bait['Density'] > 1, 'Density'].mean()
        d['Median number of sequences per cluster(Baiting on dnaclust centers)'] = df_clust_bait.loc[df_clust_bait['Density'] > 1, 'Density'].median()
        d['Num Clusters Above Min Cluster Size(Baiting on dnaclust centers)'] = len(df_clust_bait.loc[df_clust_bait['Density'] >= min_cluster_size])
      
        if len(df_clust_modeshift) > 0:
                d['Densest cluster(After shifting the centers)'] = df_clust_modeshift['Density'].max()
                d['Average number of sequences per cluster(After shifting the centers)']  = df_clust_modeshift.loc[df_clust_modeshift['Density']>1, 'Density'].mean()
                d['Number of clusters(After shifting the centers)'] = len(df_clust_modeshift.loc[df_clust_modeshift['Density'] > 1])
                d['Clustered sequences(After shifting the centers)'] = df_clust_modeshift.loc[df_clust_modeshift['Density'] > 1,'Density'].sum()
                d['Median number of sequences per cluster(After shifting the centers)']  = df_clust_modeshift.loc[df_clust_modeshift['Density']>1, 'Density'].median()
                d['Num Clusters Above Min Cluster Size(After shifting the centers)'] = len(df_clust_modeshift.loc[df_clust_modeshift['Density'] >= min_cluster_size])
      
        return d

def SCRAPT_Iteration(unclustered_seqs, unclustered_seq_count, sampling_rate, out_path, iter_id, min_cluster_size, counts, d_cutoff, num_threads, modeshift_flag)->Tuple[dict, list]:
        #####SAMPLING
        iteration_time = time.time()
        mkdir(out_path)
        f = SAMPLE(unclustered_seqs, sampling_rate, out_path)
        if f == False:
                print("Failed to sample sequences. Check")
                sys.exit(0)
        sampled_seq = out_path+'/Sample.fasta'
        
        #####CLUSTERING
        cluster_time, cluster_summary = Cluster_Sequences(sampled_seq, str(d_cutoff), out_path, counts, num_threads)
        print(cluster_summary.head())

        ######Split_Fasta_Files
        split_size = min(int(unclustered_seq_count/int(num_threads)), 1000000)
        Split_Sequences(unclustered_seqs, split_size, out_path + '/Split_Seqs/')

        
        #####BAITING
        try:
                bait_centroids = cluster_summary.loc[((cluster_summary['Density'] > 1) | (cluster_summary['Mode'] > 1)), 'Centroid'].tolist()
                #bait_centroids = cluster_summary.loc[((cluster_summary['Density'] > 1)), 'Centroid'].tolist()
        except KeyError:
                return {}, unclustered_seqs

        op_filenames = [out_path+'/bait_centroids.txt', out_path+'/bait_centroids.fna', out_path+'/dnaclust_bait']
        (bait_time, 
         bait_summary, 
         bait_singletons, 
         bait_nonsingletons) = Bait_Sequences_Split_Merge(bait_centroids, out_path+'/Split_Seqs/', unclustered_seqs, str(d_cutoff), counts, op_filenames, num_threads, out_path+'/Bait/')

        #####MODE SHIFTING
        if modeshift_flag == True:
                print(bait_summary.head())
                try:
                        mode_centroids = bait_summary.loc[(bait_summary['Density'] > 1), 'Mode_Seq'].tolist()
                        op_filenames = [out_path+'/bait_mode.txt', out_path+'/bait_mode.fna', out_path+'/dnaclust_mode_bait']
                        (mode_shift_time, 
                         mode_shift_summary, 
                         mode_shift_singletons, 
                         mode_shift_nonsingletons) = Bait_Sequences_Split_Merge(mode_centroids, out_path+'/Split_Seqs/',unclustered_seqs, str(d_cutoff), counts, op_filenames, num_threads, out_path+'/Bait_Mode/')
                        mode_shift_summary.to_csv(out_path+'/Cluster_Summary.txt', sep = '\t')
                        clustered_seqs_list = mode_shift_nonsingletons
                except KeyError:
                        mode_shift_summary = []
                        mode_shift_singletons = []
                        clustered_seqs_list = bait_nonsingletons
                        bait_summary.to_csv(out_path+'/Cluster_Summary.txt', sep = '\t')
                        mode_shift_time = 0

        else:
                mode_shift_summary = []
                mode_shift_singletons = []
                clustered_seqs_list = bait_nonsingletons
                bait_summary.to_csv(out_path+'/Cluster_Summary.txt', sep = '\t')
                mode_shift_time = 0

        #####CREATE UNCLUSTERED SEQUENCES
        fp = open(out_path+'/unclustered_seqs','w')
        lines = [c + '\n' for c in clustered_seqs_list]
        fp.writelines(lines)
        fp.close()
        unclust_extract = Extract_Sequences(unclustered_seqs, out_path+'/unclustered_seqs', 2, [out_path+'/unclustered_seqs_'+str(iter_id)+'.fna'])
        if(unclust_extract == False):
                print('Failed to extract unclustered sequences')
                sys.exit(0)   
        iter_bm = " %s minutes " % str(round((time.time() - iteration_time)/60, 2))
        remove(sampled_seq)
        if iter_id > 1:
                remove(unclustered_seqs)

        ######SUMMARIZE CLUSTERS
        rmtree(out_path + '/Split_Seqs/')
        unclustered_seqs = out_path+'/unclustered_seqs_'+str(iter_id)+'.fna'
        if len(bait_summary) == 0:
                return {}, unclustered_seqs
                
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