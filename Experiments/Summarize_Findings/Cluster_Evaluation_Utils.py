import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from os import listdir, mkdir
from os.path import isdir
from functools import partial


def Compile_Cluster_Summaries(scrapt_cluster_directory):
	iteration_results = listdir(scrapt_cluster_directory)
	df_clusters = pd.DataFrame()
	df_cluster_summary = pd.DataFrame()
	summary = []

	for i in iteration_results:
		if i.startswith('Iteration_'):
			df = pd.read_csv(scrapt_cluster_directory+i+'/Cluster_Summary.txt', sep = '\t')
			del df['Unnamed: 0']
			df['Iteration'] = int(i.replace("Iteration_",""))+1
			clustered_seq_count = df['Density'].sum()
			sig_cluster_10 = len(df[df['Density']>=10])
			sig_cluster_20 = len(df[df['Density']>=20])
			sig_cluster_30 = len(df[df['Density']>=30])
			sig_cluster_35 = len(df[df['Density']>=35])
			sig_cluster_40 = len(df[df['Density']>=40])

			sig_cluster_10_counts = df[df['Density']>=10]['Density'].sum()
			sig_cluster_20_counts = df[df['Density']>=20]['Density'].sum()
			sig_cluster_30_counts = df[df['Density']>=30]['Density'].sum()
			sig_cluster_40_counts = df[df['Density']>=40]['Density'].sum()
			sig_cluster_50_counts = df[df['Density']>=50]['Density'].sum()
			

			summary.append({'Iteration':int(i.replace("Iteration_",""))+1, 
							'Seq_Counts':clustered_seq_count, 
							'Clusters_Above_10':sig_cluster_10,
							'Clusters_Above_20':sig_cluster_30,
							'Clusters_Above_30':sig_cluster_30,
							'Clusters_Above_35':sig_cluster_35,
							'Clusters_Above_40':sig_cluster_40, 
							'Seqs_in_Cluster_Above_10':sig_cluster_10_counts,
							'Seqs_in_Cluster_Above_20':sig_cluster_20_counts,
							'Seqs_in_Cluster_Above_30':sig_cluster_30_counts,
							'Seqs_in_Cluster_Above_40':sig_cluster_40_counts,
							'Seqs_in_Cluster_Above_50':sig_cluster_50_counts,
							})
			
			df_clusters = df_clusters.append(df, ignore_index = True)
	df_cluster_summary = pd.DataFrame(summary)

	df_SCRAPT_summary = pd.read_csv(scrapt_cluster_directory+'Summary.txt',sep = "\t")
	del df_SCRAPT_summary['Unnamed: 0']
	try:
		df_SCRAPT_summary['Time(Cluster)'] = df_SCRAPT_summary['Time(Cluster)'].str.replace(" minutes","").astype(float)
	except AttributeError:
		pass
	try:
		df_SCRAPT_summary['Time(Bait)'] = df_SCRAPT_summary['Time(Bait)'].str.replace(" minutes","").astype(float)
	except AttributeError:
		pass
	try:	
		df_SCRAPT_summary['Time(Mode shift)'] = df_SCRAPT_summary['Time(Mode shift)'].str.replace(" minutes","").astype(float)
	except AttributeError:
		pass
	try:	
		df_SCRAPT_summary['Time(Total)'] = df_SCRAPT_summary['Time(Total)'].str.replace(" minutes","").astype(float)
	except AttributeError:
		pass

	df_SCRAPT_summary = df_SCRAPT_summary.set_index('Iteration')
	df_cluster_summary = df_cluster_summary.set_index('Iteration')
	df_SCRAPT_summary = df_SCRAPT_summary.join(df_cluster_summary)

	return df_clusters, df_SCRAPT_summary

def Compute_Fragmentation_Measure(cluster_list, clust_thresh=0): 
	###Add a normalization factor
	x = np.arange(0, 100, 0.5)
	cluster_list = np.sort(cluster_list)[::-1]
	cluster_list = cluster_list[cluster_list > clust_thresh]
	temp = np.cumsum(cluster_list)
	NG, cum = [], []
	norm = np.sum(cluster_list)
	for i in x:
		try:
			NG.append(round(cluster_list[np.where(temp <= norm*i/100)[0][-1]],2))
			cum.append(temp[np.where(temp <= norm*i/100)[0][-1]])
		except IndexError:
			pass
	return cum, NG

def Parse_DNACLUST_Outputs(filepath):
	dnaclust_counts = open(filepath,'r').readlines()
	dnaclust_size = []
	with open(filepath) as fileobject:
		for d in fileobject:
			d = d.rstrip()
			if len(d.split()) >= 2:
				dnaclust_size.append(len(d.split()))
	return dnaclust_size

def Parse_CDHIT_Outputs(filepath):
	cdhit_clusters = []
	with open(filepath) as file:
		counter = 0
		centroid = []
		density = []
		l = 0
		for line in file:
			a = line.rstrip()
			if a[0] == ">":
				if(l!= 0 ):
					density.append(counter)       
				centroid.append(a[1:])
				counter = 0
			else:
				counter = counter + 1 
			l = l + 1
		density.append(counter)
	return density