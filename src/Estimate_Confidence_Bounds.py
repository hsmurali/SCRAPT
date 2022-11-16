import numpy as np
import pandas as pd

def Enumerate(v):
    e_v = []
    for c in range(1, len(v[1:])+1):
        try: e_v += [c]*int(v[c])
        except ValueError: continue
    e_v = np.array(e_v)
    return e_v

def Return_Counts(v, normalize = 1):
    D = np.bincount(v)
    if normalize:
        return D/D.sum()
    return D

def SCRAPT_Simulation(observed_clusters, N, alpha, D, num_realizations=10000):
    obs_counts = Return_Counts(observed_clusters, 0)
    Freq_Mat, Dist = [], []

    V = np.arange(0, len(D))
    V[0] = 1
    cluster_counts = D*N/V
    
    cluster_counts[:len(obs_counts)] = np.maximum(obs_counts,cluster_counts[:len(obs_counts)]) 
    clusters = Enumerate(cluster_counts)
    n = int(alpha*N)
    
    P =  np.exp(-1*clusters*n/N)*(1+clusters*n/N*np.exp(clusters/N))
    for i in range(num_realizations):
        Freq_Counts = np.zeros(len(D)+1)
        R = np.random.uniform(0, 1, len(clusters))
        Clusters_Discovered = clusters[np.where(R >= P)[0]]
        np.add.at(Freq_Counts, Clusters_Discovered, 1)
        Freq_Mat.append(Freq_Counts)
    
    Freq_Mat = np.array(Freq_Mat).T.astype(int)
    for j in range(Freq_Mat.shape[0]):
        d = Return_Counts(Freq_Mat[j])
        Dist.append(d)
    return Dist           
 