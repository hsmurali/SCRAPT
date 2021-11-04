import sys
from Cluster_Utils import *
import argparse as ap
import ast

if __name__ == '__main__':
        parser = ap.ArgumentParser(description="SCRAPT: Sampling Clustering Recruiting AdaPt and iTerate."+
                                  "SCRAPT is a tool to cluster 16S genen sequences, using an iterative approach")
        parser.add_argument("-f", "--filepath", help="Path to the file containing 16S sequences. At the moment we support only fatsa file containing the 16S reads.", required=True)
        parser.add_argument("-o","--output_directory", help="Location to write the outputs to", required = True)
        parser.add_argument("-s","--Sampling_Rate", help="Initial Sampling Rate (between 0 and 100)", required=False, default = "0.1")
        parser.add_argument("-a","--adaptive", help="Flag to run SCRAPT in adaptve mode", required = False, default = "True")
        parser.add_argument("-d","--delta", help = "Adjustment constant for the adaptive sampling", required = False, default = "0.008")
        parser.add_argument("-r","--similarity", help="Similairity to run clustering with", required = False, default = "0.99")
        parser.add_argument("-m","--max_iterations", help="Maximum number of iterations to run the iterative clustering.", required = False, default = "50")
        parser.add_argument("-k","--min_cluster", help="Size of the smallest cluster to detect", required = False, default = "50")
        parser.add_argument("-t","--num_threads", help="Number of threads", required = False, default="8")
        args = parser.parse_args()

        seq_path = args.filepath
        directory = args.output_directory
        sampling_rate = float(args.Sampling_Rate)
        adaptive = args.adaptive
        delta = float(args.delta)
        num_iterations = int(args.max_iterations)
        d_cutoff = float(args.similarity)
        min_cluster_size = int(args.min_cluster)
        num_threads = args.num_threads
        prev_counts, curr_counts = -1, -1
        start_alpha = float(sampling_rate)/100.0
        cum_ctr, cum_clstr = 0, 0

        if not(isdir(directory)):
                mkdir(directory)

        f = True#DE_DUPLICATE_SEQUENCES(seq_path, out_path)
        #WRITE_COUNTS_DICT(out_path+"Counts.txt")

        if(f == False):
                print("Failed to Deduplicate sequences. Check")
                sys.exit(0)

        start_time = time.time()
        unclustered_seqs = seq_path
        counts = {} 
        with open('Counts.dict') as multiplicities:
                counts = ast.literal_eval(multiplicities.read())
        pre_processing_time = " %s minutes " % str(round((time.time() - start_time)/60, 2))
        summary_list = []
        ctr = 0

        while (ctr < num_iterations):
                if adaptive == "True":
                        if(ctr >= 2):
                                alpha_new = Pick_New_Alpha(curr_counts, prev_counts, alpha_new, delta)
                        else:
                                alpha_new = start_alpha                
                else:
                        alpha_new = start_alpha
                print(ctr, '------->', prev_counts, curr_counts, alpha_new, cum_ctr, cum_clstr)
                out_path = directory+'/Iteration_' + str(ctr)
                ctr += 1
                d, unclustered_seqs = SCRAPT_Iteration(unclustered_seqs, alpha_new, out_path, ctr, min_cluster_size, counts, d_cutoff, num_threads)
                prev_counts = curr_counts
                curr_counts =  d['Clustered sequences(After shifting the centers)'] 
                cum_ctr +=  d['Clustered sequences(After shifting the centers)']
                cum_clstr += d['Num_Clusters_Above_Min_Cluster_Size']
                summary_list.append(d)
        try:
                df_summary = pd.DataFrame(summary_list)
                df_summary.to_csv(directory+'/Summary.txt', sep = '\t')
        except FileNotfoundError:
                print('Failed to create summary')

