import sys
from Cluster_Utils import *
from Estimate_Confidence_Bounds import *
from copy import deepcopy
import argparse as ap
import ast

if __name__ == '__main__':
        parser = ap.ArgumentParser(description="SCRAPT: Sampling Clustering Recruiting AdaPt and iTerate."+
                                  "SCRAPT is a tool to cluster phylogenetic marker gene sequences, using an iterative approach.")
        requiredNamed = parser.add_argument_group('required named arguments')
        optionalNamed = parser.add_argument_group('optional named arguments')
        
        requiredNamed.add_argument("-f", "--filepath", help="Path to the file containing 16S sequences. At the moment we support only fatsa file containing the 16S reads.", required=True)
        requiredNamed.add_argument("-o","--output_directory", help="Location to write the outputs to", required = True)
        
        optionalNamed.add_argument("-s","--Sampling_Rate", help="Initial Sampling Rate (between 0 and 100). [DEFAULT = 0.1]", required=False, default = "0.1")
        optionalNamed.add_argument("-a","--adaptive", help="Flag to run SCRAPT in adaptve mode. [DEFAULT = True]", required = False, default = "True")
        optionalNamed.add_argument("-d","--delta", help = "Adjustment constant for the adaptive sampling. [DEFAULT = 0.008]", required = False, default = "0.008")
        optionalNamed.add_argument("-r","--similarity", help="Similairity threshold for clustering. [DEFAULT = 0.99]", required = False, default = "0.99")
        optionalNamed.add_argument("-n","--max_iterations", help="Maximum number of iterations to run the iterative clustering. [DEFAULT = 50]", required = False, default = "50")
        optionalNamed.add_argument("-k","--min_cluster", help="Size of the smallest cluster to detect. [DEFAULT = 50]", required = False, default = "50")
        optionalNamed.add_argument("-t","--num_threads", help="Number of threads. [DEFAULT = 8]", required = False, default="8")
        optionalNamed.add_argument("-c","--counts_dict", help="Path to the dictionary of counts. If it is not provided SCRAPT deduplicates the sequences.", required = False, default = "")
        optionalNamed.add_argument("-m","--mode_shift", help="Perform Modeshifting. [DEFAULT = True]", required = False, default = "True")
        optionalNamed.add_argument("-p","--prob", help="Confidence threshold. [DEFAULT = 0.9]", required = False, default = "0.9")
        optionalNamed.add_argument("-b","--realizations", help="Number of realizations to estimate confidence bounds. [DEFAULT = 1000]", required = False, default = "1000")
        
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
        p = float(args.prob)
        t = int(args.realizations)
        cum_ctr, cum_clstr = 0, 0
        counts_path = args.counts_dict
        if args.mode_shift.lower().lstrip().rstrip() == "true":
                modeshift = True
        elif args.mode_shift.lower().lstrip().rstrip() == "false":
                modeshift = False
        else:
                print("Invalid argument to mode_shift")
                sys.exit(0)

        if not(isdir(directory)):
                mkdir(directory)

        f = True
        if counts_path == "":
                f = DE_DUPLICATE_SEQUENCES(seq_path, directory)
                f = WRITE_COUNTS_DICT(directory+"Counts.txt")
                print("Duplicated Sequences")
                counts_path = directory+"Counts.dict"
                seq_path = directory+"deduplicated.seqs.fna"

        if(f == False):
                print("Failed to Deduplicate sequences. Check")
                sys.exit(0)

        start_time = time.time()
        unclustered_seqs = seq_path
        counts = {}  
        with open(counts_path) as multiplicities:
                counts = ast.literal_eval(multiplicities.read())
        pre_processing_time = " %s minutes " % str(round((time.time() - start_time)/60, 2))
        summary_list = []
        ctr = 0
        cmd_read_counts = "grep \"^>\" -c "+seq_path
        unclustered_seq_count = int(subprocess.check_output(cmd_read_counts, shell = True))
        seq_count = deepcopy(unclustered_seq_count)

        clusters = np.array([])
        clusters_not_frozen = np.array([])
        Frozen_frontier = np.zeros(num_iterations+1)

        while (ctr < num_iterations):
                if adaptive == "True":
                        if(ctr >= 2): alpha_new = Pick_New_Alpha(curr_counts, prev_counts, alpha_new, delta)
                        else: alpha_new = start_alpha                
                else: alpha_new = start_alpha
                print(ctr, '------->', prev_counts, curr_counts, alpha_new, cum_ctr, cum_clstr)
                out_path = directory+'/Iteration_' + str(ctr)
                ctr += 1
                d, unclustered_seqs = SCRAPT_Iteration(unclustered_seqs, unclustered_seq_count, alpha_new, out_path, ctr, min_cluster_size, counts, d_cutoff, num_threads, modeshift)
                prev_counts = curr_counts
                try:
                        if (modeshift):
                                unclustered_seq_count = unclustered_seq_count - d['Clustered sequences(After shifting the centers)']
                                curr_counts =  d['Clustered sequences(After shifting the centers)'] 
                                cum_ctr +=  d['Clustered sequences(After shifting the centers)']
                                cum_clstr += d['Num Clusters Above Min Cluster Size(After shifting the centers)']
                        else:
                                unclustered_seq_count = unclustered_seq_count - d['Clustered sequences(Baiting on dnaclust centers)'] 
                                curr_counts =  d['Clustered sequences(Baiting on dnaclust centers)'] 
                                cum_ctr +=  d['Clustered sequences(Baiting on dnaclust centers)']
                                cum_clstr += d['Num Clusters Above Min Cluster Size(Baiting on dnaclust centers)']
                except KeyError:
                        unclustered_seq_count, curr_counts, cum_clstr, cum_ctr = unclustered_seq_count, 0, cum_clstr, cum_ctr
                
                df_clusters_i = pd.read_csv(out_path+'/Cluster_Summary.txt', sep = "\t")
                clusters_i = df_clusters_i['Density'].tolist()
                clusters_not_frozen = np.append(clusters_not_frozen, clusters_i).astype(int)

                cluster_counts = Return_Counts(clusters_not_frozen)
                Dist_Mat = SCRAPT_Simulation(clusters_not_frozen, seq_count, alpha_new, cluster_counts, t)
                frozen = []
                Temp = np.bincount(clusters_not_frozen)
                for j in range(1, len(Temp)):
                        counts_tmp = int(Temp[j])
                        if (len(Dist_Mat[j]) > counts_tmp) and (counts_tmp > 0) and (sum(Dist_Mat[j][counts_tmp:]) >= p):
                                frozen.append(j)
                frozen = np.array(frozen)
                
                if len(frozen) > 0:
                        not_frozen = clusters_not_frozen[np.in1d(clusters_not_frozen, frozen, invert = True)]
                        frozen = frozen[frozen >= not_frozen.max()]
                        clusters_frozen = clusters_not_frozen[np.in1d(clusters_not_frozen, frozen)]
                        clusters_not_frozen = not_frozen
                        seq_count = seq_count - sum(clusters_frozen)
                try:
                        C = np.array(clusters_i)
                        Frozen_frontier[ctr] = min(frozen)
                        d['Largest_Cluster_UB'] = Frozen_frontier[ctr]
                except ValueError:
                        Frozen_frontier[ctr] = Frozen_frontier[ctr-1]
                        d['Largest_Cluster_UB'] = Frozen_frontier[ctr]
                summary_list.append(d)

        try:
                df_summary = pd.DataFrame(summary_list)
                df_summary.to_csv(directory+'/Summary.txt', sep = '\t')
        except FileNotfoundError:
                print('Failed to create summary')

