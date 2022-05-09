import sys
from Cluster_Utils import *
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
                d, unclustered_seqs = SCRAPT_Iteration(unclustered_seqs, unclustered_seq_count, alpha_new, out_path, ctr, min_cluster_size, counts, d_cutoff, num_threads, modeshift)
                prev_counts = curr_counts

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
                summary_list.append(d)
        try:
                df_summary = pd.DataFrame(summary_list)
                df_summary.to_csv(directory+'/Summary.txt', sep = '\t')
        except FileNotfoundError:
                print('Failed to create summary')

