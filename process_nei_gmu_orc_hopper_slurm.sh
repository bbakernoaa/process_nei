#!/bin/bash -l
#SBATCH --partition=contrib          # submit to specific partition
#SBATCH --qos=qtong                  # account information
#SBATCH --job-name=proc_nei          # name the job
#SBATCH --output=process_nei-%j.out  # write stdout to named file
#SBATCH --error=process_nei-%j.err   # write stderr to named file
#SBATCH --time=0-24:00:00            # Run for max of 00 hrs, 10 mins, 00 secs
#SBATCH --nodes=1                    # Request N nodes
#SBATCH --ntasks-per-node=48         # Request n tasks per node
#SBATCH --mem=150G                   # Request nGB RAM per core

#Environment Setings
conda activate process_nei

./run_all_area_test_xargs_gmu-nei2019_merged.sh
