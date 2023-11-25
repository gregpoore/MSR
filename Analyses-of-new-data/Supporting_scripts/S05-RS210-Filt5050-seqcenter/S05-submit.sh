#!/bin/bash
#SBATCH --job-name=mlFilt5050  # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gregory.d.poore@gmail.com   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=120gb                   # Job Memory
#SBATCH --time=36:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlFilt5050_%j.log    # Standard output and error log
#SBATCH --reservation=dtmcdonald_14

work_dir=/home/gdpoore/projects/tcga/AA_Rebuttal2023b/S04-RS210-Filt5050-seqcenter
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S05-RS210-Filt5050-seqcenter.R

