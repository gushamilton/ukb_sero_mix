#!/bin/bash
#SBATCH --job-name=aim3_scenario_35
#SBATCH --output=logs/aim3_scenario_35_%%j.log
#SBATCH --error=logs/aim3_scenario_35_%%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=16
#SBATCH --account=sscm013902

# Load R module (adjust as needed for your cluster)
module load languages/R

# Set working directory to the parallel folder
cd /user/work/fh6520/sim_mixture

# Run the scenario
echo "Starting scenario 35 at \\$(date)"
Rscript test.R 
echo "Completed scenario 35 at \\$(date)"
