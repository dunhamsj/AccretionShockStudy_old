#!/bin/bash
#SBATCH --mail-user=samuel.j.dunham@vanderbilt.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=240G
#SBATCH --time=1-00:00:00
#SBATCH --output=RelaxM.out
#SBATCH --job-name=RelaxM

module restore DataAnalysis

source ${HOME}/AccretionShockDataAnalysis/bin/activate

python3 Relaxation_MultipleRuns.py

deactivate
