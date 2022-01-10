#!/bin/bash
#SBATCH --mail-user=samuel.j.dunham@vanderbilt.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --time=0-02:00:00
#SBATCH --output=Relax.out
#SBATCH --job-name=Relax

module restore DataAnalysis

source ${HOME}/AccretionShockDataAnalysis/bin/activate

python3 Relaxation.py

deactivate
