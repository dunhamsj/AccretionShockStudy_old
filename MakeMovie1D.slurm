#!/bin/bash
#SBATCH --mail-user=dummymmm@vanderbilt.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --time=1-00:00:00
# #SBATCH --time=0-00:30:00
#SBATCH --output=Mov1D.out
#SBATCH --job-name=Mov1D

module restore DataAnalysis

source ${HOME}/AccretionShockDataAnalysis/bin/activate

python3 MakeMovie1D.py

deactivate
