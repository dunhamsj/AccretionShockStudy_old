#!/bin/bash
#SBATCH --mail-user=samuel.j.dunham@vanderbilt.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --time=0-00:20:00
#SBATCH --output=Mov2D.out
#SBATCH --job-name=Mov2D

module restore DataAnalysis

source ${HOME}/AccretionShockDataAnalysis/bin/activate

python3 MakeMovie2D.py

deactivate
