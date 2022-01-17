#!/bin/bash
#SBATCH --mail-user=samuel.j.dunham@vanderbilt.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --time=0-00:05:00
#SBATCH --output=PlotFields.out
#SBATCH --job-name=PlotFields

module restore DataAnalysis

source ${HOME}/AccretionShockDataAnalysis/bin/activate

python3 PlotFieldsAMReX2D.py

deactivate
