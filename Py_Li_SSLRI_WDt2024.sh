#!/bin/bash
#SBATCH --job-name Py_li_WDt2024
#SBATCH --output Py_li_WDt2024_output.log
#SBATCH --error Py_li_WDt2024_error.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 75G
 
module load python/3.11.6

python LIANA_SSLRI_WDt2024.py
