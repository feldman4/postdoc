#!/bin/bash
#SBATCH --job-name=akta_db_update
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --output=/home/dfeldman/logs/akta_db/%j.out
#SBATCH --error=/home/dfeldman/logs/akta_db/%j.err

# abort if something doesn't work
set -e

# submit next execution
sbatch --begin="now+1hour" /home/dfeldman/s/for/akta_db_update.sh

# update database
/home/dfeldman/s/for/akta_db_export_all.sh --after="1 day ago" --cores=1