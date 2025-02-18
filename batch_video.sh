#!/bin/bash

# Request resources:
#SBATCH -c 1     # 1 entire node
#SBATCH --time=12:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=1G      # 1 GB RAM
#SBATCH -p shared


module load python
module load ffmpeg
pip install pandas

python make_video.py
