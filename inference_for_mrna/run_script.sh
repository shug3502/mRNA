#!/bin/bash
#$ -cwd
matlab -singleCompThread -nosplash -nojvm -nodisplay -r $1
