#!/bin/bash
#
#SBATCH --partition=troll
#SBATCH --job-name=core+staszka
#SBATCH --output=runlog.txt

#SBATCH --ntasks=8

~/trollMPI/mpirun -np 8 ~/mb starter.nexus
