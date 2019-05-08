for X in $1; do
    cd $X
    sbatch $X/initRun.sh
done