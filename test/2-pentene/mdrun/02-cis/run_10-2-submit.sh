base=$PWD
for rep in {0..9}
do
    cd $base/rep_$rep
    sed -i "s/JOBNAME/pentene-C-$rep/g" jobscript.0001.start.sh
    sbatch jobscript.0001.start.sh
done
