base=$PWD
for rep in {0..9}
do
    cd $base
    mkdir rep_$rep
    cd rep_$rep
    cp ../rep_999/* ./ -r
    ./run_1_eq.sh
done
