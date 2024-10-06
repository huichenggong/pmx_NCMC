base=$PWD

for i in 0 1 
do
    cd $base/eq/$i
    gmx grompp -f eq.mdp -c ../../../em.gro -p ../../../../../topol.top  -o eq
done

cd $base
mpirun -np 2 --bind-to none gmx_mpi mdrun -v -s eq.tpr -deffnm eq -multidir eq/? 

mkdir 000000
cd 000000
mkdir 0 1
cd $base/000000/0
gmx grompp -f ../../mdp/eq0.mdp -c ../../eq/0/eq.gro -p ../../../../../topol.top -o eq
cd $base/000000/1
gmx grompp -f ../../mdp/eq1.mdp -c ../../eq/1/eq.gro -p ../../../../../topol.top -o eq

cd $base

