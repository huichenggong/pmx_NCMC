# pmx_NCMC
Add replica exchange [(ref1)](#4-reference) to pmx free energy calculation. The more general form of
replica exchange is nonequilibrium candidate Monte Carlo (NCMC) [(ref2)](#4-reference). This python 
package is an IO based implementation of replica exchange for pmx-style free 
energy calculation.

## 1. Installation
### 1.1 Create a conda/mamba env
```bash
mamba create -n pmx_NCMC python=3.11 pymbar numpy pandas matplotlib -c conda-forge
mamba activate  pmx_NCMC  # mambda works the same as conda
```
### 1.2 Pip install package
```
git clone https://github.com/huichenggong/pmx_NCMC.git
pip install .
```
### 1.X Remove environment
```bash
mamba remove -n pmx_NCMC --all
```


## 2. Example
```bash
cd test/2-pentene/
```
In stateA, it's a normal pentene. In stateB, the dihedral potential around 
the double bond is removed. We can either start from trans or cis.
```bash
mdrun/01-trans/em.gro
mdrun/02-cis/em.gro
```
Eventrually, the free energy difference between stateA and stateB should be the same, 
regardless of the starting conformation. 
```bash
cd mdrun/01-trans/
./run_10-1-prepare.sh # prepare 10 replicas
./run_10-2-submit.sh  # submit 10 replicas to the cluster
```
In case you don't want to submit them to the cluster.
```bash
cd mdrun/01-trans/rep_999
./run_1_eq.sh # prepare this replica
pmx_mdrun -h
pmx_mdrun \
    -mdp_folder mdp/ -p ../../../topol.top -folder_start 000000 \
    -cycle 20 \
    -MDRUN "mpirun -np 2 --bind-to none gmx_mpi mdrun" \
    -GROMPP "gmx grompp" \
    -maxh 5 # 20 is too little. You probably need 200 cycles to converge the delta G
pmx_mdrun \
    -mdp_folder mdp/ -p ../../../topol.top -folder_start 000019 \
    -cycle 20 \
    -MDRUN "mpirun -np 2 --bind-to none gmx_mpi mdrun" \
    -GROMPP "gmx grompp" \
    -maxh 5 # Append 20 more cycles from 000019
```

## 3. Theory
![Theory](./Fig/theory.jpg)

## 4. Reference
(1)	Ballard, A. J.; Jarzynski, C. Replica Exchange with Nonequilibrium Switches. Proc. Natl. Acad. Sci. 2009, 106 (30), 12224–12229. https://doi.org/10.1073/pnas.0900406106.  
(2)	Nilmeier, J. P.; Crooks, G. E.; Minh, D. D. L.; Chodera, J. D. Nonequilibrium Candidate Monte Carlo Is an Efficient Tool for Equilibrium Simulation. Proc. Natl. Acad. Sci. 2011, 108 (45). https://doi.org/10.1073/pnas.1106094108.  
