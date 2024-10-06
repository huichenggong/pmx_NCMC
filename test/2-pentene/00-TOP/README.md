# 1. Prepare mol.itp
```bash
antechamber -i pentene.pdb -fi pdb  -o mol.mol2 -fo mol2 -c bcc -nc 0
parmchk2 -i mol.mol2 -f mol2 -o frcmod
tleap -f tleap.in
acpype -x mol.inpcrd -p mol.prmtop 
```
Edit the `MOL.amb2gmx/MOL_GMX.top` file. Split the `[ moleculetype ]` into a new file called `mol.itp`
# 2. Remove C2-C3 dihedral in stateB
```
cp mol.itp mol-hybrid.itp
```
Add stateB to the dihedrals which are `X1-C2-C3-X2` in file `mol-hybrid.itp`

# 3. Prepare `../topol.top`
copy a templete from another amber top file.  
copy past `[ atomtypes ]` from `MOL.amb2gmx/MOL_GMX.top` to `at_type.itp` 

# 4. Solvate
```
gmx editconf -f pentene.pdb -o 01-box.pdb -bt cubic -c -d 1.1
gmx solvate -p ../topol.top -o 02-solv.pdb -cp 01-box.pdb -cs spc216.gro
```

# 5. Rotate
`03-cis.pdb`  
