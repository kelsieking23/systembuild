pdb2gmx() {
    echo 9 1 0 0 | gmx pdb2gmx -f ../structure/2mxu_fragment.pdb -o ../structure/2mxu_processed.gro -ignh -ter
}

editconf() {
    gmx editconf -f ../structure/2mxu_processed.gro -o ../structure/2mxu_box.gro -c -d 1.5 -bt cubic -box 11 11 11
}

gen_seed() {
    seed='build'
    cd ../structure/seeds/
    mkdir ${seed}
    cd ../../scripts/
    #gmx insert-molecules -f ../structure/2mxu_box.gro -ci ../structure/2mxu_processed.gro -nmol 5 -o ../structure/seeds/${seed}/hexamer.gro -seed ${seed}
    #cp ../topol.top ../structure/seeds/${seed}/topol.top 
    # python -c'import systembuild; systembuild.update_top_insert_molecules()'
    # gmx solvate -cp ../structure/seeds/${seed}/hexamer.gro -cs spc216.gro -p ../structure/seeds/${seed}/topol.top -o ../structure/seeds/${seed}/hexamer_solvate.gro
    gmx grompp -f ../mdp/ions.mdp -c ../structure/seeds/${seed}/hexamer_solvate.gro -p ../structure/seeds/${seed}/topol.top -o ../structure/seeds/${seed}/ions.tpr
    # cd ../structure/seeds/${seed}
    # #mkdir mindist
    # cd ../../../scripts/
    echo 13 | gmx genion -s ../structure/seeds/${seed}/ions.tpr -o ../structure/seeds/${seed}/hexamer_ions.gro -p ../structure/seeds/${seed}/topol.top -pname NA -nname CL -conc 0.15
    echo 0 1 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_0_1.xvg
    echo 0 2 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_0_2.xvg
    echo 0 3 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_0_3.xvg
    echo 0 4 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_0_4.xvg
    echo 0 5 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_0_5.xvg
    echo 1 2 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_1_2.xvg
    echo 1 3 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_1_3.xvg
    echo 1 4 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_1_4.xvg
    echo 1 5 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_1_5.xvg
    echo 2 3 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_2_3.xvg
    echo 2 4 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_2_4.xvg
    echo 2 5 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_2_5.xvg
    echo 3 4 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_3_4.xvg
    echo 3 5 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_3_5.xvg
    echo 4 5 | gmx mindist -f ../structure/seeds/${seed}/hexamer_ions.gro -s ../structure/seeds/${seed}/ions.tpr -n ../ndx/2mxu_hexamer_renumber.ndx -od ../structure/seeds/${seed}/mindist/mindist_4_5.xvg
    python -c'import systembuild; systembuild.check_mindist(1.5)'

}

#pdb2gmx
gen_seed
# i=0
# while [ $i -lt 10 ]
# do
# gen_seed
# i=$[$i+1]
# done
