# md_sarscov2_orf3a

Necesssary run scripts and parameters to run MD simulations in publication (preprint version) https://www.biorxiv.org/content/10.1101/2025.03.20.644376v1
The two scripts job_ligand_* contain the entire MD simulation workflow for a Protein Ligand complex using parametrization with Espaloma in OpenMM
(done with the python script https://github.com/imec-int/md_sarscov2_orf3a/blob/main/openmm2gmx_protein_ligand.py) and running, multiple energy minimyation steps
followed by equlibiration, pulling and the actual MD simulation run. The output of the pulling simulation is used as starting point for running umbrella samplings
of the protein - ligand complex. This worfklow is contained at the bottom of the scripts as a series of function calls in bash:

```bash
parametrize_system ${PARAMETRIZATION}
prep_vacuum ${OUT_NAME}_${PARAMETRIZATION}
min_0_vac ${OUT_NAME}_${PARAMETRIZATION}
min_1_vac ${OUT_NAME}_${PARAMETRIZATION}
min_2_vac ${OUT_NAME}_${PARAMETRIZATION}
prep_solvated_from_vacuum ${OUT_NAME}_${PARAMETRIZATION}
min_0 ${OUT_NAME}_${PARAMETRIZATION}
min_1 ${OUT_NAME}_${PARAMETRIZATION}
min_2 ${OUT_NAME}_${PARAMETRIZATION}
min_3 ${OUT_NAME}_${PARAMETRIZATION}
min_4 ${OUT_NAME}_${PARAMETRIZATION}
npt_5 ${OUT_NAME}_${PARAMETRIZATION}
npt_6 ${OUT_NAME}_${PARAMETRIZATION}
pull_after_npt6 ${OUT_NAME}_${PARAMETRIZATION}
npt_7_prep ${OUT_NAME}_${PARAMETRIZATION}
npt_7 ${OUT_NAME}_${PARAMETRIZATION}
```

Going through this in order: 
THe first call parametrizes the system, for which a PDB file of the complex and an SDF file of the ligand with coordinates relating to the binding pocket is required.
Then the system is put into a simulation box with periodic boundary conditions but no solvent yet. First three minimization steps are executed to relax the complex in vacuum.
Then with prep_solvated_from_vacuum the equlibiration and charge neutralization is done, followed by 5 more minimization steps. Followed by preparations and a final MD run in the NPT ensemble. 
This is based on the protocoll and series of .mdp files described here https://github.com/carlocamilloni/labtools/tree/main/mdps/atomistic

In pull_after_npt6 the necessary COM pulling between ligand and protein is performed. This generates the necesarry trajectory to pick windows from for the umbrella sampling run.
To perform the umbrella sampling the jupyter notebooks https://github.com/imec-int/md_sarscov2_orf3a/blob/main/PreUmbrellaAnalysis.ipynb and https://github.com/imec-int/md_sarscov2_orf3a/blob/main/PostUmbrellaAnalysis.ipynb
are used together with the generic job scripts for running umbrella windows: [https://github.com/imec-int/md_sarscov2_orf3a/blob/main/job_ligand_bacteriochlorin](https://github.com/imec-int/md_sarscov2_orf3a/blob/main/job_ligand_sb431542_umbrella) and
https://github.com/imec-int/md_sarscov2_orf3a/blob/main/job_ligand_bacteriochlorin_umbrella .

Finally in https://github.com/imec-int/md_sarscov2_orf3a/blob/main/run_mmbsa_SB431542 an example of how to run MMPBSA with gmx_MMPBSA can be found, using the parameters defined in mmpbsa.in.
