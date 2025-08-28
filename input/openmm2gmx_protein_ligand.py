from pdbfixer import PDBFixer
from openmm.app import PDBFile
import mdtraj as md
import os
import numpy as np
from rdkit import Chem
import mdtraj as md
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
import openmm.unit as unit
import openmm.app as app
from openmm import MonteCarloBarostat
from openmm import XmlSerializer
# ForceField stuff
from openmmforcefields.generators import GAFFTemplateGenerator
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmmforcefields.generators import EspalomaTemplateGenerator
from openmm.app import ForceField
from openmm import app, unit, LangevinIntegrator, Vec3
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter
import time
import sys
import argparse
import parmed as pmd



espaloma_generator = lambda lig: EspalomaTemplateGenerator(molecules=lig, forcefield='espaloma-0.3.2')
gaff_generator = lambda lig: GAFFTemplateGenerator(molecules=lig)
smirnoff_generator = lambda lig: SMIRNOFFTemplateGenerator(molecules=lig)

class MDConfig(object):
    def __init__(self, 
                simulation_name,
                #forcefield_files = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml'], 
                water_model = 'tip3p', 
                solvent_padding = 9.0 * unit.angstroms, 
                #ionic_strength = 0.15 * unit.molar, 
                ionic_strength = 0 * unit.molar, 
                constraints = app.HBonds, 
                hmass = 3.0 * unit.amu, 
                timestep = 2 * unit.femtoseconds,
                temperature = 300.0 * unit.kelvin, 
                friction_coeff = 1,
                pme_tol = 2.5e-04, 
                pdb_out = "output.pdb",
                log_file_name = "md_log.txt",
                log_report_rate = 1000, 
                traj_report_rate = 1000,
                enforcePeriodicBox = False,
                min_out = "system_minimized.pdb",
                nvt_steps = 10000,
                npt_steps = 10000,
                pressure = 1.0 * unit.atmosphere, 
                nonbonded_method = app.PME, 
                barostat_period = 50, 
                trajectory_out_name = 'trajectory.dcd',
                num_steps = 500000 # 500000*2fs = 1ns
                ):
        #self.forcefield_files = forcefield_files
        self.simulation_name = simulation_name
        self.water_model = water_model
        self.solvent_padding = solvent_padding
        self.ionic_strength = ionic_strength
        self.constraints = constraints
        self.hmass = hmass
        # integrator params
        self.timestep = timestep
        self.temperature = temperature
        self.friction_coeff = friction_coeff / unit.picosecond        
        self.pme_tol = pme_tol
        self.ligand_forcefield_generator = espaloma_generator
        # log config
        self.pdb_out = simulation_name + "_" + pdb_out
        self.log_file_name = simulation_name + "_" + log_file_name 
        self.log_report_rate = log_report_rate
        self.traj_report_rate = traj_report_rate
        self.enforcePeriodicBox = enforcePeriodicBox
        # min config
        self.min_out = simulation_name + "_" + min_out
        #self.min_steps = min_steps
        # nvt config
        self.nvt_steps = nvt_steps
        # npt config
        self.npt_steps = npt_steps
        self.pressure = pressure
        self.nonbonded_method = nonbonded_method
        self.barostat_period = barostat_period
        # prod config
        self.trajectory_out_name = simulation_name + "_" + trajectory_out_name
        self.num_steps = num_steps

class SimulationEnvironment(object):
    def __init__(self, 
                 config, simulation, solvated_system, 
                 modeller, forcefield, protein_pdb, 
                 lig_mol, protein_topology, ligand_topology,
                 context
                ):
        self.config = config
        self.simulation=simulation
        self.solvated_system=solvated_system
        self.modeller=modeller
        self.forcefield=forcefield
        self.protein_pdb=protein_pdb
        self.lig_mol=lig_mol
        self.protein_topology=protein_topology
        self.ligand_topology=ligand_topology
        self.context = context


def fix_protein(protein_path):
    fixer = PDBFixer(protein_path)
    fixer.findMissingResidues()
    
    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    missingResidues = dict()
    for key in keys:
        chain = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
            missingResidues[key] = fixer.missingResidues[key]
    fixer.missingResidues = missingResidues
    
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()

    return fixer

def createProteinLigandOnlySystem(protein_path, ligand_path, config):
    """
        create system given fixed protein structure and ligand structure
    """
    # load protein
    protein_pdb = PDBFile(protein_path)
    #protein_pdb = pmd.load_file(protein_path)
    protein_topology = protein_pdb.topology
    # load ligand
    lig_mol = Molecule.from_file(ligand_path)
    ligand_positions = lig_mol.conformers[0]
    # ugly fix to get from angstrom to nm
    ligand_positions = ligand_positions * 0.1 * unit.nanometers
    #ligand_positions = ligand_positions.in_units_of(unit.nanometers)
    ligand_topology = lig_mol.to_topology().to_openmm()
    # convert to mdtraj topology
    protein_md_topology = md.Topology.from_openmm(protein_topology)
    ligand_md_topology = md.Topology.from_openmm(ligand_topology)
    # merge topology
    complex_md_topology = protein_md_topology.join(ligand_md_topology)
    complex_topology = complex_md_topology.to_openmm()
    # get number of atoms
    n_atoms_total = complex_topology.getNumAtoms()
    n_atoms_protein = protein_topology.getNumAtoms()
    n_atoms_ligand = ligand_topology.getNumAtoms()
    
    assert n_atoms_total == n_atoms_protein + n_atoms_ligand, "Number of atoms after merging the protein and ligand topology does not match"
    #_logger.info(f"Total atoms: {n_atoms_total} (protein: {n_atoms_protein}, ligand: {n_atoms_ligand}")
    # complex positons
    complex_positions = unit.Quantity(np.zeros([n_atoms_total, 3]), unit=unit.nanometers)
    complex_positions[:n_atoms_protein, :] = protein_pdb.positions
    complex_positions[n_atoms_protein:n_atoms_protein+n_atoms_ligand, :] = ligand_positions
    # create force field
    ligand_forcefield = config.ligand_forcefield_generator(lig_mol)
    forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
    forcefield.registerTemplateGenerator(ligand_forcefield.generator)
    # solvation
    modeller = app.Modeller(complex_topology, complex_positions)
    modeller.addSolvent(
        forcefield, model=config.water_model, 
        padding=config.solvent_padding, 
        ionicStrength=config.ionic_strength
    )
    #modeller.deleteWater()
    # delete ions here to ?
    #
    #inds_ions = [a for a in solvated_topology.atoms() if a.name == "Na"]
    #for i, x in enumerate(solvated_topology.atoms()):
    #    if "NA" in str(x):
    #        print(f' atom {i} :: {x}')
    #        #inds_ions.append(i)
    #
    #print(inds_ions)
    ## have deleted waters and now ions 
    #modeller.delete(inds_ions)
    #
    # Get topology and position
    solvated_topology = modeller.getTopology()
    solvated_positions = modeller.getPositions()
    # Create system
    solvated_system = forcefield.createSystem(
        solvated_topology,
        removeCMMotion = True, 
        ewaldErrorTolerance = config.pme_tol, 
        constraints = config.constraints, 
        rigidWater = False, 
        hydrogenMass = config.hmass,
        nonbondedMethod = config.nonbonded_method
    )
    """
    integrator = LangevinIntegrator(
        config.temperature, 
        config.friction_coeff,
        config.timestep
    )

    simulation = Simulation(modeller.topology, solvated_system, integrator)
    context = simulation.context
    context.setPositions(modeller.positions)
    
    simulation.reporters.append(
        PDBReporter(
            config.pdb_out, 
            config.log_report_rate, 
            enforcePeriodicBox=config.enforcePeriodicBox
        )
    )
    simulation.reporters.append(
        StateDataReporter(
            config.log_file_name, config.log_report_rate, step=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, speed=True,
            density=True,
        )
    )
    """
    return SimulationEnvironment(
        config, None, solvated_system, 
        modeller, None, protein_pdb, 
        lig_mol, protein_topology, ligand_topology,
        None
    )
    
def to_gromacs_system(system, outname):
    parmed_structure = pmd.openmm.load_topology(
        system.modeller.topology, 
        system=system.solvated_system, 
        xyz=system.modeller.positions
    )
    # adding default bond types where bond type is none
    bond_type = pmd.BondType(1.0, 1.0, list=parmed_structure.bond_types)
    parmed_structure.bond_types.append(bond_type)
    for bond in parmed_structure.bonds:
        if bond.type is None:
            bond.type = bond_type
    parmed_structure.save(outname + '.top', overwrite=True)
    parmed_structure.save(outname + '.gro', overwrite=True)

def to_amber_system(system, outname):
    parmed_structure = pmd.openmm.load_topology(
        system.modeller.topology, 
        system=system.solvated_system, 
        xyz=system.modeller.positions
    )
    # adding default bond types where bond type is none
    bond_type = pmd.BondType(1.0, 1.0, list=parmed_structure.bond_types)
    parmed_structure.bond_types.append(bond_type)
    for bond in parmed_structure.bonds:
        if bond.type is None:
            bond.type = bond_type

    #parmed_structure.save(outname + '.parm7', format='amber')
    #parmed_structure.save(outname + '.rst7', format='rst7')
    parmed_structure.save(outname + '.parm7', overwrite=True)
    parmed_structure.save(outname + '.rst7', overwrite=True)

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Description of your program')
    #parser.add_argument('-f','--foo', help='Description for foo argument', required=True)
    #parser.add_argument('-b','--bar', help='Description for bar argument', required=True)
    parser.add_argument("PROTEIN", help=" Protein .pdb")
    parser.add_argument("LIGAND", help=" Target Ligand .sdf")
    parser.add_argument("LIGMODEL", help=" Small Molecule FF Generator to be used as: GAFF, SMIRNOFF, ESPALOMA")
    parser.add_argument("OUTNAME", help=" Output name used for GROMACS or AMBER system")
    parser.add_argument('-amber', action='store_true')
    args = parser.parse_args()

    prot_name = args.PROTEIN.split(".pdb")[0]
    fixer = fix_protein(prot_name + '.pdb')
    PDBFile.writeFile(fixer.topology, fixer.positions, open(prot_name +'_fixed.pdb', 'w'))
    lig_name = args.LIGAND.split(".sdf")[0]

    config=MDConfig("")
    config.constraints = None
    
    if args.LIGMODEL == "GAFF":
        config.ligand_forcefield_generator = gaff_generator
    elif args.LIGMODEL == "SMIRNOFF":
        config.ligand_forcefield_generator = smirnoff_generator
    elif args.LIGMODEL == "ESPALOMA":
        config.ligand_forcefield_generator = espaloma_generator
    else:
        print("LIGMODEL has to be 'GAFF', 'SMIRNOFF' or 'ESPALOMA'")
        exit()

    prepared_system = createProteinLigandOnlySystem(
        prot_name +'_fixed.pdb',
        lig_name + ".sdf",
        config=config
    )
    print(f'Output parametrized system for {prot_name} and {lig_name}') 
    print(f' use created PDB {prot_name}_fixed.pdb')
    if not args.amber:
        to_gromacs_system(prepared_system, args.OUTNAME)
    else:
        to_amber_system(prepared_system, args.OUTNAME)



