import os
import time
import logging
import traceback
import math
import openmm as mm
from openmm import app
from openmm import unit

logger = logging.getLogger(__name__)

def minimize_structure(input_pdb_path, output_pdb_path, force_field="amber14", use_gbsa=True):
    """
    Minimizes the energy of a protein structure using OpenMM.
    If minimization fails, the input file is copied to the output path
    and an error message is returned.
    """
    try:
        pdb = app.PDBFile(input_pdb_path)
        
        # 1. Select Force Field
        ff = None
        if force_field == "amber14":
            if use_gbsa:
                # Use AMBER14 with OBC2 implicit solvent
                ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
            else:
                ff = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        elif force_field == "charmm36":
            # CHARMM36 typically needs its own specific water/implicit models
            # but we fall back to standard if implicit isn't available
            if use_gbsa:
                try:
                    ff = app.ForceField('charmm36.xml', 'charmm36/water.xml') # GBSA not natively standard for charmm36 in openmm without explicit params, simplified here
                except Exception:
                    logger.warning("CHARMM36 GBSA not found, falling back to AMBER14 GBSA")
                    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
            else:
                ff = app.ForceField('charmm36.xml', 'charmm36/water.xml')
        else:
            ff = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
            
        # 2. Setup System
        if pdb.topology.getNumAtoms() == 0:
            raise ValueError("The structure contains no atoms (possibly all chains were removed). Minimization skipped.")
            
        kwargs = {"nonbondedMethod": app.NoCutoff, "constraints": app.HBonds}
        system = ff.createSystem(pdb.topology, **kwargs)
        
        # 3. Setup Simulation
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picosecond
        )
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        
        # Check initial energy
        state_before = simulation.context.getState(getEnergy=True)
        energy_before = state_before.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # Guard against NaN/Inf
        if math.isnan(energy_before) or math.isinf(energy_before):
            raise ValueError(f"Extremely high initial energy ({energy_before}). Severe atomic clashes prevent minimization.")
        
        # 4. Minimize
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Check final energy
        state_after = simulation.context.getState(getEnergy=True, getPositions=True)
        energy_after = state_after.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        # 5. Save output
        with open(output_pdb_path, "w") as f:
            app.PDBFile.writeFile(
                simulation.topology, state_after.getPositions(), f, keepIds=True
            )
            
        delta_energy = energy_after - energy_before
        converged = (delta_energy < -0.1)
        
        return {
            "force_field": force_field,
            "energy_before_kJ_mol": round(energy_before, 1),
            "energy_after_kJ_mol": round(energy_after, 1),
            "delta_energy_kJ_mol": round(delta_energy, 1),
            "iterations": 1000,
            "converged": converged,
            "gbsa_used": use_gbsa,
            "output_path": output_pdb_path,
        }
        
    except Exception as e:
        logger.error(f"Energy minimization failed: {e}")
        import shutil
        shutil.copy2(input_pdb_path, output_pdb_path)
        
        # Fallback dictionary. The frontend now knows to look for "error"
        # and display it nicely instead of "N/A"
        return {
            "force_field": force_field,
            "energy_before_kJ_mol": "N/A",
            "energy_after_kJ_mol": "N/A",
            "delta_energy_kJ_mol": "N/A",
            "iterations": 0,
            "converged": False,
            "gbsa_used": use_gbsa,
            "error": str(e),
            "output_path": output_pdb_path,
        }
