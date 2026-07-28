# -*- coding: utf-8 -*-
"""
HanjariCosmicStandalone (v40.2 Cosmic Precision Release)

================================================================================
LICENSE & TERMS OF USE:
- ACADEMIC & NON-COMMERCIAL USE ONLY: This software is free for academic,
  educational, and non-profit research purposes.
- Commercial utilization is strictly prohibited without a separate, prior
  written agreement with the author (Byeong-gu Han / bghan2024).
- The "OS" suffix or Native Engine denotes OpenMM Solvation Support and 
  does NOT imply Open Source Initiative (OSI) compliance.
================================================================================

Key Features:
1. Phenix / CCTBX Independent: Runs 100% natively on OpenMM backends.
2. Precision Water Mesh: Constructs 104.52° bent explicit hydration layers in memory.
3. Robust Hybrid Engine: Bypasses OpenMM template matching errors.
"""

import os
import sys
import datetime
import numpy as np
import shutil
import argparse
from scipy.spatial import cKDTree

try:
    from openmm import app, unit
    from openmm.app import PDBFile, PDBxFile, Modeller, Element, Topology
    from openmm import Vec3
except ImportError:
    print("❌ Missing OpenMM library. Please install it using 'conda install -c conda-forge openmm'.")
    sys.exit(1)

# Import the core potential builder from the same directory
from HanjariAetherEngine import HanjariAetherEngine

class HanjariCosmicStandalone:
    """
    Main coordinator class for hydration and refinement using HanjariCosmicEngine.
    Loads input structure, generates explicit hydration shell with ideal geometries,
    and runs the L-BFGS refinement core.
    """
    def __init__(self, input_file, target_ph=7.4):
        self.input_path = os.path.abspath(input_file)
        self.ext = os.path.splitext(self.input_path)[1].lower()
        
        self.base_dir = os.path.dirname(self.input_path)
        self.target_ph = target_ph
        
        self.timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        self.output_dir = os.path.join(self.base_dir, f"V40_COSMIC_{self.timestamp}")
        self.temp_dir = os.path.join(self.output_dir, "temp_workspace")
        
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.temp_dir, exist_ok=True)

    def log(self, message):
        """Log a timestamped message to standard output."""
        ts = datetime.datetime.now().strftime('%H:%M:%S')
        print(f"[{ts}] {message}")

    def step1_and_2_native_load_and_hydration(self):
        """
        Step 1: Load the 3D structure using the native OpenMM parser (PDB or PDBx/mmCIF).
        Step 2: Construct a precision explicit water hydration layer using a k-d tree spatial query.
        """
        self.log(f"Step 1: Loading structure with native OpenMM parser ({self.ext.upper()})")
        
        if self.ext in ['.cif', '.mmcif']:
            structure_obj = PDBxFile(self.input_path)
        else:
            structure_obj = PDBFile(self.input_path)
            
        self.log("Step 2: Cosmic Precision Hydration (Weaving 104.52° bent explicit water mesh)")
        
        # Extract atomic coordinates in Angstroms
        coords_angstrom = np.array([[pos.x, pos.y, pos.z] for pos in structure_obj.positions.value_in_unit(unit.angstroms)])
        
        prot_coords, polar_atoms = [], []
        for atom, pos in zip(structure_obj.topology.atoms(), coords_angstrom):
            if atom.residue.name != 'HOH':
                prot_coords.append(pos)
                elem_symbol = atom.element.symbol.upper() if atom.element else ''
                if elem_symbol in ['O', 'N']:
                    polar_atoms.append(pos)
                    
        prot_coords = np.array(prot_coords, dtype=np.float32)
        polar_atoms = np.array(polar_atoms, dtype=np.float32)
        
        # Build spatial k-d tree for fast proximity search
        tree = cKDTree(prot_coords) if len(prot_coords) > 0 else None
        added_waters = []
        
        # Generate search directions (26 unit directions in 3D grid)
        directions = [[x, y, z] for x in [-1, 0, 1] for y in [-1, 0, 1] for z in [-1, 0, 1] if not (x == 0 and y == 0 and z == 0)]
        directions = (np.array(directions) / np.linalg.norm(directions, axis=1)[:, None]) * 3.05  # Radial distance ~3.05 Å
        
        if tree is not None:
            for p_coord in polar_atoms:
                for d in directions:
                    cand = p_coord + d
                    dist, _ = tree.query(cand)
                    # Keep candidate if distance to protein is between 2.78 and 3.30 Å
                    if 2.78 <= dist <= 3.30:
                        # Prevent clustering water too close to existing water molecules
                        if len(added_waters) == 0 or np.min(np.linalg.norm(np.array(added_waters) - cand, axis=1)) > 2.70:
                            added_waters.append(cand)
                            
        self.log(f"  -> Successfully wove {len(added_waters)} explicit water molecules in-memory")
        
        water_topo = Topology()
        water_chain = water_topo.addChain(id="W")
        water_positions = []
        
        # Set ideal water geometry: r_OH = 0.09572 nm, H-O-H angle = 104.52°
        r_oh = 0.09572  # nanometers
        half_angle_rad = np.radians(104.52 / 2.0)  # 52.26 degrees
        dx = r_oh * np.sin(half_angle_rad)          # ~0.0757 nm
        dz = r_oh * np.cos(half_angle_rad)          # ~0.0586 nm

        # Instantiate atoms and bonds for each generated water molecule
        for i, coord in enumerate(added_waters):
            res = water_topo.addResidue("HOH", water_chain, id=str((i % 9999) + 1))
            o_atom = water_topo.addAtom("O", Element.getBySymbol("O"), res)
            h1_atom = water_topo.addAtom("H1", Element.getBySymbol("H"), res)
            h2_atom = water_topo.addAtom("H2", Element.getBySymbol("H"), res)
            
            o_pos = Vec3(coord[0]/10.0, coord[1]/10.0, coord[2]/10.0) # convert from Å to nm
            
            # Place hydrogens symmetrical to oxygen along Z and X axes
            h1_pos = Vec3(o_pos.x + dx, o_pos.y, o_pos.z + dz)
            h2_pos = Vec3(o_pos.x - dx, o_pos.y, o_pos.z + dz)
            
            water_positions.extend([o_pos, h1_pos, h2_pos])
            water_topo.addBond(o_atom, h1_atom)
            water_topo.addBond(o_atom, h2_atom)
            
        # Merge protein and water topologies
        modeller = Modeller(structure_obj.topology, structure_obj.positions)
        modeller.add(water_topo, water_positions * unit.nanometers)
        
        self.hydrated_topology = modeller.topology
        self.hydrated_positions = modeller.positions

    def step3_refinement(self):
        """
        Step 3: Run OpenMM L-BFGS minimization using the HanjariAetherEngine force constraints.
        """
        self.log("Step 3: Activating HanjariAetherEngine Precision Refinement")
        final_output = os.path.join(self.output_dir, "cosmic_refined_final.cif")
        
        engine = HanjariAetherEngine(max_iterations=350)
        engine.minimize_topology_and_positions(self.hydrated_topology, self.hydrated_positions, final_output)
        
        shutil.rmtree(self.temp_dir)
        print("\n" + "="*80)
        print(f" 🌌 HANJARI COSMIC ENGINE REFINEMENT COMPLETE")
        print(f" Output: {final_output}")
        print("="*80 + "\n")

    def execute(self):
        """Execute the entire pipeline."""
        self.log("🚀 HANJARI COSMIC v40.2 - NATIVE ENGINE START")
        self.step1_and_2_native_load_and_hydration()
        self.step3_refinement()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Hanjari Cosmic Precision Refinement Pipeline")
    parser.add_argument("input_file", help="Input PDB or CIF file")
    parser.add_argument("--ph", type=float, default=7.4)
    args = parser.parse_args()
    
    HanjariCosmicStandalone(args.input_file, target_ph=args.ph).execute()
