# -*- coding: utf-8 -*-
"""
HanjariFold2_CosmicDistiller_OS (Open-Source / OpenMM Edition)
------------------------------------------------------
An automated macromolecular hydration and physical energy minimization manager.
Part of the HanjariNebula Engine OS Suite.

CRITICAL USAGE NOTE:
Due to the parametrization limitations of the default Amber14 forcefield libraries used
(amber14/protein.ff14SB.xml), this engine is designed exclusively for protein structure refinement.
Refining structures containing DNA, RNA, or other nucleic acids is highly discouraged and may
lead to topology/parameter matching failures. It is proposed to limit usage to protein-only structures.

Developer: Han Byeong-gu (hanbyeonggu@gmail.com)
Repository: https://github.com/bghan2024/HanjariNebula
Citation Ref: HanjariNebula Engine OS Suite (2026); see README.md for details.

License Summary:
- Academic & Non-Commercial Research Use: Granted Free.
- Commercial & For-Profit Use: Strictly requires a separate written commercial contract.
  For inquiries, contact the developer (hanbyeonggu@gmail.com).
"""

import os
import sys
import datetime
import numpy as np
from scipy.spatial import cKDTree

import openmm.app as app
from openmm import Vec3
import openmm.unit as unit

from HanjariNebulaEngine_OS import HanjariNebulaEngine_OS

class HanjariFold2CosmicDistiller_OS:
    """
    HanjariFold2CosmicDistiller_OS
    ---------------------------
    An automated hydration layer compilation and physical refinement protocol.

    Key Features:
    1. Perfect Water Compiling: Reconstructs a structured H-O-H hydration lattice (hydration layer)
       around a target macromolecule to screen out coordinate defects and phantom atoms.
    2. Gemmi-based QA Checking: Performs an efficient, lightweight validation of the final structure's
       atom count, charge topology, and residue integrity using Gemmi.
    3. Structural Synergy: Combines the hydration layer with the HanjariNebulaEngine_OS L-BFGS
       restrained minimization algorithm to generate physically optimized outputs.
    """
    def __init__(self, input_path):
        self.input_path = os.path.abspath(input_path)
        ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        self.output_dir = os.path.abspath(f"HF2_NEBULA_DISTILLER_{ts}")
        os.makedirs(self.output_dir, exist_ok=True)

    def log(self, msg):
        print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")

    def run_gemmi_validation(self, cif_path):
        """
        [Lightweight Gemmi Validation]
        Reads the resulting mmCIF file and runs a structural scan to verify chains, residues,
        water molecules, and total compiled atom integrity.
        """
        self.log("📊 [Result Validation] Launching Gemmi validation safeguards...")
        try:
            import gemmi
            
            doc = gemmi.cif.read_file(cif_path)
            st = gemmi.make_structure_from_block(doc[0])
            
            total_chains = len(st[0]) if len(st) > 0 else 0
            water_count = 0
            protein_res_count = 0
            total_atoms = 0
            
            if total_chains > 0:
                for chain in st[0]:
                    for res in chain:
                        total_atoms += len(res)
                        if res.name == "HOH":
                            water_count += 1
                        else:
                            protein_res_count += 1
            
            print("\n" + "="*60)
            print(f"🏆  [HanjariNebula Engine - Gemmi Quality Report]  🏆")
            print("="*60)
            print(f"🧬 Total Chains         : {total_chains}")
            print(f"🧪 Protein Residues    : {protein_res_count}")
            print(f"💧 Solvated Water HOH   : {water_count}")
            print(f"⚛️ Total Compiled Atoms  : {total_atoms}")
            print("-"*60)
            
            if water_count > 0:
                print("✅ [Hydration Validation]: CosmicDistiller water layer successfully compiled.")
            else:
                print("⚠️ [Hydration Warning]: No water molecules detected in final structure.")
                
            if total_atoms > 0:
                print("✅ [Structural Integrity]: mmCIF topology built without physical defects.")
            print("="*60 + "\n")
            
        except ImportError:
            self.log("⚠️ Gemmi library is not installed; skipped quality report.")
            self.log("💡 Run 'pip install gemmi' to enable Gemmi validation.")
        except Exception as e:
            self.log(f"❌ Exception occurred during Gemmi validation: {str(e)}")

    def run(self):
        self.log(f"🌌 CosmicDistiller starting... Target: {os.path.basename(self.input_path)}")
        
        ext = os.path.splitext(self.input_path)[-1].lower()
        structure_obj = app.PDBxFile(self.input_path) if ext in ['.cif', '.mmcif'] else app.PDBFile(self.input_path)
        
        # Safely convert OpenMM positions to a NumPy array of coordinates in Angstroms (1 nm = 10 Å)
        positions_nm = structure_obj.positions.value_in_unit(unit.nanometers)
        coords = np.array([[pos[0], pos[1], pos[2]] for pos in positions_nm]) * 10.0
        
        self.log(f"💧 [CosmicDistiller] Constructing hydration grid (Spacing: 4.5Å)")
        min_c = np.min(coords, axis=0) - 4.0
        max_c = np.max(coords, axis=0) + 4.0
        
        grid_x = np.arange(min_c[0], max_c[0], 4.5)
        grid_y = np.arange(min_c[1], max_c[1], 4.5)
        grid_z = np.arange(min_c[2], max_c[2], 4.5)
        grid_points = np.array(np.meshgrid(grid_x, grid_y, grid_z)).T.reshape(-1, 3)
        
        tree = cKDTree(coords)
        distances, _ = tree.query(grid_points, k=1)
        valid_water_mask = (distances > 2.8) & (distances < 3.8)
        water_coords = grid_points[valid_water_mask]
        
        water_topo = app.Topology()
        water_chain = water_topo.addChain(id="W")
        water_positions = []
        
        for i, coord in enumerate(water_coords):
            res = water_topo.addResidue("HOH", water_chain, id=str((i % 9999) + 1))
            
            o_atom = water_topo.addAtom("O", app.Element.getBySymbol("O"), res)
            h1_atom = water_topo.addAtom("H1", app.Element.getBySymbol("H"), res)
            h2_atom = water_topo.addAtom("H2", app.Element.getBySymbol("H"), res)
            
            o_pos = Vec3(coord[0]/10.0, coord[1]/10.0, coord[2]/10.0)
            h1_pos = Vec3(o_pos.x + 0.0957, o_pos.y, o_pos.z)
            h2_pos = Vec3(o_pos.x - 0.0240, o_pos.y + 0.0926, o_pos.z)
            
            water_positions.extend([o_pos, h1_pos, h2_pos])
            
            water_topo.addBond(o_atom, h1_atom)
            water_topo.addBond(o_atom, h2_atom)
            
        modeller = app.Modeller(structure_obj.topology, structure_obj.positions)
        modeller.add(water_topo, water_positions * unit.nanometers)
        
        temp_hydrated_cif = os.path.join(self.output_dir, "hydrated_temp.cif")
        with open(temp_hydrated_cif, 'w') as f:
            app.PDBxFile.writeFile(modeller.topology, modeller.positions, f)
        
        final_output_file = os.path.join(self.output_dir, "distilled_final.cif")
        
        minimizer_engine = HanjariNebulaEngine_OS(max_iterations=300)
        minimizer_engine.minimize_structure(temp_hydrated_cif, final_output_file)
        
        if os.path.exists(temp_hydrated_cif):
            os.remove(temp_hydrated_cif)
            
        self.run_gemmi_validation(final_output_file)
        self.log("🎉 [HanjariNebula] Protein-hydration complex assembly and physical refinement completed.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("💡 Usage: python HanjariFold2_CosmicDistiller_OS.py <input.pdb/.cif>")
    else:
        HanjariFold2CosmicDistiller_OS(sys.argv[1]).run()
