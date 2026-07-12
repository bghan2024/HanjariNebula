# -*- coding: utf-8 -*-
"""
HanjariNebulaEngine_OS (Open-Source / OpenMM Edition)
---------------------------------------------
A high-precision hybrid physical refinement engine designed for structural biology pipelines.
Part of the HanjariNebula Engine OS Suite.

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
import string
import numpy as np

try:
    from openmm.app import PDBFile, PDBxFile, ForceField, CutoffNonPeriodic, Modeller, Simulation, Topology
    from openmm import LocalEnergyMinimizer, Platform, LangevinIntegrator, Context, Vec3
    import openmm as mm
    import openmm.unit as unit
except ImportError:
    print("❌ OpenMM library is missing. Install it with: conda install -c conda-forge openmm")
    sys.exit(1)

class HanjariNebulaEngine_OS:
    """
    HanjariNebulaEngine_OS
    ----------------------
    A highly precise structural refinement engine designed for structural biology pipelines.

    Key Features:
    1. Perfect Index Mapping: Uses a unique chain-residue key mapping to track and lock down heavy atoms,
       preventing coordinate scrambling and index shifts when OpenMM adds hydrogens.
    2. Stereanalyst Patch: Restrains protein heavy atoms using custom harmonic forces to prevent
       bond-tearing, C-beta deviations, or non-physical distortions during minimization.
    3. Phenix-Ready Hydrogen Adjustment: Automatically scales polar and non-polar hydrogen bonds
       to standardized chemistry lengths for downstream validation tools.
    4. mmCIF Quality Patches: Automatically reformats mmCIF output formatting (occupancies and B-factors)
       to ensure compatibility with structural validation suites.
    """
    def __init__(self, max_iterations=300):
        self.max_iterations = max_iterations
        self.forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')

    def adjust_hydrogens_for_phenix(self, topology, positions):
        """
        Dynamically adjusts hydrogen bond lengths to match the target chemistry standards
        expected by Phenix and MolProbity.
        """
        pos_nm = list(positions.value_in_unit(unit.nanometers))
        h_parent_map = {}
        for bond in topology.bonds():
            a1, a2 = bond.atom1, bond.atom2
            if a1.element is not None and a2.element is not None:
                if a1.element.symbol == 'H' and a2.element.symbol != 'H': 
                    h_parent_map[a1.index] = a2
                elif a2.element.symbol == 'H' and a1.element.symbol != 'H': 
                    h_parent_map[a2.index] = a1
                
        for h_idx, parent_atom in h_parent_map.items():
            p_heavy = np.array(pos_nm[parent_atom.index])
            p_h = np.array(pos_nm[h_idx])
            bond_vector = p_h - p_heavy
            current_dist = np.linalg.norm(bond_vector)
            if current_dist == 0: 
                continue
            
            heavy_symbol = parent_atom.element.symbol
            if heavy_symbol == 'N': 
                target_dist = 0.0860
            elif heavy_symbol == 'C': 
                target_dist = 0.0930
            elif heavy_symbol == 'O': 
                target_dist = 0.0840
            elif heavy_symbol == 'S': 
                target_dist = 0.1000
            else: 
                target_dist = max(0.01, current_dist - 0.0156)
            
            refined_h_pos = p_heavy + bond_vector * (target_dist / current_dist)
            pos_nm[h_idx] = Vec3(refined_h_pos[0], refined_h_pos[1], refined_h_pos[2])
            
        return pos_nm * unit.nanometers

    def minimize_structure(self, input_path, output_path):
        print(f"[NebulaEngine_OS] Starting physical energy minimization... (Input: {os.path.basename(input_path)})")
        
        ext = os.path.splitext(input_path)[-1].lower()
        structure_obj = PDBxFile(input_path) if ext in ['.cif', '.mmcif'] else PDBFile(input_path)
        old_topology = structure_obj.topology
        old_pos_nm = structure_obj.positions.value_in_unit(unit.nanometers)
        
        # 💡 [CRITICAL BUG FIX]: Maps original heavy atoms to coordinates using chain and residue sequence index
        # to prevent key collisions from duplicate residue IDs (e.g. from multiple hydration molecules or heterogens).
        old_coords_map = {}
        for chain_idx, old_chain in enumerate(old_topology.chains()):
            for res_idx, old_res in enumerate(old_chain.residues()):
                for old_atom in old_res.atoms():
                    if old_atom.element is not None and old_atom.element.symbol != 'H' and old_res.name != 'HOH':
                        key = (chain_idx, res_idx, old_atom.name)
                        old_coords_map[key] = old_pos_nm[old_atom.index]

        print("  🧹 [NebulaGuard] Sanitizing topology and rebuilding representation...")
        new_topology = Topology()
        new_positions = []
        atom_mapping = {}
        
        used_chain_ids = set()
        available_letters = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [f"C{i}" for i in range(1, 100)]
        letter_idx = 0
        
        for old_chain in old_topology.chains():
            proposed_id = old_chain.id.strip() if (old_chain.id and old_chain.id.strip()) else "A"
            if proposed_id in used_chain_ids:
                while letter_idx < len(available_letters):
                    candidate = available_letters[letter_idx]
                    letter_idx += 1
                    if candidate not in used_chain_ids:
                        proposed_id = candidate
                        break
            
            used_chain_ids.add(proposed_id)
            new_chain = new_topology.addChain(proposed_id)
            
            # Map each residue uniquely by avoiding duplicate-ID merging bugs
            for old_res in old_chain.residues():
                new_res = new_topology.addResidue(old_res.name, new_chain, old_res.id)
                seen_atom_names = set()
                for old_atom in old_res.atoms():
                    if old_atom.element is not None and old_atom.element.symbol == 'H' and old_res.name != 'HOH':
                        continue
                    if old_atom.name not in seen_atom_names:
                        new_atom = new_topology.addAtom(old_atom.name, old_atom.element, new_res)
                        atom_mapping[old_atom] = new_atom
                        new_positions.append(old_pos_nm[old_atom.index])
                        seen_atom_names.add(old_atom.name)
                    else:
                        orig_atom = [a for a in new_res.atoms() if a.name == old_atom.name][0]
                        atom_mapping[old_atom] = orig_atom

        existing_new_bonds = set()
        for old_bond in old_topology.bonds():
            a1 = atom_mapping.get(old_bond.atom1)
            a2 = atom_mapping.get(old_bond.atom2)
            if a1 is not None and a2 is not None and a1 != a2:
                id1, id2 = min(a1.index, a2.index), max(a1.index, a2.index)
                if (id1, id2) not in existing_new_bonds:
                    new_topology.addBond(a1, a2)
                    existing_new_bonds.add((id1, id2))
                    
        modeller = Modeller(new_topology, new_positions * unit.nanometers)
        modeller.addHydrogens(forcefield=self.forcefield)
        
        system = self.forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.2, constraints=None)
        
        # Heavy-atom restraint force (Anti-Drift Force)
        restraint_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint_force.addGlobalParameter("k", 50.0 * unit.kilojoules_per_mole/unit.nanometer**2)
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")
        
        # Restrain all non-hydrogen, non-water heavy atoms to their original positions
        for final_chain_idx, final_chain in enumerate(modeller.topology.chains()):
            for final_res_idx, final_res in enumerate(final_chain.residues()):
                for final_atom in final_res.atoms():
                    if final_atom.element is not None and final_atom.element.symbol != 'H' and final_res.name != 'HOH':
                        match_key = (final_chain_idx, final_res_idx, final_atom.name)
                        if match_key in old_coords_map:
                            restraint_force.addParticle(final_atom.index, old_coords_map[match_key])

        if restraint_force.getNumParticles() > 0: 
            system.addForce(restraint_force)

        integrator = LangevinIntegrator(0, 0, 0)
        sim = Simulation(modeller.topology, system, integrator, Platform.getPlatformByName('CPU'))
        sim.context.setPositions(modeller.positions)
        
        print(f"  ⚡ Running L-BFGS energy minimization (Max iterations: {self.max_iterations})...")
        LocalEnergyMinimizer.minimize(sim.context, 20.0, self.max_iterations)
        
        state = sim.context.getState(getPositions=True)
        phenix_ready_positions = self.adjust_hydrogens_for_phenix(modeller.topology, state.getPositions())
        
        with open(output_path, 'w') as f:
            PDBxFile.writeFile(modeller.topology, phenix_ready_positions, f)
            
        print("  🎯 [NebulaGuard] Aligning mmCIF metadata layout...")
        with open(output_path, 'r') as f:
            lines = f.readlines()
            
        columns = {}
        header_idx = 0
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('_atom_site.'):
                columns[stripped] = header_idx
                header_idx += 1
            elif stripped.startswith('ATOM') or stripped.startswith('HETATM'):
                break
            elif stripped == 'loop_':
                header_idx = 0
                columns = {}

        occ_idx = columns.get('_atom_site.occupancy', 13)
        b_idx = columns.get('_atom_site.B_iso_or_equiv', 14)
        
        patched_lines = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                tokens = line.split()
                if occ_idx < len(tokens): 
                    tokens[occ_idx] = '1.00'
                if b_idx < len(tokens): 
                    tokens[b_idx] = '20.00'
                line = " ".join(tokens) + "\n"
            patched_lines.append(line)
            
        with open(output_path, 'w') as f:
            f.writelines(patched_lines)
                
        print(f"[NebulaEngine_OS] Physical refinement completed. Output file: {output_path}")
        return output_path
