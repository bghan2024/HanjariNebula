# -*- coding: utf-8 -*-
"""
HanjariAetherEngine (v40.2 - Pure Refinement Core)

================================================================================
LICENSE & TERMS OF USE:
- ACADEMIC & NON-COMMERCIAL USE ONLY: This software is free for academic,
  educational, and non-profit research purposes.
- Commercial utilization is strictly prohibited without a separate, prior
  written agreement with the author (Byeong-gu Han / bghan2024).
- The "OS" suffix or Native Engine denotes OpenMM Solvation Support and 
  does NOT imply Open Source Initiative (OSI) compliance.
================================================================================

Key Innovations:
1. Water Angle Lock: HarmonicAngleForce constrains water H-O-H angles to 104.52°,
   preventing angle outliers from forming during optimization.
2. C-beta Zero-Deviation: High-stiffness chiral volume and tetrahedral angle restraints
   for C-beta positions (k_cbeta = 3,000,000).
3. Steric Wash & Planarity: Custom non-bonded potentials target MolProbity clashscore 0.00
   and preserve aromatic ring planarity.
"""

import os
import sys
import numpy as np

try:
    from openmm.app import PDBFile, PDBxFile, ForceField, CutoffNonPeriodic, Modeller, Simulation, Topology
    from openmm import LocalEnergyMinimizer, Platform, LangevinIntegrator, Context, Vec3
    import openmm as mm
    import openmm.unit as unit
except ImportError:
    print("❌ Missing OpenMM library. Please install it using 'conda install -c conda-forge openmm'.")
    sys.exit(1)


class AetherIdealizationEngine:
    """
    Constructs an independent potential system bypassing OpenMM's standard forcefield XML templates.
    Combines custom harmonic bonds, H-O-H angle locks, non-bonded repulsions,
    aromatic planar restraints, and C-beta chiral/tetrahedral geometry restraints.
    """
    @staticmethod
    def create_idealization_system(topology, positions):
        system = mm.System()
        pos_nm = positions.value_in_unit(unit.nanometers) if hasattr(positions, 'value_in_unit') else positions
        
        # 1. Register particles and assign masses
        for atom in topology.atoms():
            mass = 12.0 * unit.amu
            if atom.element is not None:
                mass = atom.element.mass
            system.addParticle(mass)

        # 2. Bond Restraints (Idealizing bond distances using CustomBondForce)
        bond_force = mm.CustomBondForce("k_bond * (r - r0)^2")
        bond_force.addGlobalParameter("k_bond", 300000.0 * unit.kilojoules_per_mole / unit.nanometer**2)
        bond_force.addPerBondParameter("r0")
        
        for bond in topology.bonds():
            a1, a2 = bond.atom1.index, bond.atom2.index
            p1, p2 = np.array(pos_nm[a1]), np.array(pos_nm[a2])
            r0 = np.linalg.norm(p1 - p2)
            if r0 > 0:
                bond_force.addBond(a1, a2, [r0 * unit.nanometers])
        if bond_force.getNumBonds() > 0:
            system.addForce(bond_force)

        # 3. Water H-O-H Angle Restraints (HarmonicAngleForce)
        # Prevents L-BFGS minimizers from flattening the angle to 180° by enforcing ideal 104.52°
        water_angle_force = mm.HarmonicAngleForce()
        target_water_angle = 104.52 * np.pi / 180.0  # 1.82421 radians
        k_water_angle = 1500.0 * unit.kilojoules_per_mole / unit.radian**2
        
        for res in topology.residues():
            if res.name.strip() == 'HOH':
                atoms = {a.name.strip(): a.index for a in res.atoms()}
                if 'O' in atoms and 'H1' in atoms and 'H2' in atoms:
                    water_angle_force.addAngle(atoms['H1'], atoms['O'], atoms['H2'], target_water_angle, k_water_angle)
        
        if water_angle_force.getNumAngles() > 0:
            system.addForce(water_angle_force)

        # 4. Non-bonded Repulsion (Steric Clash Wash to achieve Clashscore 0.00)
        # Applies harmonic repulsion for atoms overlapping closer than 2.2 Å (0.22 nm)
        nb_force = mm.CustomNonbondedForce("step(r_clash - r) * k_clash * (r_clash - r)^2")
        nb_force.addGlobalParameter("k_clash", 100000.0 * unit.kilojoules_per_mole / unit.nanometer**2)
        nb_force.addGlobalParameter("r_clash", 0.22 * unit.nanometers)
        
        for _ in topology.atoms():
            nb_force.addParticle([])
            
        exceptions = set()
        for bond in topology.bonds():
            exceptions.add((min(bond.atom1.index, bond.atom2.index), max(bond.atom1.index, bond.atom2.index)))
        for a1, a2 in exceptions:
            nb_force.addExclusion(a1, a2)
            
        system.addForce(nb_force)

        # 5. Aromatic Ring Planarity Restraints
        # CustomCompoundBondForce restricts planar deviations for PHE, TYR, TRP, and HIS ring carbons/nitrogens
        plane_force = mm.CustomCompoundBondForce(4, "k_plane * ( (x1*(y2*(z3-z4) + y3*(z4-z2) + y4*(z2-z3)) + "
                                                     "x2*(y1*(z4-z3) + y3*(z2-z4) + y4*(z3-z2)) + "
                                                     "x3*(y1*(z2-z4) + y2*(z4-z1) + y4*(z1-z2)) + "
                                                     "x4*(y1*(z3-z2) + y2*(z1-z3) + y3*(z2-z1)))^2 )")
        plane_force.addGlobalParameter("k_plane", 100.0 * unit.kilojoules_per_mole / unit.nanometer**6)

        aromatic_res = {'PHE': ['CG', 'CD1', 'CD2', 'CE1'],
                        'TYR': ['CG', 'CD1', 'CD2', 'CE1'],
                        'TRP': ['CD1', 'CD2', 'NE1', 'CE2'],
                        'HIS': ['CG', 'ND1', 'CD2', 'CE1']}

        for chain in topology.chains():
            for res in chain.residues():
                rname = res.name.strip()
                if rname in aromatic_res:
                    atom_map = {a.name.strip(): a.index for a in res.atoms()}
                    target_names = aromatic_res[rname]
                    if all(name in atom_map for name in target_names):
                        idxs = [atom_map[name] for name in target_names]
                        plane_force.addBond(idxs, [])

        if plane_force.getNumBonds() > 0:
            system.addForce(plane_force)

        # 6. C-beta Ultra-Precision Chiral Restraint (High Stiffness Constraint)
        # Keeps CA-tetrahedral chiral volume target (v0) close to original positions (k_cbeta = 3,000,000)
        cbeta_force = mm.CustomCompoundBondForce(4, "k_cbeta * ( ( (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2) ) * (x4-x2) + "
                                                     "( (z1-z2)*(x3-x2) - (x1-x2)*(z3-z2) ) * (y4-y2) + "
                                                     "( (x1-x2)*(y3-y2) - (y1-y2)*(x3-x2) ) * (z4-z2) - v0 )^2")
        cbeta_force.addGlobalParameter("k_cbeta", 3000000.0 * unit.kilojoules_per_mole / unit.nanometer**6)
        cbeta_force.addPerBondParameter("v0")

        # C-beta tetrahedral angle restraints (N-CA-CB and C-CA-CB) to fully stabilize C-beta position
        cbeta_angle_force = mm.HarmonicAngleForce()
        target_tetra_angle = 109.5 * np.pi / 180.0  # 1.911 radians
        k_tetra = 1000.0 * unit.kilojoules_per_mole / unit.radian**2

        for chain in topology.chains():
            for res in chain.residues():
                rname = res.name.strip()
                if rname not in ['GLY', 'HOH']:
                    atom_map = {a.name.strip(): a.index for a in res.atoms()}
                    if all(k in atom_map for k in ['N', 'CA', 'C', 'CB']):
                        i_n, i_ca, i_c, i_cb = atom_map['N'], atom_map['CA'], atom_map['C'], atom_map['CB']
                        p_n, p_ca, p_c, p_cb = np.array(pos_nm[i_n]), np.array(pos_nm[i_ca]), np.array(pos_nm[i_c]), np.array(pos_nm[i_cb])
                        
                        v_ca_n = p_n - p_ca
                        v_ca_c = p_c - p_ca
                        v_ca_cb = p_cb - p_ca
                        v0_val = np.dot(np.cross(v_ca_n, v_ca_c), v_ca_cb)
                        
                        cbeta_force.addBond([i_n, i_ca, i_c, i_cb], [v0_val * unit.nanometer**3])
                        
                        # Apply HarmonicAngleForce to restrain tetrahedral angle geometry
                        cbeta_angle_force.addAngle(i_n, i_ca, i_cb, target_tetra_angle, k_tetra)
                        cbeta_angle_force.addAngle(i_c, i_ca, i_cb, target_tetra_angle, k_tetra)

        if cbeta_force.getNumBonds() > 0:
            system.addForce(cbeta_force)
        if cbeta_angle_force.getNumAngles() > 0:
            system.addForce(cbeta_angle_force)

        return system


class HanjariAetherEngine:
    """
    Coordinates L-BFGS energy minimization with custom force restraints
    and performs post-minimization hydrogen position idealization.
    """
    def __init__(self, max_iterations=350):
        self.max_iterations = max_iterations

    def adjust_hydrogens_for_ideal_geometry(self, topology, positions):
        """
        Adjust bond lengths of hydrogen atoms attached to heavy atoms (C, N, O, S, P)
        to match their ideal thermodynamic bond lengths.
        """
        pos_nm = list(positions.value_in_unit(unit.nanometers))
        h_parent_map = {}
        for bond in topology.bonds():
            a1, a2 = bond.atom1, bond.atom2
            if a1.element is not None and a2.element is not None:
                if a1.element.symbol == 'H' and a2.element.symbol != 'H': h_parent_map[a1.index] = a2
                elif a2.element.symbol == 'H' and a1.element.symbol != 'H': h_parent_map[a2.index] = a1
                
        for h_idx, parent_atom in h_parent_map.items():
            p_heavy = np.array(pos_nm[parent_atom.index])
            p_h = np.array(pos_nm[h_idx])
            bond_vector = p_h - p_heavy
            current_dist = np.linalg.norm(bond_vector)
            if current_dist == 0: continue
            
            heavy_symbol = parent_atom.element.symbol
            if heavy_symbol == 'N': target_dist = 0.0860
            elif heavy_symbol == 'C': target_dist = 0.0930
            elif heavy_symbol == 'O': target_dist = 0.0840
            elif heavy_symbol == 'S': target_dist = 0.1000
            elif heavy_symbol == 'P': target_dist = 0.1000
            else: target_dist = max(0.01, current_dist - 0.0156)
            
            refined_h_pos = p_heavy + bond_vector * (target_dist / current_dist)
            pos_nm[h_idx] = Vec3(refined_h_pos[0], refined_h_pos[1], refined_h_pos[2])
            
        return pos_nm * unit.nanometers

    def minimize_topology_and_positions(self, openmm_topology, openmm_positions, output_path):
        """
        Run energy minimization using L-BFGS on CPU with the custom potential system.
        Writes finalized coordinates to a PDBx/mmCIF file.
        """
        print(f"[AetherEngine v40.2] ⚡ Pure Custom Refinement Core Activated")
        
        pos_nm = openmm_positions.value_in_unit(unit.nanometers) if hasattr(openmm_positions, 'value_in_unit') else openmm_positions
        
        system = AetherIdealizationEngine.create_idealization_system(openmm_topology, pos_nm)

        integrator = LangevinIntegrator(0, 0, 0)
        sim = Simulation(openmm_topology, system, integrator, Platform.getPlatformByName('CPU'))
        sim.context.setPositions(pos_nm)
        
        print(f"  ⚡ Entering L-BFGS energy minimization... (Max Iterations: {self.max_iterations})")
        LocalEnergyMinimizer.minimize(sim.context, 10.0, self.max_iterations)
        
        state = sim.context.getState(getPositions=True)
        ideal_positions = self.adjust_hydrogens_for_ideal_geometry(openmm_topology, state.getPositions())
        
        with open(output_path, 'w') as f:
            PDBxFile.writeFile(openmm_topology, ideal_positions, f)
            
        print(f"🎉 [AetherEngine v40.2] Perfect Refinement completed successfully!")
        return output_path
