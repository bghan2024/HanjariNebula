# -*- coding: utf-8 -*-
"""
HanjariNebulaEngine_GPU (OpenMM Edition - GPU Accelerated)
----------------------------------------------------------------------
[KOR] 구조생물학 파이프라인을 위한 고정밀 하이브리드 물리적 구조 정밀화(Refinement) 엔진입니다.
      OpenMM 기반으로 동작하며, CUDA 및 OpenCL을 통한 GPU 가속을 자동 지원합니다.
[ENG] A high-precision hybrid physical refinement engine designed for structural biology pipelines.
      Powered by OpenMM, featuring automatic GPU acceleration via CUDA and OpenCL.

Developer: Han Byeong-gu (hanbyeonggu@gmail.com)
Repository: https://github.com/bghan2024/HanjariNebula
Citation Ref: HanjariNebula Engine Suite (2026); see README.md for details.

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

class HanjariNebulaEngine_GPU:
    """
    HanjariNebulaEngine_GPU
    ----------------------
    [KOR] 구조적 정밀화를 수행하는 핵심 엔진 클래스입니다. 주요 기능은 다음과 같습니다.
    1. Perfect Index Mapping: 수소 추가 시 발생하는 원자 인덱스 뒤바뀜을 방지하기 위해 
       체인-잔기 키 맵을 생성하여 중원자(Heavy atoms)를 추적 및 고정합니다.
    2. Stereanalyst Patch: 단백질의 중원자에 조화 구속력(Harmonic Restraints)을 부여하여 
       최소화 중 공유 결합 손상이나 비물리적 구조 왜곡을 방지합니다.
    3. Phenix-Ready Hydrogen Adjustment: 페닉스(Phenix) 및 몰프로비티(MolProbity) 유효성 검증 표준에 
       부합하도록 수소 원자의 결합 길이를 조정합니다.
    4. mmCIF Quality Patches: 결과물 검증 오류를 예방하기 위해 mmCIF 파일의 점유율(Occupancy)과 
       온도 인자(B-factor) 메타데이터 포맷을 보정합니다.
    5. Adaptive GPU Acceleration: 시스템에 설치된 하드웨어를 감지하여 CUDA -> OpenCL -> CPU 순서로 
       최적의 컴퓨팅 플랫폼을 자동 할당합니다.

    [ENG] Main engine class that runs structural refinement. Key features:
    1. Perfect Index Mapping: Uses chain-residue mapping to prevent atom index scrambles.
    2. Stereanalyst Patch: Applies harmonic restraints on heavy atoms to prevent non-physical distortion.
    3. Phenix-Ready Hydrogen Adjustment: Re-scales hydrogen bonds to match Phenix chemistry standards.
    4. mmCIF Quality Patches: Re-formats occupancy/B-factor columns in mmCIF.
    5. Adaptive GPU Acceleration: Automatically selects CUDA, OpenCL, or CPU platform dynamically.
    """
    def __init__(self, max_iterations=300):
        self.max_iterations = max_iterations
        # Load Amber14 forcefield for proteins and TIP3P for water molecules
        self.forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')

    def _get_best_platform(self):
        """
        [KOR] 사용 가능한 최적의 컴퓨팅 플랫폼을 감지하여 반환합니다 (CUDA -> OpenCL -> CPU).
              GPU가 설치되어 있다면 CUDA나 OpenCL 플랫폼이 자동 선택되어 연산 속도가 극대화됩니다.
        [ENG] Detects and returns the best available computing platform (CUDA -> OpenCL -> CPU).
              If a compatible GPU is installed, CUDA or OpenCL is prioritized for maximum speed.
        """
        for platform_name in ['CUDA', 'OpenCL', 'CPU']:
            try:
                platform = Platform.getPlatformByName(platform_name)
                print(f"  🚀 [NebulaEngine] Active Compute Platform: {platform_name} (GPU/Hardware Accelerated)")
                return platform
            except Exception:
                continue
        print("  ⚠️ [NebulaEngine] No GPU platform found. Defaulting to CPU.")
        return Platform.getPlatformByName('CPU')

    def adjust_hydrogens_for_phenix(self, topology, positions):
        """
        [KOR] 수소 원자의 위치를 미세 조정하여 Phenix 및 MolProbity에서 권장하는 
              표준 기하학적 결합 길이(N-H, C-H, O-H, S-H)로 맞춥니다.
        [ENG] Dynamically adjusts hydrogen bond lengths to match standard chemistry bounds
              expected by Phenix and MolProbity.
        """
        pos_nm = list(positions.value_in_unit(unit.nanometers))
        h_parent_map = {}
        
        # Build map of hydrogen atoms and their parent heavy atoms
        for bond in topology.bonds():
            a1, a2 = bond.atom1, bond.atom2
            if a1.element is not None and a2.element is not None:
                if a1.element.symbol == 'H' and a2.element.symbol != 'H': 
                    h_parent_map[a1.index] = a2
                elif a2.element.symbol == 'H' and a1.element.symbol != 'H': 
                    h_parent_map[a2.index] = a1
                
        # Scale bonds to standard lengths
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
            
        return unit.Quantity(pos_nm, unit.nanometers)

    def minimize_structure(self, input_path, output_path):
        """
        [KOR] 입력된 구조 파일(PDB 또는 mmCIF)을 읽고, 수소 원자를 보강한 후 
              중원자 구속 조건을 가해 물리적 에너지 최소화(L-BFGS)를 실행합니다.
        [ENG] Reads structural input, adds hydrogens, applies heavy-atom restraints, 
              and executes L-BFGS energy minimization.
        """
        print(f"[NebulaEngine_GPU] Starting physical energy minimization... (Input: {os.path.basename(input_path)})")
        
        ext = os.path.splitext(input_path)[-1].lower()
        structure_obj = PDBxFile(input_path) if ext in ['.cif', '.mmcif'] else PDBFile(input_path)
        old_topology = structure_obj.topology

        # Check for nucleic acid residues and issue a warning
        # [KOR] 핵산 잔기 검출 및 forcefield 호환성 경고 출력
        nucleic_residues = {'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'RA', 'RC', 'RG', 'RU'}
        has_nucleic = any(r.name.strip().upper() in nucleic_residues for r in old_topology.residues())
        if has_nucleic:
            print("⚠️ [Warning] Nucleic acid residues (DNA/RNA) detected in the structure.")
            print("   The loaded Amber14 forcefield is optimized exclusively for proteins.")
            print("   Refining DNA/RNA structures may cause forcefield parameter matching failures.")
            print("   It is highly recommended to use this suite for protein-only refinement.")
            
        old_pos_nm = structure_obj.positions.value_in_unit(unit.nanometers)
        
        # Maps original heavy atoms to coordinates using chain and residue sequence index
        # [KOR] 기존 구조의 중원자 좌표를 체인 및 잔기 인덱스를 기준으로 매핑하여 저장
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
        
        # Resolve any chain ID duplicates (OpenMM strict requirement)
        used_chain_ids = set()
        available_letters = list(string.ascii_uppercase) + list(string.ascii_lowercase) + [f"C{i}" for i in range(1, 100)]
        letter_idx = 0
        
        # Rebuild topology excluding hydrogens (they will be freshly added later)
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

        # Re-add bonds to the new topology
        existing_new_bonds = set()
        for old_bond in old_topology.bonds():
            a1 = atom_mapping.get(old_bond.atom1)
            a2 = atom_mapping.get(old_bond.atom2)
            if a1 is not None and a2 is not None and a1 != a2:
                id1, id2 = min(a1.index, a2.index), max(a1.index, a2.index)
                if (id1, id2) not in existing_new_bonds:
                    new_topology.addBond(a1, a2)
                    existing_new_bonds.add((id1, id2))
                    
        # Construct Modeller object and add missing hydrogens
        modeller = Modeller(new_topology, unit.Quantity(new_positions, unit.nanometers))
        modeller.addHydrogens(forcefield=self.forcefield)
        
        # Create physical system
        system = self.forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.2, constraints=None)
        
        # Heavy-atom restraint force (Anti-Drift Force / Stereanalyst Patch)
        # [KOR] 중원자 위치 이탈 방지용 외부 고조파 복원력 설정
        restraint_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint_force.addGlobalParameter("k", 50.0 * unit.kilojoules_per_mole/unit.nanometer**2)
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")
        
        for final_chain_idx, final_chain in enumerate(modeller.topology.chains()):
            for final_res_idx, final_res in enumerate(final_chain.residues()):
                for final_atom in final_res.atoms():
                    if final_atom.element is not None and final_atom.element.symbol != 'H' and final_res.name != 'HOH':
                        match_key = (final_chain_idx, final_res_idx, final_atom.name)
                        if match_key in old_coords_map:
                            coord = old_coords_map[match_key]
                            # [KOR] Vec3 객체를 리스트 [x, y, z] 형태로 변환하여 OpenMM 입력을 통일
                            restraint_force.addParticle(final_atom.index, [coord.x, coord.y, coord.z])

        if restraint_force.getNumParticles() > 0: 
            system.addForce(restraint_force)

        # [KOR] LangevinIntegrator에 물리 단위를 명시하여 초기화 에러 방지 (300K, 1/ps 마찰, 0.002ps timestep)
        # [ENG] LangevinIntegrator initialized with proper physical units to prevent initialization error.
        integrator = LangevinIntegrator(300.0 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
        
        # ⚡ GPU / CPU Platform Dynamic Selection
        platform = self._get_best_platform()
        sim = Simulation(modeller.topology, system, integrator, platform)
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

        # Safe index fallback if column headers are not detected
        occ_idx = columns.get('_atom_site.occupancy', 13)
        b_idx = columns.get('_atom_site.B_iso_or_equiv', 14)
        
        patched_lines = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                tokens = line.split()
                # Ensure we don't cause IndexError by verifying token length
                if occ_idx < len(tokens): 
                    tokens[occ_idx] = '1.00'
                if b_idx < len(tokens): 
                    tokens[b_idx] = '20.00'
                line = " ".join(tokens) + "\n"
            patched_lines.append(line)
            
        with open(output_path, 'w') as f:
            f.writelines(patched_lines)
                
        print(f"[HanjariNebulaEngine_GPU] Physical refinement completed. Output file: {output_path}")
        return output_path

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("💡 Usage: python HanjariNebulaEngine_GPU.py <input.pdb/.cif> <output.cif>")
    else:
        # Default behavior to run independently
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else "refined_output.cif"
        engine = HanjariNebulaEngine_GPU()
        engine.minimize_structure(input_file, output_file)
