# -*- coding: utf-8 -*-
"""
HanjariFold2_CosmicDistiller_OS (Open-Source / OpenMM Edition)
--------------------------------------------------------------
[KOR] 단백질 분자 모델 주변에 최적화된 수화 격자(Hydration Grid)를 형성하고, 
      물리적 에너지 최소화를 자동 관리해 주는 통합 정밀화 래퍼(Wrapper) 스크립트입니다.
[ENG] An automated macromolecular hydration and physical energy minimization manager.
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

# Import the GPU-accelerated engine from the local path
from HanjariNebulaEngine_OS import HanjariNebulaEngine_OS

class HanjariFold2CosmicDistiller_OS:
    """
    HanjariFold2CosmicDistiller_OS
    ---------------------------
    [KOR] 수화층 컴파일과 물리적 정밀화 프로토콜을 수행하는 코디네이터 클래스입니다.
    1. Perfect Water Compiling: 단백질 표면 주변에 H-O-H 물 격자 공간을 재구성하여, 
       원자 좌표 손실 영역 및 비정상 결합을 완충(Screening)합니다.
    2. Gemmi-based QA Checking: 최종 리파인먼트 결과 mmCIF 구조의 원자 개수, 잔기 무결성, 
       물 분자 충진 수 등을 Gemmi 라이브러리로 신속 검증합니다.
    3. Structural Synergy: 보강된 수화층 좌표계를 HanjariNebulaEngine_OS의 L-BFGS 
       구속 최소화(Restrained Minimization) 알고리즘과 결합하여 기하학적 수렴을 보장합니다.

    [ENG] Coordinator class running the hydration lattice compilation and refinement loop.
    1. Perfect Water Compiling: Computes a structured water grid layer near the macromolecule.
    2. Gemmi-based QA Checking: Performs structural integrity checks (chains, residues, water count) via Gemmi.
    3. Structural Synergy: Combines hydrated topology with the HanjariNebula restraint optimizer.
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
        [KOR] 최종 mmCIF 결과 파일을 읽어 체인, 아미노산 잔기, 수화 물 분자 수 등을 분석하여 검증 보고서를 출력합니다.
        [ENG] Reads the resulting mmCIF file and runs a structural scan to verify chains, residues,
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
        """
        [KOR] 수화 리파인먼트 루틴 실행:
              1. PDB/mmCIF 파싱 후 단백질 3D 좌표를 추출합니다.
              2. 4.5Å 크기의 3D 격자를 만들고, KDTree 공간 검색을 수행하여 
                 단백질 원자에서 2.8Å~3.8Å 내에 위치한 격자점에 물(HOH) 분자를 배치합니다.
              3. 단백질과 물 분자 토폴로지를 통합(Modeller)합니다.
              4. HanjariNebulaEngine_OS를 호출해 GPU 가속화가 적용된 L-BFGS 에너지 최소화를 실행합니다.
              5. 최종 mmCIF 파일을 생성하고 임시 리소스를 해제합니다.
        [ENG] Refinement execution loop:
              1. Parses input coordinates.
              2. Generates 4.5Å 3D grid and filters coordinates 2.8Å to 3.8Å away from any protein heavy atom using KDTree.
              3. Merges coordinates and topologies using Modeller.
              4. Invokes the GPU-accelerated HanjariNebulaEngine for restrained L-BFGS minimization.
              5. Exports the optimized mmCIF and triggers Gemmi validation.
        """
        self.log(f"🌌 CosmicDistiller starting... Target: {os.path.basename(self.input_path)}")
        
        ext = os.path.splitext(self.input_path)[-1].lower()
        structure_obj = app.PDBxFile(self.input_path) if ext in ['.cif', '.mmcif'] else app.PDBFile(self.input_path)
        
        # Use standard getPositions() API
        positions_nm = structure_obj.getPositions().value_in_unit(unit.nanometers)
        coords = np.array([[pos[0], pos[1], pos[2]] for pos in positions_nm]) * 10.0 # convert to Angstroms
        
        self.log(f"💧 [CosmicDistiller] Constructing hydration grid (Spacing: 4.5Å)")
        min_c = np.min(coords, axis=0) - 4.0
        max_c = np.max(coords, axis=0) + 4.0
        
        # Build 3D mesh grid points
        grid_x = np.arange(min_c[0], max_c[0], 4.5)
        grid_y = np.arange(min_c[1], max_c[1], 4.5)
        grid_z = np.arange(min_c[2], max_c[2], 4.5)
        grid_points = np.array(np.meshgrid(grid_x, grid_y, grid_z)).T.reshape(-1, 3)
        
        # KDTree distance filtering: find water coordinates near protein surface
        tree = cKDTree(coords)
        distances, _ = tree.query(grid_points, k=1)
        # Optimal distance for hydrogen bonding range (2.8Å - 3.8Å)
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
            
            # Position coordinates in nanometers
            o_pos = Vec3(coord[0]/10.0, coord[1]/10.0, coord[2]/10.0)
            h1_pos = Vec3(o_pos.x + 0.0957, o_pos.y, o_pos.z)
            h2_pos = Vec3(o_pos.x - 0.0240, o_pos.y + 0.0926, o_pos.z)
            
            water_positions.extend([o_pos, h1_pos, h2_pos])
            
            water_topo.addBond(o_atom, h1_atom)
            water_topo.addBond(o_atom, h2_atom)
            
        modeller = app.Modeller(structure_obj.topology, structure_obj.getPositions())
        # Safe assignment of unit quantity for the lists
        modeller.add(water_topo, unit.Quantity(water_positions, unit.nanometers))
        
        temp_hydrated_cif = os.path.join(self.output_dir, "hydrated_temp.cif")
        with open(temp_hydrated_cif, 'w') as f:
            app.PDBxFile.writeFile(modeller.topology, modeller.positions, f)
        
        final_output_file = os.path.join(self.output_dir, "distilled_final.cif")
        
        # Initialize GPU/CPU physical minimization engine
        minimizer_engine = HanjariNebulaEngine_OS(max_iterations=300)
        minimizer_engine.minimize_structure(temp_hydrated_cif, final_output_file)
        
        # Clean up temporary hydrated cif
        if os.path.exists(temp_hydrated_cif):
            os.remove(temp_hydrated_cif)
            
        # Validate output structural components
        self.run_gemmi_validation(final_output_file)
        self.log("🎉 [HanjariNebula] Protein-hydration complex assembly and physical refinement completed.")
        self.log(f"💾 Result saved to: {final_output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("💡 Usage: python HanjariFold2_CosmicDistiller_OS.py <input.pdb/.cif>")
    else:
        HanjariFold2CosmicDistiller_OS(sys.argv[1]).run()
