# HanjariNebula Engine Suite

[![Compute](https://img.shields.io/badge/GPU%20Acceleration-CUDA%20%7C%20OpenCL-success.svg)](#gpu-acceleration-details)

HanjariNebula Engine is a high-precision macromolecular physical refinement and hydration compilation framework designed for modern structural biology pipelines. It utilizes **OpenMM** to reconstruct structured hydration grids and minimize energy states using restrained mechanics.

---

## Key Features

1. **Perfect Water Compiling (Hydration Lattice)**
   - Reconstructs a structured $H$-$O$-$H$ water grid layer around target proteins.
   - Resolves atom coordinate vacancies and acts as a buffer grid to smooth physical constraints.
2. **Stereanalyst Patch (Harmonic Heavy-Atom Restraints)**
   - Protects protein backbone and sidechain structures from non-physical distortions during minimization using harmonic force restraints ($k = 50.0\text{ kJ/mol/nm}^2$).
3. **Phenix-Ready Hydrogen Adjustment**
   - Re-scales and standardizes $N$-$H$, $C$-$H$, $O$-$H$, and $S$-$H$ bonds to meet MolProbity and Phenix downstream geometry validation standards.
4. **mmCIF Compatibility Quality Patches**
   - Corrects occupancy (to 1.00) and B-factors (to 20.00) in output mmCIF data to avoid validation workflow interruptions.
5. **Gemmi-Based Quality Assessment**
   - Verifies the final structure's chain definitions, residue counts, and solvated water ratios dynamically.

---

## GPU Acceleration Details

HanjariNebula Engine dynamically scans your system's hardware environment and selects the most optimal computation platform:

$$\text{Computation Priority: } \mathbf{CUDA} \longrightarrow \mathbf{OpenCL} \longrightarrow \mathbf{CPU}$$

- **CUDA Platform (Highly Recommended)**: Automatically activated if an NVIDIA GPU driver and compatible CUDA Toolkit are present. This provides massive parallel computing gains for physical molecular energy evaluation.
- **OpenCL Platform**: Activated on systems with AMD, Intel, or other GPU hardware supporting OpenCL standards.
- **CPU Platform**: Fallback platform used if no GPU drivers or acceleration runtime libraries are detected.

---

## Installation

We recommend using `conda` or `mamba` to manage your environment, as OpenMM and Gemmi are best installed via the `conda-forge` channel.

```bash
# 1. Create a virtual environment
conda create -n hanjari_nebula python=3.10 -y
conda activate hanjari_nebula

# 2. Install OpenMM and core dependencies
conda install -c conda-forge openmm scipy numpy gemmi -y

# 3. Install HanjariNebula Engine directly from GitHub
pip install git+[https://github.com/bghan2024/HanjariNebula.git](https://github.com/bghan2024/HanjariNebula.git)
(Note: Verify your GPU drivers and CUDA toolkit are correctly configured to enable CUDA/OpenCL acceleration.)

---

## Usage

You can run either the molecular refinement engine independently or the entire hydration-minimization pipeline.

### 1. Run Complete Hydration & Refinement Pipeline
This automatically builds the hydration grid (water layer) around your protein and minimizes the system:
```bash
python HanjariFold2_CosmicDistiller_GPU.py <input_protein.pdb_or_cif>
```
- **Output**: Generates a run directory (e.g., `HF2_NEBULA_DISTILLER_YYYYMMDD_HHMMSS/`) containing:
  - `distilled_final.cif`: The final, physically optimized protein-hydration complex.

### 2. Run Restrained Minimization Only
To perform L-BFGS energy minimization with heavy-atom restraints without adding a water layer:
```bash
python HanjariNebulaEngine_GPU.py <input.pdb_or_cif> <output.cif>
```

---

## Critical Usage Constraints

> [!WARNING]
> Due to parametrization constraints within the default Amber14 libraries (`amber14/protein.ff14SB.xml` and `amber14/tip3p.xml`), **this engine is optimized exclusively for proteins**. 
> Refining nucleic acid coordinates (DNA/RNA) is highly discouraged as it may trigger topology-matching or parameter missing failures. Please limit inputs to protein-only structures.

---

## Developers & Citation

- **Developer**: Han Byeong-gu ([hanbyeonggu@gmail.com](mailto:hanbyeonggu@gmail.com))
- **Repository**: [https://github.com/bghan2024/HanjariNebula](https://github.com/bghan2024/HanjariNebula)

### Citation
```text
HanjariNebula Engine Suite (2026); GitHub Repository: https://github.com/bghan2024/HanjariNebula
```

---

## License

- **Academic & Non-Commercial Research Use**: Granted Free of charge.
- **Commercial & For-Profit Use**: Strictly requires a separate written commercial contract. Please reach out to [hanbyeonggu@gmail.com](mailto:hanbyeonggu@gmail.com) for inquiries.
