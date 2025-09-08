# Docker-Nextflow-GROMACS Pipeline

A containerized workflow for running molecular dynamics simulations using GROMACS through Nextflow pipeline. This project combines Docker, Nextflow, and GROMACS to create a reproducible and scalable environment for molecular dynamics simulations.

## Prerequisites

- Docker
- NVIDIA GPU with CUDA support
- NVIDIA Container Toolkit
- Git

## Features

- GROMACS 2024.4 with CUDA support
- Nextflow workflow for automated MD simulations
- Both single and double precision GROMACS installations
- Complete MD pipeline including:
  - Structure preparation
  - System solvation
  - Energy minimization
  - NVT equilibration
  - NPT equilibration
  - Production MD

## Quick Start

1. Clone the repository:
```bash
git clone https://github.com/RenatoJVA/Docker-Nextflow-Gromacs.git
cd Docker-Nextflow-Gromacs
```

2. Build the Docker image:
```bash
docker build \
  --build-arg USER_UID=$(id -u) \
  --build-arg USER_GID=$(id -g) \
  --build-arg USERNAME=$(whoami) \
  -t gmx .
```

3. Run the container:
```bash
docker run --gpus all -it --rm -v ${PWD}:/workspace gmx
```

4. Inside the container, run the Nextflow pipeline:
```bash
nextflow run main.nf
```

## Pipeline Parameters

You can customize the simulation by modifying these parameters:

- `pdb_id`: PDB ID to download (default: "1AKI")
- `water`: Water model to use (default: "spce")
- `output`: Output directory name (default: "results")
- `mdp_dir`: Directory containing MDP parameter files (default: "./mdp")

Example with custom parameters:
```bash
nextflow run main.nf --pdb_id "2LYZ" --water "tip3p" --output "lysozyme_results"
```

## Pipeline Steps

1. Download PDB structure
2. Process structure with pdb2gmx
3. Define simulation box
4. Solvate the system
5. Add ions for neutralization
6. Energy minimization
7. NVT equilibration
8. NPT equilibration
9. Production MD simulation

## Directory Structure

```
.
├── dockerfile         # Docker configuration file
├── main.nf           # Nextflow pipeline script
├── README.md         # This file
└── mdp/              # GROMACS parameter files
    ├── ions.mdp      # Ion addition parameters
    ├── md.mdp        # Production MD parameters
    ├── minim.mdp     # Energy minimization parameters
    ├── npt.mdp       # NPT equilibration parameters
    └── nvt.mdp       # NVT equilibration parameters
```

## Hardware Requirements

- NVIDIA GPU with CUDA support
- Minimum 16GB RAM recommended
- Storage space depending on simulation size

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

If you encounter any issues or have questions, please open an issue in the GitHub repository.