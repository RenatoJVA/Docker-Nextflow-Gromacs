# Docker-Nextflow-Gromacs


docker build -t gromacs-nextflow .

docker run -it --rm -v /work/MD_NEXTFLOW:/workspace gromacs-nextflow\n
docker run -it --rm -v ${PWD}:/workspace gromacs-nextflow


windows: 
docker run -it --rm -v ${PWD}:/workspace gromacs-nextflow
