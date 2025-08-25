nextflow.enable.dsl=2

params.pdb_id   = "1AKI"
params.water    = "spce"
params.output   = "results"
params.mdp_dir  = "${baseDir}/mdp"

// ==========================
// Procesos
// ==========================
process download_pdb {
    publishDir ".", mode: 'copy'

    output:
    tuple val(params.pdb_id), path("${params.pdb_id}")

    script:
    """
    mkdir -p ${params.pdb_id}
    wget https://files.rcsb.org/download/${params.pdb_id}.pdb
    grep -v HOH ${params.pdb_id}.pdb > ${params.pdb_id}_clean.pdb
    mv ${params.pdb_id}.pdb ${params.pdb_id}/
    mv ${params.pdb_id}_clean.pdb ${params.pdb_id}/
    # Copiar archivos mdp uno por uno
    for mdp in ions.mdp md.mdp minim.mdp npt.mdp nvt.mdp; do
        cp -L ${params.mdp_dir}/\$mdp ${params.pdb_id}/ || echo "Warning: Could not copy \$mdp"
    done
    """
}

process pdb2gmx {
    publishDir ".", mode: 'copy'
    
    input:
    tuple val(pdb_id), path(workdir)
    
    output:
    tuple val(pdb_id), path(workdir), path("${workdir}/${pdb_id}_processed.gro"), path("${workdir}/topol.top"), path("${workdir}/posre.itp")
    
    script:
    """
    # Nos aseguramos de estar en el directorio correcto
    cd ${workdir}
    
    # Ejecutamos pdb2gmx desde el directorio de trabajo y especificamos las rutas completas
    gmx pdb2gmx -f ${pdb_id}_clean.pdb -o ${pdb_id}_processed.gro -p topol.top -i posre.itp -water ${params.water} <<EOF
1
EOF
    """
}

process define_box {
    publishDir ".", mode: 'copy'

    input:
    tuple val(pdb_id), path(workdir), path(gro), path(top), path(itps)

    output:
    tuple val(pdb_id), path(workdir), path("${workdir}/boxed.gro"), path(top), path(itps)

    script:
    """
    cd ${workdir}
    gmx editconf -f ${gro.simpleName}.gro -o boxed.gro -c -d 1.0 -bt cubic
    """
}

process solvate {
    publishDir ".", mode: 'copy'

    input:
    tuple val(pdb_id), path(workdir), path(gro), path(top), path(itps)

    output:
    tuple val(pdb_id), path(workdir), path("${workdir}/solvated.gro"), path(top), path(itps)

    script:
    """
    cd ${workdir}
    gmx solvate -cp ${gro.name} -cs spc216.gro -o solvated.gro -p topol.top
    """
}

process ions {
    publishDir ".", mode: 'copy'

    input:
    tuple val(pdb_id), path(workdir), path(gro), path(top), path(itps)

    output:
    tuple val(pdb_id), path(workdir), path("${workdir}/ions.gro"), path(top), path(itps)

    script:
    """
    cd ${workdir}
    gmx grompp -f ions.mdp -c ${gro.name} -p topol.top -o ions.tpr
    echo "SOL" | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral
    """
}

process minimization {
    publishDir ".", mode: 'copy'

    input:
    tuple val(pdb_id), path(workdir), path(gro), path(top), path(itps)

    output:
    tuple val(pdb_id), path(workdir), path("${workdir}/em.gro"), path(top), path(itps)

    script:
    """
    cd ${workdir}
    gmx grompp -f minim.mdp -c ${gro.name} -p topol.top -o em.tpr
    gmx mdrun -deffnm em -nt 8
    """
}

process nvt {
    publishDir ".", mode: 'copy'

    input:
    tuple val(pdb_id), path(workdir), path(gro), path(top), path(itps)

    output:
    tuple val(pdb_id), path(workdir), path("${workdir}/nvt.gro"), path(top), path(itps)

    script:
    """
    cd ${workdir}
    gmx grompp -f nvt.mdp -c ${gro.name} -r ${gro.name} -p topol.top -o nvt.tpr
    gmx mdrun -nice 0 -v -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu -pin on -pinoffset 0 -pinstride 1
    """
}

process npt {
    publishDir ".", mode: 'copy'

    input:
    tuple val(pdb_id), path(workdir), path(gro), path(top), path(itps)

    output:
    tuple val(pdb_id), path(workdir), path("${workdir}/npt.gro"), path(top), path(itps)

    script:
    """
    cd ${workdir}
    gmx grompp -f npt.mdp -c ${gro.name} -r ${gro.name} -p topol.top -o npt.tpr
    gmx mdrun -nice 0 -v -deffnm npt -nt 8 -ntmpi 1 -ntomp 8 -nb gpu -pin on -pinoffset 0 -pinstride 1
    """
}

process md {
    publishDir ".", mode: 'copy'

    input:
    tuple val(pdb_id), path(workdir), path(gro), path(top), path(itps)

    output:
    tuple path("${workdir}/md.gro"), path("${workdir}/md.edr"), path("${workdir}/md.log"), path("${workdir}/md.trr")

    script:
    """
    cd ${workdir}
    gmx grompp -f md.mdp -c ${gro.name} -r ${gro.name} -t ${gro.name} -p topol.top -o md.tpr
    gmx mdrun -nice 0 -v -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pin on -pinoffset 0 -pinstride 1
    """
}

process copy_results {
    input:
    path("*.gro")
    path("*.log")
    path("*.edr")
    path("*.trr")

    output:
    path("${params.output}")

    script:
    """
    mkdir -p ${params.output}
    cp *.gro *.log *.edr *.trr ${params.output}/
    """
}

// ==========================
// Workflow
// ==========================
workflow {
    ch_setup = download_pdb()
    
    // Todos los procesos trabajarán en el directorio de la proteína
    ch_gmx   = pdb2gmx(ch_setup)
    ch_box   = define_box(ch_gmx)
    ch_solv  = solvate(ch_box)
    ch_ions  = ions(ch_solv)
    ch_em    = minimization(ch_ions)
    ch_nvt   = nvt(ch_em)
    ch_npt   = npt(ch_nvt)
    md(ch_npt)
}
