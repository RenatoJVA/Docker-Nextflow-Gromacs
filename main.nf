nextflow.enable.dsl=2

params.pdb_id   = "1AKI"
params.water    = "spce"
params.output   = "results"
params.mdp      = "mdp"

// ==========================
// Procesos
// ==========================
process download_pdb {
    publishDir "input", mode: 'copy'

    output:
    path("${params.pdb_id}.pdb")
    path("${params.pdb_id}_clean.pdb")

    script:
    """
    wget https://files.rcsb.org/download/${params.pdb_id}.pdb -O ${params.pdb_id}.pdb
    grep -v HOH ${params.pdb_id}.pdb > ${params.pdb_id}_clean.pdb
    """
}

process pdb2gmx {
    input:
    path pdb_file

    output:
    tuple path("${params.pdb_id}_processed.gro"), path("topol.top"), path("*.itp")

    script:
    """
    gmx pdb2gmx -f $pdb_file -o ${params.pdb_id}_processed.gro -water ${params.water} <<EOF
1
EOF
    """
}

process define_box {
    input:
    tuple path(gro), path(top), path(itps)

    output:
    tuple path("boxed.gro"), path(top), path(itps)

    script:
    """
    gmx editconf -f $gro -o boxed.gro -c -d 1.0 -bt cubic
    """
}

process solvate {
    input:
    tuple path(gro), path(top), path(itps)

    output:
    tuple path("solvated.gro"), path(top), path(itps)

    script:
    """
    gmx solvate -cp $gro -cs spc216.gro -o solvated.gro -p $top
    """
}

process ions {
    input:
    tuple path(gro), path(top), path(itps)
    path ions_mdp = "${params.mdp}/ions.mdp"

    output:
    tuple path("ions.gro"), path(top), path(itps)

    script:
    """
    gmx grompp -f $ions_mdp -c $gro -p $top -o ions.tpr
    echo "SOL" | gmx genion -s ions.tpr -o ions.gro -p $top -pname NA -nname CL -neutral
    """
}

process minimization {
    input:
    tuple path(gro), path(top), path(itps)
    path minim_mdp = "${params.mdp}/minim.mdp"

    output:
    tuple path("em.gro"), path(top), path(itps)

    script:
    """
    gmx grompp -f $minim_mdp -c $gro -p $top -o em.tpr
    gmx mdrun -deffnm em
    """
}

process nvt {
    input:
    tuple path(gro), path(top), path(itps)
    path nvt_mdp = "${params.mdp}/nvt.mdp"

    output:
    tuple path("nvt.gro"), path(top), path(itps)

    script:
    """
    gmx grompp -f $nvt_mdp -c $gro -p $top -o nvt.tpr
    gmx mdrun -deffnm nvt
    """
}

process npt {
    input:
    tuple path(gro), path(top), path(itps)
    path npt_mdp = "${params.mdp}/npt.mdp"

    output:
    tuple path("npt.gro"), path(top), path(itps)

    script:
    """
    gmx grompp -f $npt_mdp -c $gro -p $top -o npt.tpr
    gmx mdrun -deffnm npt
    """
}

process md {
    input:
    tuple path(gro), path(top), path(itps)
    path md_mdp = "${params.mdp}/md.mdp"

    output:
    path("md.gro")
    path("md.edr")
    path("md.log")
    path("md.trr")

    script:
    """
    gmx grompp -f $md_mdp -c $gro -p $top -o md.tpr
    gmx mdrun -deffnm md
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
    ch_pdb   = download_pdb()
    ch_gmx   = pdb2gmx(ch_pdb[1])
    ch_box   = define_box(ch_gmx)
    ch_solv  = solvate(ch_box)
    ch_ions  = ions(ch_solv)
    ch_em    = minimization(ch_ions)
    ch_nvt   = nvt(ch_em)
    ch_npt   = npt(ch_nvt)
    ch_md    = md(ch_npt)
    copy_results(ch_md)
}
