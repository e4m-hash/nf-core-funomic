process GET_PWY_ABUNDANCE {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.2.1 bioconda::bioconductor-keggrest=1.44.0 conda-forge::r-base=4.3.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://amancevice/pandas:2.2.1'
        : 'docker.io/amancevice/pandas:2.2.1'}"

    input:
    tuple val(meta), path(filtered_blast)
    path(protein_db)        // directory: jgi_ann, ncbi_uniprot_species, taxonomy_for_function

    output:
    tuple val(meta), path("*_pwy_abundance.csv"),      emit: pwy_abundance
    tuple val(meta), path("*_pwyClass_abundance.csv"),  emit: pwy_class_abundance
    tuple val(meta), path("*_pwyType_abundance.csv"),   emit: pwy_type_abundance
    tuple val(meta), path("*_full_annotation.csv"),     emit: full_annotation
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_pwy_abundance.py \\
        -i  ${filtered_blast} \\
        -s  ${protein_db} \\
        -o1 ${prefix}_pwy_abundance.csv \\
        -o2 ${prefix}_pwyClass_abundance.csv \\
        -o3 ${prefix}_pwyType_abundance.csv \\
        -o4 ${prefix}_full_annotation.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pwy_abundance.csv
    touch ${prefix}_pwyClass_abundance.csv
    touch ${prefix}_pwyType_abundance.csv
    touch ${prefix}_full_annotation.csv
    """
}
