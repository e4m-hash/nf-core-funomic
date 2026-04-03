process BUGLIST {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.2.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://amancevice/pandas:2.2.1'
        : 'docker.io/amancevice/pandas:2.2.1'}"

    input:
    tuple val(meta), path(assembly_hits)
    path(taxonomy)

    output:
    tuple val(meta), path("*_buglist_stratified.txt"), emit: stratified
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    buglist.py \\
        -i ${assembly_hits} \\
        -t ${taxonomy} \\
        -o ${prefix}_buglist_stratified.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_buglist_stratified.txt
    """
}
