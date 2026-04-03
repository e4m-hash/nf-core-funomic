process COVERAGEFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.2.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://amancevice/pandas:2.2.1'
        : 'docker.io/amancevice/pandas:2.2.1'}"

    input:
    tuple val(meta), path(q30_sam)

    output:
    tuple val(meta), path("*.c80.list"), emit: id_list
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ""
    
    """
    coverageFilter.py \\
        -i ${q30_sam} \\
        -o ${prefix}.c80.list \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.c80.list
    """
}
