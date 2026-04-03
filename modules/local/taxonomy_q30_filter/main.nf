process TAXONOMY_Q30_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.q30.sam"),        emit: q30_sam
    tuple val(meta), path("*.q30.sorted.bam"), emit: q30_sorted_bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -b -q 30 ${sam} \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.q30.sorted.bam

    samtools view -h -o ${prefix}.q30.sam ${prefix}.q30.sorted.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.q30.sam
    touch ${prefix}.q30.sorted.bam
    """
}
