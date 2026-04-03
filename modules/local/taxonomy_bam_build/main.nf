process TAXONOMY_BAM_BUILD {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(sorted_bam)
    path(id_list)

    output:
    tuple val(meta), path("*.idxstats.txt"), emit: idxstats
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract reads whose QNAME is in the coverage-pass list
    # grep may return exit 1 (no matches) → handle gracefully
    samtools view -h -N ${id_list} ${sorted_bam} > ${prefix}.c80.sam

    if [ -s ${prefix}.c80.sam ]; then

        samtools view -b -o ${prefix}.c80.bam ${prefix}.c80.sam
        samtools index ${prefix}.c80.bam
        samtools idxstats \\
            --threads ${task.cpus} \\
            ${prefix}.c80.bam \\
            > ${prefix}.idxstats.txt
    else
        # No coverage-passing reads: produce empty idxstats
        touch ${prefix}.idxstats.txt
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idxstats.txt
    """
}
