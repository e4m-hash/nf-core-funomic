process MERGE_READS {
    tag "$meta.id"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.joined.fastq"), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Decompress and concatenate all reads (R1 + R2 for paired-end,
    # single file for single-end) into one FASTQ for DIAMOND blastx
    cat <(zcat ${reads[0]} 2>/dev/null || cat ${reads[0]}) \
        <(zcat ${reads[1]} 2>/dev/null || cat ${reads[1]}) \
        > ${prefix}.joined.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.joined.fastq
    """
}
