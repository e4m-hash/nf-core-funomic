process ASSEMBLY_HITS {
    tag "$meta.id"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(idxstats)

    output:
    tuple val(meta), path("*_assembly_hits.txt"), emit: hits

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F '\\t' '\$3 > 0 {
        split(\$1, id_parts, "|");
        print id_parts[1] "\\t" \$3
    }' ${idxstats} > ${prefix}_assembly_hits.txt

    touch ${prefix}_assembly_hits.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_assembly_hits.txt
    """
}