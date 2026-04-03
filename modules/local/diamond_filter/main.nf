process DIAMOND_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(blast_out)

    output:
    tuple val(meta), path("*.func.filtered.out"), emit: filtered

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    sed 's/[[:space:]]\\+/\\t/g' "${blast_out}" | \
    awk -F'\\t' '\$4>=20' | \
    LC_ALL=C sort -t\$'\\t' -k1,1 -k12,12nr | \
    awk -F'\\t' '!seen[\$1]++' > "${prefix}.func.filtered.out" || true
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch "${prefix}.func.filtered.out"
    """
}