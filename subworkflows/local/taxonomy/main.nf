include { BOWTIE2_ALIGN           as BOWTIE2_ALIGN_TAXONOMY     } from '../../../modules/nf-core/bowtie2/align/main'
include { TAXONOMY_Q30_FILTER                                   } from '../../../modules/local/taxonomy_q30_filter/main'
include { COVERAGEFILTER                                        } from '../../../modules/local/coveragefilter/main'
include { TAXONOMY_BAM_BUILD                                    } from '../../../modules/local/taxonomy_bam_build/main'
include { ASSEMBLY_HITS                                         } from '../../../modules/local/assembly_hits/main'
include { BUGLIST                                               } from '../../../modules/local/buglist/main'

workflow TAXONOMY {

    take:
    reads   // channel: [ val(meta), [ reads ] ]

    main:
    def format_version_yaml = { process, tool, version ->
        "\"${process}\":\n    ${tool}: ${version}"
    }

    ch_multiqc_files = channel.empty()

    ch_taxonomic_index = Channel.fromPath(
            ["${params.taxonomic_db}/*.bt2", "${params.taxonomic_db}/*.bt2l"])
        .collect()
        .map { files -> [[id: 'taxonomic_db'], files] }

    ch_taxonomy_txt = Channel.value(
        file("${params.taxonomic_db}/taxonomy.txt", checkIfExists: true))

    BOWTIE2_ALIGN_TAXONOMY(
        reads,
        ch_taxonomic_index,
        [[id: 'null'], []],
        false,  // save_unaligned
        false   // sort_bam → samtools view → SAM via args2 in modules.config
    )

    TAXONOMY_Q30_FILTER(
        BOWTIE2_ALIGN_TAXONOMY.out.sam
    )

    COVERAGEFILTER(TAXONOMY_Q30_FILTER.out.q30_sam)

    TAXONOMY_Q30_FILTER.out.q30_sorted_bam
        .join(COVERAGEFILTER.out.id_list, by: 0)
        .multiMap { meta, bam, id_list ->
            bam_ch:  [ meta, bam ]
            list_ch: id_list
        }
        .set { ch_bam_build_inputs }

    TAXONOMY_BAM_BUILD(
        ch_bam_build_inputs.bam_ch,
        ch_bam_build_inputs.list_ch
    )

    ASSEMBLY_HITS(TAXONOMY_BAM_BUILD.out.idxstats)

    BUGLIST(
        ASSEMBLY_HITS.out.hits,
        ch_taxonomy_txt
    )

    emit:
    buglist   = BUGLIST.out.stratified          // channel: [ val(meta), path(stratified_txt) ]
    mqc       = ch_multiqc_files
}
