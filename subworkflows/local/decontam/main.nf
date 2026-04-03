include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_DECONTAM } from '../../../modules/nf-core/bowtie2/align/main.nf'

workflow DECONTAM {
    take:
    reads // [[meta], [reads]]

    main:
    ch_multiqc_files = channel.empty()

    ch_bacterial_index = Channel.fromPath(["${params.bacterial_db}/*.bt2", "${params.bacterial_db}/*.bt2l"])
        .collect()
        .map { files -> [[id: 'bacterial_db'], files] }

    BOWTIE2_ALIGN_DECONTAM(
        reads,
        ch_bacterial_index,
        [[id: 'null'], []],
        true,   // save_unaligned: write unmapped reads as fastq.gz
        false   // view
    )
    ch_clean_reads = BOWTIE2_ALIGN_DECONTAM.out.fastq

    emit:
    reads    = ch_clean_reads // channel: [ val(meta), [ reads ] ]
    mqc      = ch_multiqc_files
}