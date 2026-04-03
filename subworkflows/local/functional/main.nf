include { MERGE_READS       } from '../../../modules/local/merge_reads/main'
include { DIAMOND_BLASTX    } from '../../../modules/nf-core/diamond/blastx/main'
include { DIAMOND_FILTER    } from '../../../modules/local/diamond_filter/main'
include { GET_PWY_ABUNDANCE } from '../../../modules/local/get_pwy_abundance/main'

workflow FUNCTIONAL {

    take:
    reads   // channel: [ val(meta), [ reads ] ]  (decontaminated)

    main:
    ch_multiqc_files = channel.empty()

    // FunOMIC-P protein database directory
    // Contains: FunOMIC.P.v2_nr95.dmnd, jgi_ann_05-2022_reordered.tab,
    //           ncbi_uniprot_species.txt, taxonomy_for_function.csv
    ch_protein_db_dir  = Channel.value(
        file(params.protein_db, checkIfExists: true))

    // DIAMOND database file (.dmnd) for BLASTX alignment
    ch_dmnd_db = Channel.value(
        [[id: 'protein_db'],
         file("${params.protein_db}/FunOMIC.P.v2_nr95.dmnd", checkIfExists: true)])

    MERGE_READS(reads)

    DIAMOND_BLASTX(
        MERGE_READS.out.fastq,
        ch_dmnd_db,
        'txt',  // output format extension (tabular BLAST format 6)
        []      // no extra blast columns
    )

    DIAMOND_FILTER(DIAMOND_BLASTX.out.txt)

    GET_PWY_ABUNDANCE(
        DIAMOND_FILTER.out.filtered,
        ch_protein_db_dir
    )

    emit:
    pwy_abundance       = GET_PWY_ABUNDANCE.out.pwy_abundance       // [ meta, pwy_abundance.csv ]
    pwy_class_abundance = GET_PWY_ABUNDANCE.out.pwy_class_abundance // [ meta, pwyClass_abundance.csv ]
    pwy_type_abundance  = GET_PWY_ABUNDANCE.out.pwy_type_abundance  // [ meta, pwyType_abundance.csv ]
    full_annotation     = GET_PWY_ABUNDANCE.out.full_annotation      // [ meta, full_annotation.csv ]
    mqc                 = ch_multiqc_files
}
