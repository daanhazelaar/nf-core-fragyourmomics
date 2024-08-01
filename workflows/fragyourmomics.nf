/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_fragyourmomics_pipeline'

include { PREPARE_REFERENCE_GENOME  } from '../subworkflows/local/prepare_reference_genome'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { BEDTOOLS_INTERSECT } from '../modules/nf-core/bedtools/intersect/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FRAGYOURMOMICS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    PREPARE_REFERENCE_GENOME()


    SAMTOOLS_VIEW(
        ch_samplesheet.map{meta, bam, bai, assay -> return[ meta, bam, bai ]},
        PREPARE_REFERENCE_GENOME.out.fasta.map{[ [:], it]},
        []
    )

    // MODULE: SAMTOOLS_INDEX
    SAMTOOLS_SORT (
        SAMTOOLS_VIEW.out.bam,
        PREPARE_REFERENCE_GENOME.out.fasta.map{[ [:], it]}
    )

    // MODULE: SAMTOOLS_INDEX
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )

    // MODULE: BEDTOOLS_INTERSECT
    BEDTOOLS_INTERSECT(
        SAMTOOLS_SORT.out.bam
            .combine(Channel.value(file(params.bed))),
        PREPARE_REFERENCE_GENOME.out.chrom_sizes.map{[ [:], it]}
    )




    emit:
    multiqc_report = Channel.empty()             // channel: /path/to/multiqc_report.html
    versions       = Channel.empty()             // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
