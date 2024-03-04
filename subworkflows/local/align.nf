//
// Map to reference
//

include { ALIGN_BWAMEM2              } from './alignment/align_bwamem2'
include { ALIGN_SENTIEON             } from './alignment/align_sentieon'
include { SAMTOOLS_VIEW              } from '../../modules/nf-core/samtools/view/main'
include { ALIGN_MT                   } from './alignment/align_MT'
include { ALIGN_MT as ALIGN_MT_SHIFT } from './alignment/align_MT'
include { CONVERT_MT_BAM_TO_FASTQ    } from './mitochondria/convert_mt_bam_to_fastq'

workflow ALIGN {
    take:
        ch_reads                 // channel: [mandatory] [ val(meta), [path(reads)]  ]
        ch_genome_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_bwamem2index   // channel: [mandatory] [ val(meta), path(index) ]
        val_platform             // string:  [mandatory] illumina or a different technology

    main:
        ch_versions       = Channel.empty()
        ch_bwamem2_bam    = Channel.empty()
        ch_sentieon_bam   = Channel.empty()
        ch_bwamem2_bai    = Channel.empty()
        ch_sentieon_bai   = Channel.empty()

        ALIGN_BWAMEM2 (             // Triggered when params.aligner is set as bwamem2
            ch_reads,
            ch_genome_bwamem2index,
            ch_genome_fasta,
            ch_genome_fai,
            val_platform
            )
        ch_bwamem2_bam = ALIGN_BWAMEM2.out.marked_bam
        ch_bwamem2_bai = ALIGN_BWAMEM2.out.marked_bai
        ch_versions   = ch_versions.mix(ALIGN_BWAMEM2.out.versions)

        ch_genome_marked_bam = Channel.empty().mix(ch_bwamem2_bam)
        ch_genome_marked_bai = Channel.empty().mix(ch_bwamem2_bai)
        ch_genome_bam_bai    = ch_genome_marked_bam.join(ch_genome_marked_bai, failOnMismatch:true, failOnDuplicate:true)


    emit:
        genome_marked_bam  = ch_genome_marked_bam  // channel: [ val(meta), path(bam) ]
        genome_marked_bai  = ch_genome_marked_bai  // channel: [ val(meta), path(bai) ]
        genome_bam_bai     = ch_genome_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]
        versions           = ch_versions           // channel: [ path(versions.yml) ]
}
