/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRaredisease.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK MANDATORY PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def mandatoryParams = [
    "aligner",
    "analysis_type",
    "fasta",
    "input",
    "intervals_wgs",
    "intervals_y",
    "platform",
    "variant_catalog",
    "variant_caller"
]
def missingParamsCount = 0


if (missingParamsCount>0) {
    error("\nSet missing parameters and restart the run. For more information please check usage documentation on github.")
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config              ? Channel.fromPath( params.multiqc_config )  : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo                ? Channel.fromPath( params.multiqc_logo )    : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description )  : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { SMNCOPYNUMBERCALLER                   } from '../modules/nf-core/smncopynumbercaller/main'
include { ENSEMBLVEP_FILTERVEP as FILTERVEP_MT } from '../modules/nf-core/ensemblvep/filtervep'
include { ENSEMBLVEP_FILTERVEP as FILTERVEP_SNV } from '../modules/nf-core/ensemblvep/filtervep'
include { ENSEMBLVEP_FILTERVEP as FILTERVEP_SV  } from '../modules/nf-core/ensemblvep/filtervep'
include { TABIX_BGZIPTABIX as BGZIPTABIX_MT    } from '../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIPTABIX_SNV    } from '../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIPTABIX_SV     } from '../modules/nf-core/tabix/bgziptabix'

//
// SUBWORKFLOWS
//

include { ALIGN                                 } from '../subworkflows/local/align'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_MT    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV   } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SV    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_GENOME_SNVS                  } from '../subworkflows/local/annotate_genome_snvs'
include { ANNOTATE_MT_SNVS                      } from '../subworkflows/local/annotate_mt_snvs'
include { ANNOTATE_STRUCTURAL_VARIANTS          } from '../subworkflows/local/annotate_structural_variants'
include { CALL_REPEAT_EXPANSIONS                } from '../subworkflows/local/call_repeat_expansions'
include { CALL_SNV                              } from '../subworkflows/local/call_snv'
include { CALL_STRUCTURAL_VARIANTS              } from '../subworkflows/local/call_structural_variants'
include { GENS                                  } from '../subworkflows/local/gens'
include { PREPARE_REFERENCES                    } from '../subworkflows/local/prepare_references'
include { QC_BAM                                } from '../subworkflows/local/qc_bam'
include { RANK_VARIANTS as RANK_VARIANTS_MT     } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SNV    } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SV     } from '../subworkflows/local/rank_variants'
include { SCATTER_GENOME                        } from '../subworkflows/local/scatter_genome'
include { PEDDY_CHECK                           } from '../subworkflows/local/peddy_check'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAREDISEASE {

    ch_versions = Channel.empty()

    // Initialize read, sample, and case_info channels
    ch_input = Channel.fromPath(params.input)
    Channel.fromSamplesheet("input")
        .tap { ch_original_input }
        .map { meta, fastq1, fastq2 -> meta.id }
        .reduce([:]) { counts, sample -> //get counts of each sample in the samplesheet - for groupTuple
            counts[sample] = (counts[sample] ?: 0) + 1
            counts
        }
        .combine( ch_original_input )
        .map { counts, meta, fastq1, fastq2 ->
            new_meta = meta + [num_lanes:counts[meta.id],
                        read_group:"\'@RG\\tID:"+ fastq1.toString().split('/')[-1] + "\\tPL:ILLUMINA\\tSM:"+meta.id+"\'"]
            if (!fastq2) {
                return [ new_meta + [ single_end:true ], [ fastq1 ] ]
            } else {
                return [ new_meta + [ single_end:false ], [ fastq1, fastq2 ] ]
            }
        }
        .tap{ ch_input_counts }
        .map { meta, fastqs -> fastqs }
        .reduce([:]) { counts, fastqs -> //get line number for each row to construct unique sample ids
            counts[fastqs] = counts.size() + 1
            return counts
        }
        .combine( ch_input_counts )
        .map { lineno, meta, fastqs -> //append line number to sampleid
            new_meta = meta + [id:meta.id+"_T"+lineno[fastqs]]
            return [ new_meta, fastqs ]
        }
        .set { ch_reads }

    ch_samples   = ch_reads.map { meta, fastqs -> meta}
    ch_pedfile   = ch_samples.toList().map { makePed(it) }
    ch_case_info = ch_samples.toList().map { create_case_channel(it) }

    // Initialize file channels for PREPARE_REFERENCES subworkflow
    ch_genome_fasta             = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
    ch_genome_fai               = params.fai            ? Channel.fromPath(params.fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                        : Channel.empty()

    // Prepare references and indices.
    PREPARE_REFERENCES (
        ch_genome_fasta
    )
    .set { ch_references }

    
    ch_genome_bwamem2index      = params.bwamem2                           ? Channel.fromPath(params.bwamem2).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.genome_bwamem2_index
    ch_versions                 = ch_versions.mix(ch_references.versions)

    //
    // ALIGNING READS, FETCH STATS, AND MERGE.
    //
    ALIGN (
        ch_reads,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_bwamem2index,
        params.platform
    )
    .set { ch_mapped }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)

    //
    // MODULE: Pipeline reporting
    //

    // The template v2.7.1 template update introduced: ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    // This caused the pipeline to stall
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRaredisease.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRaredisease.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    if (!params.skip_fastqc) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.hs_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.ped.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.csv.map{it[1]}.collect().ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def makePed(samples) {

    def case_name  = samples[0].case_id
    def outfile  = file("${params.outdir}/pipeline_info/${case_name}" + '.ped')
    outfile.text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\t')
    def samples_list = []
    for(int i = 0; i<samples.size(); i++) {
        sample_name        =  samples[i].sample
        if (!samples_list.contains(sample_name)) {
            outfile.append('\n' + [samples[i].case_id, sample_name, samples[i].paternal, samples[i].maternal, samples[i].sex, samples[i].phenotype].join('\t'));
            samples_list.add(sample_name)
        }
    }
    return outfile
}

// Function to get a list of metadata (e.g. case id) for the case [ meta ]
def create_case_channel(List rows) {
    def case_info    = [:]
    def probands     = []
    def upd_children = []
    def father       = ""
    def mother       = ""

    for (item in rows) {
        if (item.phenotype == "2") {
            probands.add(item.sample)
        }
        if ( (item.paternal!="0") && (item.paternal!="") && (item.maternal!="0") && (item.maternal!="") ) {
            upd_children.add(item.sample)
        }
        if ( (item.paternal!="0") && (item.paternal!="") ) {
            father = item.paternal
        }
        if ( (item.maternal!="0") && (item.maternal!="") ) {
            mother = item.maternal
        }
    }

    case_info.father       = father
    case_info.mother       = mother
    case_info.probands     = probands.unique()
    case_info.upd_children = upd_children.unique()
    case_info.id           = rows[0].case_id

    return case_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
