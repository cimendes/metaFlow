startReads = Channel.fromFilePairs(params.fastq)

// TRIMMOMATIC CHANNELS //
IN_trimmomatic_opts = Channel
                .value([params.trimSlidingWindow,
                        params.trimLeading,
                        params.trimTrailing,
                        params.trimMinLength])

/** INTEGRITY_COVERAGE - MAIN
This process will check the integrity, encoding and get the estimated
coverage for each FastQ pair. Corrupted FastQ files will also be detected
and filtered here.
*/
process integrity_coverage {

    tag { fastq_id }

    // This process can only use a single CPU
    cpus 1

	input:
	set fastq_id, file(fastq_pair) from startReads
	val gsize from Channel.value(0)
	val cov from Channel.value(0)
	// This channel is for the custom options of the integrity_coverage.py
	// script. See the script's documentation for more information.
	val opts from Channel.value('')

	output:
	set fastq_id, file('*_phred') into SIDE_phred
	set fastq_id, file(fastq_pair) into integrityReads

	script:
	template "integrity_coverage.py"
}

/** FASTQC - MAIN
This process will perform the fastQC analysis for each sample. In this run,
the output files (summary and data) of FastQC are sent to the output channel
as pair_1* and pair_2* files.
*/
process fastQC {

        tag { fastq_id }

        input:
        set fastq_id, file(fastq_pair) from integrityReads
        val ad from Channel.value('None')

        output:
        set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into MAIN_fastqc_out
        //set fastq_id, val("fastqc"), file(".status") into STATUS_fastqc

        script:
        template "fastqc.py"
}

/** FASTQC_REPORT - MAIN
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report {

    tag { fastq_id }

    cpus 1
    publishDir 'reports/fastqc/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out
    val opts from Channel.value("--ignore-tests")
    //substituir por string vazia se quiser correr testes

    output:
    set fastq_id, file(fastq_pair), '.status', 'optimal_trim' into MAIN_fastqc_trim
    file '*_trim_report' into LOG_trim
    file "*_status_report" into LOG_fastqc_report
    file "${fastq_id}_*_summary.txt" optional true

    script:
    template "fastqc_report.py"
}

/** TRIMMOMATIC - MAIN
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic {

    tag { fastq_id }

    input:
    //magia negra
    set fastq_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim.phase(SIDE_phred).map{ [it[0][0], it[0][1], file(it[0][3]).text, it[1][1]] }
    val opts from IN_trimmomatic_opts

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into MAIN_trimmomatic_out, SIDE_bowtie_in
    set fastq_id, val("trimmomatic"), file(".status") into STATUS_trimmomatic
    file '*_trimlog.txt' optional true into LOG_trimmomatic

    //when:
    //params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"
}

