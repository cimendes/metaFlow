startReads = Channel.fromFilePairs(params.fastq)

// TRIMMOMATIC CHANNELS //
IN_trimmomatic_opts = Channel
                .value([params.trimSlidingWindow,
                        params.trimLeading,
                        params.trimTrailing,
                        params.trimMinLength])

// SPADES CHANNELS //
IN_spades_opts = Channel
                .value([params.spadesMinCoverage,
                        params.spadesMinKmerCoverage])
IN_spades_kmers = Channel
                .value(params.spadesKmers)

/** 1. INTEGRITY_COVERAGE - MAIN
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
	val gsize from Channel.value(1)
	val cov from Channel.value(1)
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

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out
    val opts from Channel.value("--ignore-tests")
    //substituir por string vazia se quiser correr testes

    output:
    set fastq_id, file(fastq_pair), '.status', 'optimal_trim' into MAIN_fastqc_trim
    //file '*_trim_report' into LOG_trim
    //file "*_status_report" into LOG_fastqc_report
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
    set fastq_id, "${fastq_id}_*P*" into MAIN_trimmomatic_out
    //set fastq_id, val("trimmomatic"), file(".status") into STATUS_trimmomatic
    //file '*_trimlog.txt' optional true into LOG_trimmomatic

    script:
    template "trimmomatic.py"
}

/** BOWTIE2 - MAIN
This process will execute Bowtie2 with the hg19 indexed genome.
The unmapped reads will be saved and their headers renamed.
*/
process bowtie {

    tag { fastq_id }


    input:
    set fastq_id, file(fastq_pair) from MAIN_trimmomatic_out

    output:
    set fastq_id , "${fastq_id}*.headersRenamed_*.fq.gz" into UNMAPPED_out

    script:
    """
    bowtie2 -x /index_hg19/hg19 -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p 3 > ${fastq_id}.bam

    samtools view -buh -f 12 -o ${fastq_id}_samtools.bam -@ 2 ${fastq_id}.bam

    samtools fastq -1 ${fastq_id}_unmapped_1.fq -2 ${fastq_id}_unmapped_2.fq ${fastq_id}_samtools.bam

    python /NGStools/renamePE_samtools/renamePE_samtoolsFASTQ.py -1 ${fastq_id}_unmapped_1.fq -2 ${fastq_id}_unmapped_2.fq

    gzip *.headersRenamed_*.fq
    """
}

UNMMAPPED_forAssembly = Channel.create()
UNMMAPPED_forBowtie = Channel.create()
UNMAPPED_forCARDrgi = Channel.create()
UNMAPPED_out.into{ UNMMAPPED_forAssembly ; UNMMAPPED_forBowtie ; UNMAPPED_forCARDrgi}

/** METASPADES - MAIN
This process will execute metaSPAdes
*/
process metaspades {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from UNMMAPPED_forAssembly
    val opts from IN_spades_opts
    val kmers from IN_spades_kmers

    output:
    set fastq_id, "${fastq_id}_contigs.fasta" into MAIN_spades_out
    //set fastq_id, val("spades"), file(".status") into STATUS_spades

    script:
    //template "spades.py"
    """

    metaspades.py -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o .

    mv contigs.fasta ${fastq_id}_contigs.fasta

    sed -i 's/>/>${fastq_id}_/g' ${fastq_id}_contigs.fasta
    """
}

MAIN_spades_out_mapping = Channel.create()
MAIN_spades_out_card_rgi = Channel.create()
MAIN_spades_out.into{MAIN_spades_out_mapping ; MAIN_spades_out_card_rgi}

/** BOWTIE ASSEMBLY - MAIN
This process will execute Bowtie on the assembly
with the filtered read data
*/
process bowtie_assembly {

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from MAIN_spades_out_mapping
    set fastq_id, file(fastq_pair) from UNMMAPPED_forBowtie

    output:
    set fastq_id, file('*.bam') into MAIN_bowtie_assembly

    script:
    """
    bowtie2-build ${assembly} ${fastq_id}

    bowtie2 -x ${fastq_id} -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p 3 > ${fastq_id}.bam
    """
}

/** CARD RGI IMPLEMENTATION - MAIN
This process will execute CARD's RGI on the
assembly.
*/
process card_rgi_assembly{

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from MAIN_spades_out_card_rgi

    output

    set fastq_id, ""${fastq_id}_card_rgi.txt" into RGI_assembly

    script:
    """
    rgi main --input_sequence ${assembly} --output_file ${fastq_id}_card_rgi --input_type contig --alignment_tool DIAMOND --low_quality --include_loose -d wgs --clean
    """
}

/** CARD RGI IMPLEMENTATION - MAIN
This process will execute CARD's RGI on the
reads.
*/
process card_rgi_reads{

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from UNMAPPED_forCARDrgi

    script
}