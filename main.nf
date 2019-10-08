#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/nanorna
========================================================================================
 nf-core/nanorna Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/nanorna
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""
    ============================================================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nanorna --directRNA "file1 file2" --ref "reference_genome. " --bed12 "file "[--ram int] [--threads int]"
    Mandatory arguments
      --ref                         The reference genome file
      --outdir                      Absolute path to directory to output data (default is within nextflow process folder)
      --bed12                       Bed12 file to use with the --bed-junc flag in order to improve Alignment
                                    n.b. You can generate the bed12 file using "k8 ./paftools.js gff2bed genomeAnnotationFile.gtf"
                                    (requires the k8 javascript shell to run)
    And at least one of the following
      --directRNA                   A fastq file containing direct RNA. A gzipped file type can also be passed as an argument
      --cDNA                        A fastq (or gzipped) file containing data from either PCR-cDNA or direct cDNA.
      --custom                      Enter custom command line flags for minimap2 alignment.
                                    No error checking is performed on command line input for custom option, so please manually test on a small dataset first.
                                    Enter command line input surrounded in single quotesfollowed by desired files, all space separated.
                                    e.g. --custom "'-ax -k26' file1 file2"
                                    This will be inputted here:
                                    "minimap2-2.17 {your flags} --junc-bed bed12file --ref reffile"
    Optional
      --ram                         An integer specifying the GB of RAM to use (default is 16)
      --threads                     An integer specifying the number of threads (default is 3)
      --aligner                     Software to use for alignment: graphmap2 or minimap2 (default is minimap2)

    Help
      --help                        Using this flag will bring up usage and then exit (no alignment will occur)
    ==========================================================================================
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

if (!params.outdir || params.outdir == true)
  exit 1, "Please specify output directory (complete path)"


if (!params.ref || params.ref == true)
  exit 1, "Please specify reference genome (complete path)"
if (params.ref.endsWith('.fa') || params.ref.endsWith('.fa.gz')) {
  refgenome = file(params.ref, checkIfExists: true)
} else {
  exit 1, "Incorrect file type for reference genome: fa or fa.gz required."
}


// if no alignment flag called
if (!params.directRNA && !params.cDNA && !params.custom){
  exit 1, "No input (directRNA, cDNA or custom)."
}

// if flag used but no files given
// True is stored in param if it is called but no variable assigned
if (params.directRNA == true || params.cDNA == true || params.custom == true){
  exit 1, "Parameter used but no files named."
}

// check aligner
if (params.aligner != "graphmap2" && params.aligner != "minimap2")
  exit 1, "Please enter a valid aligner: graphmap2 or minimap2"

if (params.aligner == "graphmap2" && params.index){
  if (params.index.endsWith(".gmidx") == false)
    exit 1, "gmidx file required for graphmap2 index - please amend or enter no index"
} else if (params.aligner == "minimap2" && params.index){
  if (params.index.endsWith(".mmi") == false)
    exit 1, "mmi file required for minimap2 index - please amend or enter no index"
}

minimap2call = ""
index = ""
if (params.directRNA && !params.cDNA && !params.custom){
  minimap2call = "-ax splice -uf -k14"
  directRNAList = params.directRNA.split(" ").collect{file(it)}
  file_ch = Channel.fromPath(directRNAList, checkIfExists: true)
} else if (params.cDNA && !params.directRNA && !params.custom){
  cDNAList = params.cDNA.split(" ").collect{file(it)}
  file_ch = Channel.fromPath(cDNAList, checkIfExists: true)
  minimap2call = "-ax splice -uf"
} else if (params.custom && !params.cDNA && !params.directRNA){
  list = params.custom.split("'")
  println(list)
  customCall = list[1]
  files = list[2]
  customList = files.split(" ")
  // remove empty elements
  customList = customList - [""]
  // transform to file paths
  customList = customList.collect{file(it)}
  // add to chennel, confirm files exist
  file_ch = Channel.fromPath(customList, checkIfExists: true)
  minimap2call = params.customCall
} else {
  exit 1, "Invalid option: --directRNA OR --cDNA OR --custom"
}

// println("chromosome sizes = $params.chrSizes, skipvis = $params.skipvis")
if (!params.chrSizes  && !params.skipvis){
  process getChrSizes {
    publishDir "${params.outdir}/chrSizes", mode: 'copy'

    input:
    file genome from refgenome

    output:
    file("${genome.simpleName}.chrSizes.txt") into chrsize_ch1, chrsize_ch2
    file("${genome.baseName}.chrSizes.ucsc.txt") into ucsc_chrsize_ch1, ucsc_chrsize_ch2

    script:
    """
    samtools faidx ${genome} > ${genome}.fai
    cut -f1,2 ${genome}.fai > ${genome.simpleName}.chrSizes.txt
    bin/formatUCSC.pl ${genome.simpleName}.chrSizes.txt > \
    ${genome.baseName}.chrSizes.ucsc.txt
    """
  }
} else if (params.chrSizes && !params.skipvis){
  chrSizes = file(params.chrSizes, checkIfExists: true)
}

juncbed = ""
if (params.juncbed == true)
    exit 1, "Please specify a path to use --juncbed"
if (params.juncbed && params.aligner == "graphmap2")
  exit 1, "--juncbed option only available with minimap2"
if (params.juncbed && params.aligner == "minimap2"){
  juncfile = file(params.juncbed, checkIfExists: true)
  juncbed = "--junc-bed $juncfile"
}

graphmapIndex = ""
if (params.index){
  index = file(params.index, checkIfExists: true)
  mmi_ch = Channel.from(index)
  graphmapIndex = " -i $index"
} else if (params.aligner == "minimap2"){
  process minimapIndex {
    input:
    file genome from refgenome

    output:
    file "*.mmi" into mmi_ch

    when:
    !params.index

    """
    minimap2 ${minimap2call} -I ${params.ram}G -t ${params.threads} -d \
    ${genome.baseName}.mmi ${genome}
    """
  }
}

if (params.aligner == "minimap2"){

  process minimapAlign {
    input:
    file rawdata from file_ch
    file mmi from mmi_ch

    output:
    file "*.sam" into sam_ch

    """
    minimap2 ${minimap2call} -I ${params.ram}G -t ${params.threads} ${juncbed} \
    $mmi ${rawdata} > ${rawdata.baseName}.sam
    """
  }
} else if (params.aligner == "graphmap2"){
  process graphmapAlign {
    echo true
    input:
    file rawdata from file_ch
    file genome from refgenome

    output:
    file "*.sam" into sam_ch

    // TODO
    """
    graphmap2 align -x rnaseq -r $rawdata -d $genome -o ${rawdata.baseName}.sam \
    -t ${params.threads} $graphmapIndex
    """
  }
}

process sam2bam {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file sam from sam_ch

  output:
  file "*.bam" into indexbam_ch, bam2bed12_ch, bam2bedgraph_ch, bam2bigwig_ucsc_ch

  script:
  """
  samtools view -Sb $sam | samtools sort -o ${sam.baseName}.bam -
  """
}

process indexbam {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file bam from indexbam_ch

  output:
  file "*.bai" into bai_ch

  script:
  """
  samtools index $bam
  """
}

if (params.skipvis != true){
  process bam2bed12 {
    input:
    file bam from bam2bed12_ch

    output:
    file "${bam.baseName}.bed12" into makeUCSCbed12_ch, bed12_ch
    file "${bam.baseName}.sorted.bed12" into sorted_bed12_ch

    script:
    """
    bamToBed -bed12 -cigar -i $bam > ${bam.baseName}.bed12
    bedtools sort -i ${bam.baseName}.bed12 > ${bam.baseName}.sorted.bed12
    """
  }

  process makeUCSCbed12 {
    input:
    file bed12 from makeUCSCbed12_ch

    output:
    file "${bed12.baseName}.ucsc.bed12" into ucsc_bed12_ch
    file "${bed12.baseName}.sorted.ucsc.bed12" into sorted_ucsc_bed12_ch

    script:
    """
    bin/formatUCSC.pl $bed12 > ${bed12.baseName}.ucsc.bed12
    bedtools sort -i ${bed12.baseName}.ucsc.bed12 > ${bed12.baseName}.sorted.ucsc.bed12
    """
  }

  // TODO: get this working with the setup-nanorna as input ~ no hardcoded chrSizes.txt
  process bam2bedgraph {
    input:
    file bam from bam2bedgraph_ch

    output:
    file "${bam.baseName}.bedgraph" into bedgraph_UCSC_ch, bedgraph_ch

    script:
    """
    genomeCoverageBed -split -bg -ibam $bam > ${bam.baseName}.bedgraph
    bedSort ${bam.baseName}.bedgraph ${bam.baseName}.bedgraph
    """
  }

  process bedgraph2bw {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    file bedgraph from bedgraph_ch
    file chrSizes from chr_size_ch1

    output:
    file "*.bw" into bw_ch

    """
    bedGraphToBigWig ${bedgraph} $chrSizes ${bedgraph.baseName}.bw
    """
  }

  // TODO: get this working with the setup-nanorna as input ~ no hardcoded chrSizes.txt
  process bedgraph2UCSCbw {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file bedgraph from bedgraph_UCSC_ch
    file ucscChrSize from ucsc_chrsize_ch1

    output:
    file "*.ucsc.bw" into ucsc_bw_ch

    script:
    """
    bin/formatUCSC.pl $bedgraph > ${bedgraph.baseName}.ucsc.bedgraph
    bedSort ${bedgraph.baseName}.ucsc.bedgraph ${bedgraph.baseName}.ucsc.bedgraph
    bedGraphToBigWig ${bedgraph.baseName}.ucsc.bedgraph $ucscChrSize ${bedgraph.baseName}.ucsc.bw
    """
  }
  //TODO
  process bed12ToBigBed {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file bed12 from sorted_bed12_ch
    file chrSizes from chr_size_ch2

    output:
    file "*.bb" into bb_ch

    script:
    """
    bedToBigBed $bed12 $chrSizes ${bed12.simpleName}.bb
    """
  }

  process UCSCbed12ToBigBed {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file bed12 from sorted_ucsc_bed12_ch
    file ucscChrSize from ucsc_chrsize_ch2

    output:
    file "*.bb" into bb_ucsc_ch

    script:
    """
    bedToBigBed $bed12 $ucscChrSize ${bed12.simpleName}.ucsc.bb
    """
  }
}

custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['Fasta Ref']        = params.fasta
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-nanorna-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/nanorna Workflow Summary'
    section_href: 'https://github.com/nf-core/nanorna'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}



/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}



/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}



/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/nanorna] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/nanorna] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/nanorna] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/nanorna] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/nanorna] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/nanorna] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/nanorna]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/nanorna]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/nanorna v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
