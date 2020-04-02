params.outputDir = "output"

// ~~~~~ START WORKFLOW ~~~~~ //
log.info "~~~~~~~ VEP Pipeline ~~~~~~~"
log.info "* Project dir:        ${workflow.projectDir}"
log.info "* Launch dir:         ${workflow.launchDir}"
log.info "* Work dir:           ${workflow.workDir.toUriString()}"
log.info "* Profile:            ${workflow.profile ?: '-'}"
log.info "* Script name:        ${workflow.scriptName ?: '-'}"
log.info "* Script ID:          ${workflow.scriptId ?: '-'}"
log.info "* Container engine:   ${workflow.containerEngine?:'-'}"
log.info "* Workflow session:   ${workflow.sessionId}"
log.info "* Nextflow run name:  ${workflow.runName}"
log.info "* Nextflow version:   ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})"
log.info "* Launch command:\n${workflow.commandLine}\n"

Channel.fromPath( file(params.ref_fa) )
.into{ ref_fa;
    ref_fa2
}

Channel.fromPath( file(params.ref_fai) )
.into{ ref_fai;
ref_fai2
}

Channel.fromPath( file(params.ref_dict) )
.into{ ref_dict;
    ref_dict2
}

Channel.fromPath( file(params.gtf) )
.into{ gtf;
    gtf2
}

Channel.fromPath("variants/**.vcf").map { item ->
    def sampleID = "${item.getName()}".replaceFirst(/.vcf$/, "")
    return([sampleID, item])
}.set { input_vcfs }

process download_ref {
    // http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
    // ftp://ftp.ensembl.org/pub/release-99/variation/indexed_vep_cache
    // GRCh37 : hg19
    // GRCh38 : hg38
    storeDir "${params.VEP_refDir}"

    output:
    file("homo_sapiens_vep_99_GRCh38") into vep_ref_dir

    script:
    """
    curl -O ftp://ftp.ensembl.org/pub/release-99/variation/indexed_vep_cache/homo_sapiens_vep_99_GRCh38.tar.gz && \
    mkdir -p homo_sapiens_vep_99_GRCh38 && \
    mv homo_sapiens_vep_99_GRCh38.tar.gz homo_sapiens_vep_99_GRCh38/ && \
    (
        cd homo_sapiens_vep_99_GRCh38 && \
        tar xzf homo_sapiens_vep_99_GRCh38.tar.gz
        )
    """
}
vep_ref_dir.map{ item ->
    def assembly = "GRCh38"
    return([item, assembly])
}.set{ vep_ref_dir_assembly }

process vep {
    // http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic
    tag "${sampleID}"
    publishDir "${params.outputDir}/VEP/raw", mode: 'copy'

    input:
    set val(sampleID), file(vcf), file(ref_dir), val(assembly), file(refFasta) from input_vcfs.combine(vep_ref_dir_assembly)
        .combine(ref_fa)
 

    output:
    set val(sampleID), file("${output_file}") into vcf_annotated
    file("${output_html}")

    script:
    prefix = "${sampleID}"
    output_file = "${prefix}.vep.vcf"
    output_html = "${vcf}".replaceFirst(/.vcf$/, ".vep.vcf_summary.html")
    """
    vep \
    --offline \
    --cache \
    --dir "${ref_dir}" \
    --assembly "${assembly}" \
    --fasta "${refFasta}" \
    --force_overwrite \
    --species homo_sapiens \
    -i "${vcf}" \
    --format vcf \
    -o "${output_file}" \
    --vcf
    """
}
