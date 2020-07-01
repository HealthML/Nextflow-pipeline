params.outputDir = "/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/test_pipeline"
import java.nio.file.Paths
// export PATH="/home/Aliki.Zavaropoulou/miniconda3/bin:$PATH" , run first in the bash shell or add it in the .bashrc so it can source activate the environment
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

Channel.fromPath( file(params.gtf_tbi) )
.into{ gtf_tbi;
gtf_tbi2
}

Channel.fromPath( file(params.ref_dict) )
.into{ ref_dict;
    ref_dict2
}

Channel.fromPath( file(params.gtf) )
.into{ gtf;
    gtf2
}

Channel
    .fromPath("${params.fam}/ukb_50k_exome_seq_filtered_for_VEP_ID.txt", checkIfExists:true )
    .ifEmpty { exit 1, "Fam file NOT found: ${params.fam}" }
    //.println()
    .set { fam_for_plink2 }

//fam_for_plink2.subscribe { println it }

Channel
    .fromFilePairs("${params.plink_input}/ukb_SPB_50k_exome_seq.{bed,bim,fam}",size:3) {
        file -> file.baseName
    }
    .filter { key, files -> key in params.pops }
    .set { plink_data }


params.pops = "ukb_SPB_50k_exome_seq"
dir = params.plink_input

Channel
    .fromFilePairs("${params.plink_input}/ukb_SPB_50k_exome_seq.{bed,bim,fam}",size:3) {
        file -> file.baseName
    }
    .filter { key, files -> key in params.pops }
    .set { plink_data }

// plink_data.subscribe { println "$it" }

Channel.fromPath("${params.outputDir}/ukb_SPB_50k_exome_seq_filtered_vcf/**.vcf").map { item ->
    def sampleID = "${item.getName()}".replaceFirst(/.vcf$/, "")
    return([sampleID, item])
}.set { input_vcfs }

//input_vcfs.subscribe { println () }

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

process pling_1 {
    //publishDir "${params.outputDir}/ukb_FE_50k_exome_seq_filtered"
    storeDir "${params.outputDir}/ukb_FE_50k_exome_seq_filtered"

    input:
    set pop, file(pl_files) from plink_data
    output:
    file "ukb_FE_50k_exome_seq_filtered.{bed,fam,bim}" into pling1_results

    script:
    output_file="ukb_FE_50k_exome_seq_filtered"

     """
        plink2 \
        --bfile $pop \
        --hwe 0.00001 \
        --make-bed \
        --out ${output_file} \
     """
}

process pling_2 {
    //publishDir "${params.outputDir}/ukb_SPB_50k_exome_seq_filtered_vcf"
    storeDir "${params.outputDir}/ukb_SPB_50k_exome_seq_filtered_vcf"

    input:
    file(pling1) from pling1_results
    file(fam1) from fam_for_plink2

    output:
    file "ukb_FE_50k_exome_seq_filtered.vcf" into pling2_results

    script:
    output_file="ukb_FE_50k_exome_seq_filtered"

     """
     plink2 \
     --bfile "${params.outputDir}/ukb_FE_50k_exome_seq_filtered/ukb_FE_50k_exome_seq_filtered" \
     --keep-fam ${params.fam}/ukb_50k_exome_seq_filtered_for_VEP_ID.txt \
     --recode vcf-iid bgz --out ${output_file}
     gunzip ${output_file}.vcf.gz
     """
}

process vep {

    // http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basie
    tag "${sampleID}"
    storeDir "${params.outputDir}/VEP"

    input:
    set val(sampleID),file(vcf), file(ref_dir), val(assembly), file(refFasta), file(GTF), file(GTF_tbi) from input_vcfs.combine(vep_ref_dir_assembly)
        .combine(ref_fa)
        .combine(gtf)
        .combine(gtf_tbi)

    output:
    set val(sampleID), file("${output_file}") into vcf_annotated
    file("${output_html}")

    script:
    prefix = "${sampleID}"
    output_file = "${prefix}.vep.txt"
    output_html = "${vcf}".replaceFirst(/.vcf$/, ".vep.txt_summary.html")

    """
    vep \
    --offline \
    --cache  \
    --dir "${ref_dir}" \
    --assembly "${assembly}" \
    --fasta "${refFasta}" \
    --gtf "${GTF}" \
    --force_overwrite \
    --species homo_sapiens \
    -i "${vcf}" \
    --stats_file "${output_html}" \
    --tab \
    -o stdout | \
    filter_vep -filter "Consequence is stop_lost or Consequence is start_lost or Consequence is splice_donor_variant or Consequence is frameshift_variant or Consequence is splice_acceptor_variant or Consequence is stop_gained" \
    -o "${output_file}"
    """
}


process seak_analysis {
    publishDir "${params.outputDir}/seaktsv"
    conda '/home/Aliki.Zavaropoulou/miniconda3/envs/py3'

    input:
    set val(sampleID), file(vcf_vep) from vcf_annotated
    // file(pling2_vcf) from pling2_results

    output:
    file "LOF_filtered.tsv"  into vcf_filtered

    script:
    filtered_VEP = "LOF_filtered.tsv"
    intermediate_vcf = "/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/test_pipeline/ukb_SPB_filteredvariants.vep.vcf"
    //target_phenotype = 'ApoA'   ${vcf_vep}$   /home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/vep_SPB_out/vep_ensembl/v2/output/vep/ukb_SPB.vep.vcf
    """
    awk '\$NF ~ /IMPACT=HIGH/' /home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/vep_SPB_out/vep_ensembl/v2/output/vep/ukb_SPB.vep.vcf > ${intermediate_vcf}
    python /home/Aliki.Zavaropoulou/pipeline/Nextflow-pipeline/seak/filter_VEP.py \
    -i "${intermediate_vcf}" \
    -o "${filtered_VEP}"
    python /home/Aliki.Zavaropoulou/pipeline/Nextflow-pipeline/seak/run_test_proteinlof.py \
    -pheno 'ApoA' \
    -i "${filtered_VEP}"
    """
}
// awk '\$NF ~ /IMPACT=HIGH/' $vcf_vep > ${intermediate_vcf}
// to install seak library: cd seak_call   (new line)   git clone https://github.com/HealthML/seak.git   (new line)   python setup.py install
