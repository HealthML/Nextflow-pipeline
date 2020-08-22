import java.nio.file.Paths

// ~~~~~ START WORKFLOW ~~~~~ //

println """\
  UKbiobank - N F   P I P E L I N E
===================================
transcriptome: ${params.gtf}
reference    : ${params.ref_fa}
outdir       : ${params.outputDir}
Project dir  : ${workflow.projectDir}
Work dir     : ${workflow.workDir.toUriString()}
Container engine: ${workflow.containerEngine?:'-'}
Nextflow version: ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
Launch dir      : ${workflow.launchDir}
Nextflow run name: ${workflow.runName}
Sipt name      : ${workflow.scriptName ?: '-'}
Launch command: \n${workflow.commandLine}\n
      """
.stripIndent()
//docker.runOptions = '-u $(id -u):$(id -g)'

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
<<<<<<< HEAD
    .fromFilePairs("${params.plink_input}/ukb_SPB_50k_exome_seq.{bed,bim,fam}",size:3) {
        file -> file.baseName
    }
    .filter { key, files -> key in params.pops }
    .set { plink_data }


params.pops = "ukb_SPB_50k_exome_seq"
dir = params.plink_input
sampleID= "ukb_SPB_50k_exome_seq"

Channel
    .fromFilePairs("${params.plink_input}/ukb_SPB_50k_exome_seq.{bed,bim,fam}",size:3) {
=======
    .fromFilePairs("${params.dir}/${params.pops}.{bed,bim,fam}",size:3) {
>>>>>>> dev_nf
        file -> file.baseName
    }
    .filter { key, files -> key in params.pops }
    .set { plink_data }

Channel.fromPath("${params.dir}/**.bed").map{item ->
       def sampleID = "${item.getName()}".replaceFirst(/.bed$/, "")
       return(sampleID)}.set { sampleid }

<<<<<<< HEAD
params.vcf = "/home/Alva.Rani/UKbiobank/derived/projects/kernels_VEP/vep_SPB_out/vep_ensembl/v2/next"
Channel.fromPath("${params.outputDir}/ukb_SPB_50k_exome_seq_filtered_vcf/**.vcf").map { item ->
    def sampleID = "${item.getName()}".replaceFirst(/.vcf$/, "")
    return([sampleID, item])
}.set { input_vcfs }

//input_vcfs.subscribe { println () }
=======
//sampleid.subscribe { println it }
>>>>>>> dev_nf

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
<<<<<<< HEAD
    publishDir "${params.outputDir}/ukb_SPB_50k_exome_seq_filtered"

    input:
    set pop, file(pl_files) from plink_data

    output:
    file "ukb_SPB_50k_exome_seq_filtered.{bed,fam,bim}" into pling1_results
=======
  //  publishDir "${params.outputDir}/ukb_SPB_50k_exome_seq_filtered"
    publishDir "${params.outputDir}/ukb_SPB_50k_exome_seq_filtered", mode: 'link'

    input:
    set pop, file(pl_files) from plink_data
>>>>>>> dev_nf

    output:
   // file "${pop}_filtered.{bed,fam,bim}" into pling1_results
   set file("${pop}_filtered.bed"), file("${pop}_filtered.bim"), file("${pop}_filtered.fam") into pling1_results
    script:
<<<<<<< HEAD
    output_file="ukb_SPB_50k_exome_seq_filtered"
=======
    output_file = "${pop}_filtered"
    base        = pl_files[0].baseName
>>>>>>> dev_nf

        """
        plink2 \
        --bfile $pop \
        --hwe 0.00001 \
        --make-bed \
        --out ${output_file} \
     """
}

pling1_results
   .collect()
   .flatten()
    .map { file -> tuple(file.baseName, file)}
    .groupTuple(by: 0)
    .map { input -> tuple(input[0], input[1][0], input[1][1], input[1][2])}
    .set { pl1 }

process pling_2 {
    // publishDir "${params.outputDir}/ukb_SPB_50k_exome_seq_filtered_vcf"
    publishDir "${params.outputDir}/ukb_SPB_50k_exome_seq_filtered_vcf", mode:'copy'

    input:
    //set file(bed), file(bim), file(fam) from pling1_results.view()
     set val(pop1),file(bed), file(bim), file(fam),file(fam1) from pl1.combine(fam_for_plink2)
     //file(fam1) from fam_for_plink2

    output:
<<<<<<< HEAD
    file "ukb_SPB_50k_exome_seq_filtered.vcf.gz" into pling2_results

    script:
    output_file="ukb_SPB_50k_exome_seq_filtered"
    base       = pling1[0].baseName

     """
     plink2 \
     --bfile $base \
=======
    file("${pop1}.vcf.gz") into pling2_results
   //file ("${output_file}") into pling2_results

    script:
    //prefix        = ['ukb_SPB_50k_exome_seq','ukb_FE_50k_exome_seq']
    //output_file   = "${prefix}.vcf.gz"
    base            = bed.baseName
    output_file     = "${pop1}"

     """
     plink2 \
     --bfile $pop1  \
>>>>>>> dev_nf
     --keep-fam ${params.fam}/ukb_50k_exome_seq_filtered_for_VEP_ID.txt \
     --recode vcf-iid bgz --out ${output_file}
     """
}
process vep {
<<<<<<< HEAD
    // http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basie
    tag "${sampleID}"
    publishDir "${params.outputDir}/VEP"

    input:
        set file(ref_dir), val(assembly), file(refFasta), file(GTF), file(GTF_tbi) from vep_ref_dir_assembly.combine(ref_fa)
        .combine(gtf)
        .combine(gtf_tbi)
         file(vcf) from pling2_results

    output:
     file("${output_file}") into vcf_annotated
     file("${output_html}")

    script:
  //  prefix = "${sampleID}"
    output_file = "ukb_SPB.vep.vcf"
    output_html = "ukb_SPB.vep.vcf_summary.html"

    """
    vep \
    --fasta "${refFasta}" \
    --format vcf --force_overwrite \
    --dir "${ref_dir}" \
    --assembly "${assembly}" \
    --gtf "${GTF}" \
    --force_overwrite \
    --species homo_sapiens \
    --input_file ${vcf} \
    --stats_file "${output_html}" \
    --output_file "${output_file}"
=======
          // http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basie
          tag "Runing VEP on $vcf"
          // publishDir "${params.outputDir}/VEP"
          publishDir "${params.outputDir}/VEP",mode:'copy'

          input:
              set file(ref_dir), val(assembly),file(vcf) from vep_ref_dir_assembly.combine(pling2_results)

          output:
          file("${output_file}") into vcf_annotated
          file("${output_html}")

script:
          // prefix    = "${sampleID}"
          //output_file = "ukb_SPB.vep.vcf"
          //output_html = "ukb_SPB.vep.vcf_summary.html"
  base          = vcf.baseName
  output_file   = "${base}.vep.vcf"
  output_html   = "${base}_summary.html"

  """
  vep \
  -- cache \
  --format vcf \
  --offline \
  --dir "${ref_dir}" \
  --assembly "${assembly}" \
  --force_overwrite \
  --species homo_sapiens \
  --input_file ${base}.gz \
  --sift b \
  --polyphen b \
  --stats_file "${output_html}" \
  --output_file "${output_file}"
  """
    }


//pling2_results.println { "plink2: $it" }

process seak_analysis_1 {
    tag "Running seak on $intermediate_vcf"
    publishDir "${params.outputDir}/seaktsv", mode:'link'
    // conda '/home/Aliki.Zavaropoulou/miniconda3/envs/py3'

    input:
    file(vcf_vep) from vcf_annotated.view()
    // file(pling2_vcf) from pling2_results

    output:
    file("${filtered_VEP}")  into vcf_filtered
    //file("${intermediate_vcf_o}")  into vcf_filtered

    script:

    intermediate_vcf     =  vcf_vep[0].baseName
    intermediate_vcf_o   = "${intermediate_vcf}_filteredvar.vep.vcf"
    filtered_VEP         = "${intermediate_vcf}_LOF_filtered.tsv"

    """
    awk '\$NF ~ /IMPACT=HIGH/' ${intermediate_vcf}.vcf > ${intermediate_vcf_o}
    python /filter_VEP.py \
      -i ${intermediate_vcf}_filteredvar.vep.vcf \
      -o ${filtered_VEP}
>>>>>>> dev_nf
    """
}


//  python /run_test_proteinlof.py \
//  -pheno 'ApoA' \
//  -g ${intermediate_vcf} \
//  -ukbdir ${params.outputDir} \
//  -i ${filtered_VEP}
