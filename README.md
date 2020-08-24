# Nextflow-pipeline

The current pipeline analysis the exon-seq data generated from either of Regeneron’s own pipeline `(SPB)` or Functionally
Equivalent `(FE)` piplines from UKbiobank. Initially, the raw input data (`.bed`,`.fam`, `.bai`) are filtered and converted to vcf files using two `plink2` processes. Then, the variants (`vcf` file) is annotated using the ensembleVariant Effect Predictor (`vep`) tool. The corresponding result is then processed using an in-house tool called SEAK. The pipeline has a total of four processes. The tools used on all processes are containerized in the [docker image ](https://github.com/HealthML/Nextflow-pipeline/blob/master/container/Dockerfile)


# Installation
```
git clone https://github.com/HealthML/Nextflow-pipeline.git
cd Nextflow-pipeline
git checkout dev_nf
```

# Nextflow
`make install`

# VEP: Docker
To install VEP using Docker, run the Makefile command in the container directory.

```
cd container
make docker-build
```

# How to run the pipeline


In order to run the pipeline for the data generated from Regeneron’s own pipeline `(SPB)` or Functionally
Equivalent `(FE)` pipleine from UKbiobank using the VEP's cache references, please use the following command. For example, if you wanna run the samples from FE pipeline try the follwoing command on the terminal.

`./nextflow run main.nf --samples FE `

The pipeline downloads automatically hg38 fasta file. However, for the current pipeline I am using reference genome (`.fa`), the annoation file (`.gtf`) and their corresponding indexed files (`.fai` & `.tbi` files). For runing the pipeline using these references, please run the following command on the terminal.

```
./nextflow run main.nf --ref_fa /home/Alva.Rani/UKbiobank/derived/projects/kernels_VEP/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /home/Alva.Rani/data/reference/Homo_sapiens.GRCh38.97.gtf.gz --gtf_tbi /home/Alva.Rani/data/reference/Homo_sapiens.GRCh38.97.gtf.gz.tbi
```


If you can access the VM server and the above mentioned folder, there is index for the reference genome.

Otherwise, you can also run the whole pipeline by using the following one liner,

`./nextflow run main.nf`


