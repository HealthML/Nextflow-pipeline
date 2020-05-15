# Nextflow-pipeline
Pipeline for annotating variants in .vcf files using Variant Effect Predictor (VEP).

# Installtaion
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

The pipeline downloads automatically hg38 fasta file. However, for the current pipeline I am using reference genome (`.fa`), the annoation file (`.gtf`) and their corresponding indexed files (`.fai` & `.tbi` files). For runing the pipeline using these references, please run the following command on the terminal.

```
nextflow run main.nf --ref_fa /home/Alva.Rani/UKbiobank/derived/projects/kernels_VEP/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /home/Alva.Rani/data/reference/Homo_sapiens.GRCh38.97.gtf.gz --gtf_tbi /home/Alva.Rani/data/reference/Homo_sapiens.GRCh38.97.gtf.gz.tbi

```
If you can access the VM server and the above mentioned folder, there is index for the reference genome.

Otherwise, you can also run the whole pipeline by using the following one liner,

`./nextflow run main.nf`


