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
./nextflow run main.nf --ref_fa /home/Aliki.Zavaropoulou/vep_data/input/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --gtf /home/Aliki.Zavaropoulou/vep_data/input/Homo_sapiens.GRCh38.97.gtf.gz --gtf_tbi /home/Aliki.Zavaropoulou/vep_data/input/Homo_sapiens.GRCh38.97.gtf.gz.tbi

```
If you can access the VM server and the above mentioned folder, there is index for the reference genome.

Otherwise, you can also run the whole pipeline by using the following one liner,

`./nextflow run main.nf`


# About the pipeline steps

The flow of the tasks is encoded in processes, which are executed in the same oder as they appear in main.nf top to down.
Firstly, process download_ref downloads automatically the human reference genome hg38 fasta file from the webpage Ensembl.
Secondly, process pling_1 uses the tool PLINK and takes as input SNP data given in .bed, .bim, .fam format and processes those according to quality control criteria according to Hardy-Weinberg equilibrium.
Thirdly, process pling_2 uses the tool PLINK and takes as input the processed SNP data and excludes samples from the pipeline given a .txt file. The output is a zipped vcf file.
In the fourth place, process vep uses the tool VEP and takes as input the vcf file combining files for the annotation and indexing of the reference genome. Output is the functional annotation of variants.
In the fifth place, process seak_analysis performs genotype-phenotype association test using kernel functions applying the package seak available in python. Here, the output of the process vcf has to be transformed into a specific format. Afterwards comes the association test on a user-selected phenotype and the beforehand produced functional annotation.
