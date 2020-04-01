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

# Reference Files

The pipeline downloads automatically hg38 fasta file. However, for the current pipeline I am using reference genome which I manually downloaded and indexed (`.fai`) and converted to `.dict`. And later provided these reference files as the command line arguments.
In addition, to that I have used the `GTF` file as one of the paramaters for the VEP tool.



