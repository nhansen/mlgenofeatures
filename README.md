# Generating simulated and actual k-mer count-based features with mlgenofeatures

The mlgenofeatures Snakemake workflow can be used to create files of
k-mer counts from actual or simulated short read datasets from a diverse set of
mutated haplotype reference sequences.

These k-mer counts can then be used as input to scripts in the mlgenotype
python package to train random forest classifiers to recognize large structural
variants in real whole genome datasets.

## Software dependencies

The mlgenofeatures pipeline requires:

* Snakemake (https://snakemake.readthedocs.io/en/stable/)
* ART version ART-MountRainer-2016-06-05 (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
* meryl version 1.3 (https://github.com/marbl/meryl)
* Python with the pysam library (https://pysam.readthedocs.io/en/latest/)


