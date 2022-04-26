import os

genome_file = snakemake.output[0]
assembly_dir = "/cluster/ifs/projects/AlphaThal/haplotype_assemblies/"
syscommand = "gunzip -c " + assembly_dir + "HG01109.maternal.f1_assembly_v2.fa.gz " + assembly_dir + "HG01109.paternal.f1_assembly_v2.fa.gz > " + genome_file
print(syscommand)
os.system(syscommand)

