import os
import re

genome_file = snakemake.output[0]
match = re.search(r'^genomes/(\S{2})\_(\S{4})\_(\S{2})\_(\S{4})\.fasta$', genome_file)
if match:
    genome1 = snakemake.config["genomefiles"][match.group(1)]
    genome2 = snakemake.config["genomefiles"][match.group(3)]
    print('found', match.group(1), match.group(2), match.group(3), match.group(4))
else:
    print("Couldn\'t parse filename " + genome_file)

assembly_dir = snakemake.config["genomedir"] + "/"
syscommand = "gunzip -c " + assembly_dir + genome1 + " " + assembly_dir + genome2 + " > " + genome_file
print(syscommand)
os.system(syscommand)

