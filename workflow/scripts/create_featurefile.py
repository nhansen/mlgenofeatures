import os
import re
import gzip
from collections import namedtuple

def retrieve_genotype(match):
    genolist = []
    genolist.append(match.group(2))
    genolist.append(match.group(4))

    return genolist

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    # kmer count file has kmer counts from the simulated fastq file
    kmer_count_file = snakemake.input[0]

    # all kmer file establishes the order of the columns, so all sequence sets will generate the same columns in their feature files
    all_kmer_file = snakemake.input[1]

    # feature file will contain a column for each kmer with counts
    feature_file = snakemake.output[0]

    match = re.search(r'^features/(\S{2})\_(\S{4})\_(\S{2})\_(\S{4})\.(\d+)x\.(\d+)\.features.txt$', feature_file)

    if match:
        genotypelist = retrieve_genotype(match)
        genotypestring = "_".join(genotypelist)

        print("Found genotype " + genotypestring)
   
        kmercounts = {}
 
        print("Opening " + kmer_count_file)
        with open(kmer_count_file, "r") as inkmercount_f:
            for line in inkmercount_f:
                line = line.strip()
                countmatch = re.search(r'^(\S+)\s(\d+)$', line)
                if countmatch:
                    currentkmer = countmatch.group(1)
                    currentcount = countmatch.group(2)
                    kmercounts[currentkmer] = currentcount
                    #print("Recording kmer " + currentkmer + " with count " + currentcount)
                else:
                    print("Can\'t find kmer in line:\n" + line)

        print("Successfully read in simulated kmer counts from file " + kmer_count_file)
           
        countlist = [] 
        print("Opening " + all_kmer_file)
        with open(all_kmer_file, "r") as inallkmers_f:
            for line in inallkmers_f:
                line = line.strip()
                kmermatch = re.search(r'^(\S+)\s(\d+)$', line)
                if kmermatch:
                    printkmer = kmermatch.group(1)
                    if kmercounts.get(printkmer) == None:
                        countlist.append("0")
                        print("No count for kmer " + printkmer)
                    else:
                        countlist.append(str(kmercounts[printkmer]))
                else:
                    print("Can\'t find kmer in line:\n" + line)

        countlist.append(genotypestring)

        with open(feature_file, "w") as outfeature_f:
            featurestring = "\t".join(countlist)
            print(featurestring, file=outfeature_f)
    else:
        print("Couldn\'t parse filename " + feature_file)
    
