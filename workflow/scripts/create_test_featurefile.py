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

    # calculate total counts for this sample:
    # append normalized kmer counts to a list for printing:

    kmertotal = 0
    numkmers = 0   
    countlist = [] 
    with open(all_kmer_file, "r") as inallkmers_f:
        for line in inallkmers_f:
            line = line.strip()
            kmermatch = re.search(r'^(\S+)\s(\d+)$', line)
            if kmermatch:
                numkmers = numkmers + 1
                printkmer = kmermatch.group(1)
                if kmercounts.get(printkmer) == None:
                    countlist.append("0")
                else:
                    countlist.append(str(kmercounts[printkmer]))
                    kmertotal = kmertotal + int(kmercounts[printkmer])
            else:
                print("Can\'t find kmer in line:\n" + line)

        print("Normalizing")
        normcountlist = []
        for count in countlist:
            normcountlist.append(str(round(50.0*int(count)*numkmers/kmertotal)))
        normcountlist.append("Unknown")

    with open(feature_file, "w") as outfeature_f:
        featurestring = "\t".join(normcountlist)
        print(featurestring, file=outfeature_f)
