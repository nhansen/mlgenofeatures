import os
import re
import gzip
import pysam
from collections import namedtuple

def print60(seq, outfh):
    for i in range(0, len(seq), 60):
        if i+60 <= len(seq):
            print(seq[i:i+60], file=outfh)
        else:
            print(seq[i:len(seq)], file=outfh)

def gatherdata(match):
    Samplegenodata = namedtuple("Samplegenodata", ["samples", "genotypes", "genomefiles", "mutationfiles"])

    samplelist = []
    genolist = []
    genofilelist = []
    mutfilelist = []
    samplelist.append(match.group(1))
    samplelist.append(match.group(3))
    genolist.append(match.group(2))
    genolist.append(match.group(4))
    genofilelist.append(snakemake.config["genomefiles"][samplelist[0]])
    genofilelist.append(snakemake.config["genomefiles"][samplelist[1]])
    mutfilelist.append(snakemake.config["mutationfiles"][samplelist[0]])
    mutfilelist.append(snakemake.config["mutationfiles"][samplelist[1]])

    return Samplegenodata(samplelist, genolist, genofilelist, mutfilelist)

def retrievemutdata(mutfile, geno):
    Mutationdata = namedtuple("Mutationdata", ["mutation", "mutentry", "mutstart", "mutend", "regionoffset", "regionstart", "regionend", "regionstrand"])

    # Determine whether a mutation is needed for this genotype:
    mutation = "none"
    mutentry = ""
    mutstart = 0
    mutend = 0
    regionoffset = 0
    regionstrand = ''
    inmuts = open(mutfile, "r")
    for line in inmuts:
        fields = line.split()
        if fields[3] == geno:
            mutation = fields[5]
            mutentry = fields[0] 
            mutstart = int(fields[1])
            mutend = int(fields[2])
        if fields[3] == "REGION":
            mutentry = fields[0] 
            regionoffset = int(fields[1])
            regionstart = regionoffset
            regionend = int(fields[2])
            regionstrand = fields[4]
    inmuts.close()

    # mutstart is zero-based start within AlphaThal region, regionoffset is zero-based
    # start of AlphaThal region in contig, so sum of mutstart and regionoffset is the
    # zero-based start of the deletion, in coordinates along the contig. First deleted
    # base is mutstart + 1 (one-based) or mutstart (zero-based)
    # mutend is one based, so after summing, mutend is one-based end of deletion in 
    # coordinates within the contig

    mutstart = mutstart + regionoffset
    mutend = mutend + regionoffset

    # if there is a deletion or duplication, the regionend value will be altered:
    if mutation == "deletion":
        regionend = regionend - (mutend - mutstart)
    if mutation == "duplication":
        regionend = regionend + (mutend - mutstart)

    print("Found mutation " + mutation + " with position " +  str(mutstart) + "-" + str(mutend) + " in region from " + str(regionstart) + "-" + str(regionend) )

    return Mutationdata(mutation, mutentry, mutstart, mutend, regionoffset, regionstart, regionend, regionstrand)

def truncateseq(inputseq, fiveprimeflank, threeprimeflank, capturestart, captureend, strand):
    print("Applying 5p buffer " + str(fiveprimeflank))
    print("Applying 3p buffer " + str(threeprimeflank))
    if strand == "+":
        startbuf = fiveprimeflank
        endbuf = threeprimeflank
    if strand == "-":
        startbuf = threeprimeflank
        endbuf = fiveprimeflank
    seqstart = capturestart - startbuf
    seqend = captureend + endbuf
    seqlength = len(inputseq)
    if seqstart < 0:
        seqstart = 0
        print("Only able to include " + str(capturestart) + " of five prime flank to the region of interest")
    if seqend > seqlength:
        seqend = seqlength
        print("Only able to include " + str(seqend - captureend) + " of three prime flank to the region of interest")
    print("Attempting to extract from  " + str(seqstart) + " to " + str(seqend))
    truncatedseq = inputseq[seqstart:seqend]

    return truncatedseq

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f    
    genome_file = snakemake.output[0]
    match = re.search(r'^genomes/(\S{2})\_(\S{4})\_(\S{2})\_(\S{4})\.fasta$', genome_file)
    if match:
    
        sampledata = gatherdata(match)
        samples = sampledata.samples
        genotypes = sampledata.genotypes
        genomefiles = sampledata.genomefiles
        mutationfiles = sampledata.mutationfiles
        fivepgenomeregionbuffer = int(snakemake.config["fivepgenomeregionbuffer"])
        threepgenomeregionbuffer = int(snakemake.config["threepgenomeregionbuffer"])
    
        outgenome = open(genome_file, "w")
        for genomeindex in range(2):
            inputgenomefile = snakemake.config["genomedir"] + "/" + genomefiles[genomeindex]
            mutationfile = snakemake.config["mutationdir"] + "/" + mutationfiles[genomeindex]
            thisgeno = genotypes[genomeindex]
    
            mutdata = retrievemutdata(mutationfile, thisgeno)
            mutation = mutdata.mutation
            mutentry = mutdata.mutentry
            mutstart = mutdata.mutstart
            mutend = mutdata.mutend
            regionstrand = mutdata.regionstrand
            regionoffset = mutdata.regionoffset
            regionstart = mutdata.regionstart
            regionend = mutdata.regionend
    
            currentseqid = ''
            currentseq = ''
            startbuf = 0
            endbuf = 0

            fastagenome = pysam.Fastafile(inputgenomefile)
            contigseq = fastagenome.fetch(mutentry)
            print("Length of entry " + mutentry + ":" + str(len(contigseq)))
            if mutation == "deletion":
                currentseq = contigseq[0:mutstart] + contigseq[mutend:len(contigseq)]
            elif mutation == "duplication":
                currentseq = contigseq[0:mutend] + contigseq[mutstart:mutend] + contigseq[mutend:len(contigseq)]
            else:
                currentseq = contigseq[0:len(contigseq)]
            print("Length of mutated entry " + mutentry + ":" + str(len(currentseq)))
            if snakemake.config["fivepgenomeregionbuffer"] != "None":
                currentseq = truncateseq(currentseq, fivepgenomeregionbuffer, threepgenomeregionbuffer, regionstart, regionend, regionstrand)
            print(">" + mutentry, file=outgenome)
            print60(currentseq, outgenome)
    
            #ingenome = gzip.open(inputgenomefile, "r")
            #for line in ingenome:
                #line = str(line,'utf-8')
                #match = re.search(r'^>\s*(\S+)', line)
                #if match:
                    ## print current entry if there is one
                    #if len(currentseq) > 0:
                        ## handle case where entry contains our target region:
                        #if currentseqid == mutentry:
                            #print("Found mutation entry " + mutentry + " with length " + str(len(currentseq)))
                            #if mutation == "deletion":
                                #oneseq = currentseq.replace("\n", "")
                                #print("Length of entry " + mutentry + ":" + str(len(oneseq)))
                                #currentseq = oneseq[0:mutstart] + oneseq[mutend:len(oneseq)]
                                #print("Length of mutated entry " + mutentry + ":" + str(len(currentseq)))
                            #if mutation == "duplication":
                                #oneseq = currentseq.replace("\n", "")
                                #print("Length of entry " + mutentry + ":" + str(len(oneseq)))
                                #currentseq = oneseq[0:mutend] + oneseq[mutstart:mutend] + oneseq[mutend:len(oneseq)]
                                #print("Length of mutated entry " + mutentry + ":" + str(len(currentseq)))
                            #if snakemake.config["fivepgenomeregionbuffer"] != "None":
                                #currentseq = truncateseq(currentseq, fivepgenomeregionbuffer, threepgenomeregionbuffer, regionstart, regionend, regionstrand)
    #
                            #print(">" + currentseqid, file=outgenome)
                            #print60(currentseq, outgenome)
                        #else:
                            #if snakemake.config["fivepgenomeregionbuffer"] == "None":
                                #print(">" + currentseqid, file=outgenome)
                                #print(currentseq, file=outgenome, end="")
                            #
                    #currentseqid = match.group(1)
                    ##print("Found", currentseqid)
                    #currentseq = ''
                #else:
                    #currentseq = currentseq + line.strip()
            #ingenome.close()
            #if len(currentseq) > 0:
                #if currentseqid == mutentry:
                    #if mutation == "deletion":
                        ## delete appropriate sequence
                        #oneseq = currentseq.replace("\n", "")
                        #print("Length of entry" + mutentry + ":" + str(len(oneseq)))
                        #print("0 to " + str(mutstart) + ", " + str(len(oneseq)))
                        #currentseq = oneseq[0:mutstart] + oneseq[mutend:len(oneseq)]
                    #if mutation == "duplication":
                        #oneseq = currentseq.replace("\n", "")
                        #print("Length of entry" + mutentry + ":" + str(len(oneseq)))
                        #currentseq = oneseq[0:mutend] + oneseq[mutstart:mutend] + oneseq[mutend:len(oneseq)]
                        #print("Length of mutated entry" + mutentry + ":" + str(len(currentseq)))
                    #if snakemake.config["fivepgenomeregionbuffer"] != "None":
                        #currentseq = truncateseq(currentseq, fivepgenomeregionbuffer, threepgenomeregionbuffer, regionstart, regionend, regionstrand)
                    #print(">" + currentseqid, file=outgenome)
                    #print60(currentseq, outgenome)
                #else:
                    #if snakemake.config["fivepgenomeregionbuffer"] == "None":
                        #print(">" + currentseqid, file=outgenome)
                        #print(currentseq, file=outgenome, end="")
                        #
        outgenome.close()
    else:
        print("Couldn\'t parse filename " + genome_file)
    
