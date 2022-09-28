import os
import re
import gzip

def print60(seq, outfh):
    for i in range(0, len(seq), 60):
        if i+60 <= len(seq):
            print(seq[i:i+60], file=outfh)
        else:
            print(seq[i:len(seq)], file=outfh)
    
genome_file = snakemake.output[0]
match = re.search(r'^genomes/(\S{2})\_(\S{4})\_(\S{2})\_(\S{4})\.fasta$', genome_file)
if match:
    samples = []
    genotypes = []
    genomefiles = []
    mutationfiles = []
    samples.append(match.group(1))
    samples.append(match.group(3))
    genotypes.append(match.group(2))
    genotypes.append(match.group(4))
    genomefiles.append(snakemake.config["genomefiles"][samples[0]])
    genomefiles.append(snakemake.config["genomefiles"][samples[1]])
    mutationfiles.append(snakemake.config["mutationfiles"][samples[0]])
    mutationfiles.append(snakemake.config["mutationfiles"][samples[1]])

    outgenome = open(genome_file, "w")
    for genomeindex in range(2):
        inputgenomefile = snakemake.config["genomedir"] + "/" + genomefiles[genomeindex]
        mutationfile = snakemake.config["mutationdir"] + "/" + mutationfiles[genomeindex]
        thisgeno = genotypes[genomeindex]

        # Determine whether a mutation is needed for this genotype:
        mutation = "none"
        mutentry = ""
        mutstart = 0
        mutend = 0
        regionoffset = 0
        inmuts = open(mutationfile, "r")
        for line in inmuts:
            fields = line.split()
            if fields[3] == thisgeno:
                mutation = fields[5]
                mutentry = fields[0] 
                mutstart = int(fields[1])
                mutend = int(fields[2])
            if fields[3] == "REGION":
                regionoffset = regionoffset + int(fields[1])
        inmuts.close()

        # mutstart is zero-based start within AlphaThal region, regionoffset is zero-based
        # start of AlphaThal region in contig, so sum of mutstart and regionoffset is the
        # zero-based start of the deletion, in coordinates along the contig. First deleted
        # base is mutstart + 1 (one-based) or mutstart (zero-based)
        # mutend is one based, so after summing, mutend is one-based end of deletion in 
        # coordinates within the contig

        mutstart = mutstart + regionoffset
        mutend = mutend + regionoffset

        print("Found mutation " + mutation + "with start" +  str(mutstart) + "-" + str(mutend) )

        currentseqid = ''
        currentseq = ''
        ingenome = gzip.open(inputgenomefile, "r")
        for line in ingenome:
            line = str(line,'utf-8')
            match = re.search(r'^>\s*(\S+)', line)
            if match:
                #print current entry if there is one
                if len(currentseq) > 0:
                    if currentseqid == mutentry:
                        if mutation == "deletion":
                            oneseq = currentseq.replace("\n", "")
                            print("Length of entry" + mutentry + ":" + str(len(oneseq)))
                            currentseq = oneseq[0:mutstart] + oneseq[mutend:len(oneseq)]
                            print("Length of mutated entry" + mutentry + ":" + str(len(currentseq)))
                        if mutation == "duplication":
                            oneseq = currentseq.replace("\n", "")
                            print("Length of entry" + mutentry + ":" + str(len(oneseq)))
                            currentseq = oneseq[0:mutend] + oneseq[mutstart:mutend] + oneseq[mutend:len(oneseq)]
                            print("Length of mutated entry" + mutentry + ":" + str(len(currentseq)))
                        print(">" + currentseqid, file=outgenome)
                        print60(currentseq, outgenome)
                    else:
                        print(">" + currentseqid, file=outgenome)
                        print(currentseq, file=outgenome, end="")
                        
                currentseqid = match.group(1)
                print("Found", currentseqid)
                currentseq = ''
            else:
                line.strip()
                currentseq = currentseq + line
        ingenome.close()
        if len(currentseq) > 0:
            if currentseqid == mutentry:
                if mutation == "deletion":
                    # delete appropriate sequence
                    oneseq = currentseq.replace("\n", "")
                    print("Length of entry" + mutentry + ":" + str(len(oneseq)))
                    print("0 to " + str(mutstart) + ", " + str(len(oneseq)))
                    currentseq = oneseq[0:mutstart] + oneseq[mutend:len(oneseq)]
                print(">" + currentseqid, file=outgenome)
                print60(currentseq, outgenome)
            else:
                print(">" + currentseqid, file=outgenome)
                print(currentseq, file=outgenome, end="")
                    
    outgenome.close()
else:
    print("Couldn\'t parse filename " + genome_file)

