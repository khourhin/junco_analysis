#! /usr/bin/env python

import sqlite3 as lite
import logging as log
import os


def get_genotype(all_b, geno_field, GQi, GQ_thres):

    # In case no data on genotype
    if geno_field.startswith("./."):
        return "NA"
    else:
        geno_f = geno_field.split(":")
        GQ = int(geno_f[GQi])

        if GQ >= GQ_thres:
            # Assume the genotype is the first field in the genotype field
            genotype = geno_f[0].split("/")
            # Form the genotype using the 0/1 notation as index
            genotype = [all_b[int(x)] for x in genotype]
            return "/".join(genotype)
        else:
            return "FILT"


def get_SNP(line, GQ_thres):
    """
    Parse a SNP line from a vcf file
    """
    line = line.strip().split()
    chro = line[0]
    pos = line[1]
    ref = line[3]

    # Join all alternatives bases into a single string
    alt = "".join(line[4].split(","))
    all_b = ref + alt

    # Get the position of GQ in the INFO field
    GQi = line[8].split(":").index("GQ")

    geno_fields = line[9:]
    genotypes = [get_genotype(all_b, x, GQi, GQ_thres) for x in geno_fields]

    return (chro, pos, ref, all_b, "\t".join(genotypes))


def get_all_SNPs(vcf, GQ_thres, outfile):
    """
    Create a iterator over all SNPs from a vcf file
    """
    with open(vcf, "r") as f:
        for line in f:

            # Get the sample names
            if line.strip().startswith("#CHROM"):
                header = ["CHROM", "POS", "GENEID", "GENE_NAME", "REF", "ALT"]
                header += line.strip().split()[9:]
                with open(outfile, "w") as fout:
                    fout.write("\t".join(header) + "\n")
            # Skip comments
            if not line.strip().startswith("#"):
                yield get_SNP(line, GQ_thres)


def getGeneId(snp, gtf_db):
    """
    Fetch the gene name for the SNPs 
    from the gtf_db created by SNPs_annotator.py
    """
    chro = snp[0]
    pos = snp[1]

    with lite.connect(gtf_db) as con:
        con.text_factory = str  # to not have unicode formatting
        cur = con.execute("""SELECT chromo, start, end, gene_id, gene_name from GTF
        WHERE chromo=? AND start < ? AND end > ?""", (chro, pos, pos))
        data = cur.fetchone()
        if data:
            chromo, start, end, gene_id, gene_name = data
            return gene_id, gene_name
        else:
            return "NA", "NA"


def print_genotypes(vcf, GQ_thres, out_prefix, gtf_db=None):

    outfile = out_prefix + "_genotypes.csv"

    if os.path.isfile(outfile):
        raise IOError(
            "The outfile '%s' already exists and will not be overwritten" % outfile)

    log.info("CREATING genotype file: %s" % outfile)

    with open(outfile, "a") as fout:
        for snp in get_all_SNPs(vcf, GQ_thres, outfile):

            if gtf_db:
                gAnnot = getGeneId(snp, gtf_db)
                snp = snp[:2] + gAnnot + snp[2:]

            fout.write("\t".join(snp) + "\n")
