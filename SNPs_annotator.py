#! /usr/bin/env python

import sqlite3 as lite
import logging as log
import work_on_vcf as wvcf
import os
import click


def parseGtfAttributes(attrs_str):

    # Parsing the attributes of the gtf
    attrs_str = attrs_str.replace('"', '').strip()
    # Have to remove trailing ";", or final extra empty list in attrs
    attrs = attrs_str.strip(";").split(";")
    attrs = [i.strip().split(" ") for i in attrs]
    attrs_d = {k: v for k, v in attrs}

    if "gene_name" not in attrs_d:
        attrs_d["gene_name"] = "NO_NAME"

    return attrs_d


def create_gtf_db(gtf):
    """
    Create a database from a gtf file
    Return the name of the db
    """

    # outdb=os.path.basename(gtf)
    outdb = os.path.splitext(gtf)[0]
    outdb = outdb + ".db"

    if os.path.isfile(outdb):
        log.warning(
            "The database file %s already exists. Using it for further analysis." % outdb)
        return outdb

    log.info("Creating db: " + outdb)
    con = lite.connect(outdb)
    cur = con.cursor()

    cur.execute("""CREATE TABLE GTF(chromo TEXT, source TEXT, feature TEXT,
                start INT, end INT, score TEXT, strand TEXT,
                frame TEXT, gene_id TEXT, gene_name TEXT)""")

    # Skip the header (starts with #) and have only "gene entries
    with open(gtf, "r") as gtf:

        avail_att_s = set()
        for line in gtf:
            if line.startswith("#"):
                continue

            entry = line.split('\t')
            if entry[2] == "gene":

                # Get the parsed attributes
                attrs_d = parseGtfAttributes(entry[8])

                # To have a list of any attributes names present
                avail_att_s = avail_att_s | set(attrs_d.keys())

                # Adding to already parsed entry (without the full
                # attribute list i.e [:-1]
                entry = entry[:-1] + [attrs_d["gene_id"]] + \
                    [attrs_d["gene_name"]]
                cur.execute("""INSERT INTO GTF (chromo,source,feature,start,end,score,
                    strand,frame,gene_id,gene_name) VALUES (?,?,?,?,?,?,?,?,?,?);""", entry)

    log.info("Attributes present in the GTF: " + ", ".join(avail_att_s))
    log.warning(
        "Currently, this function only add gene_id and gene_name (more attributes available)")
    con.commit()
    con.close()

    return outdb


def import_snps_location(vcf):
    """
    Create a zip from a vcf file with SNPs informations 
    """
    # Make the zip of snps:
    chros = []
    pos = []
    with open(vcf, "r") as f:
        for line in f:
            # Pass the header
            if line.startswith("#"):
                continue

            line = line.strip().split("\t")
            chros.append(line[0])
            pos.append(line[1])

    return zip(chros, pos)


def import_go_terms(go_terms_file):
    """
    Return a dictionnary from the GO terms annotations from Ensembl biomart.
    The go_term file should be a csv with col1: genes IDs, col3: GO Accessions
    col4: GO name
    Only keep the biological processes
    """
    go_dict = {}
    with open(go_terms_file, "r") as f:
        next(f)  # pass the header

        for line in f:
            entry = line.strip().split(",")
            gene = entry[0]
            GO_id = entry[2]

            if gene in go_dict:
                go_dict[gene].append(GO_id)
            else:
                go_dict[gene] = [GO_id]

    return go_dict


def getFsts(fst_f):
    """
    From a FST file out of VCFtools ($1: chromo, $2:pos, $3: fsts),
    return a dictionnary:
    key: zip(chromo, pos)
    value: fst
    """
    fst_d = {}

    with open(fst_f) as f:
        next(f)  # pass the header

        for line in f:
            line = line.strip().split("\t")
            chro = line[0]
            pos = line[1]
            fst = float(line[2])

            fst_d[(chro, pos)] = fst

    return fst_d


def fetch_snps_location(db, snps_zip, go_dict, out_prefix, fst_f=False):
    """
    Fetch the annotations information from a list of snps positions
    """
    outfile = out_prefix + "_annot.csv"
    con = lite.connect(db)
    con.text_factory = str  # to not have unicode formatting

    count = 1

    if os.path.isfile(outfile):
        raise IOError(
            "The outfile '%s' already exists and will not be overwritten" % outfile)

    log.info("CREATING snps annotation file: %s" % outfile)

    if fst_f:
        fst_d = getFsts(fst_f)

    with open(outfile, "w") as fout:

        if fst_f:
            fout.write(
                "SNP\tCHRO\tPOS\tGID\tGNAME\tG_CHRO(DBUG)\tSTART\tEND\tFST\tGOs\n")
        else:
            fout.write(
                "SNP\tCHRO\tPOS\tGID\tGNAME\tG_CHRO(DBUG)\tSTART\tEND\tGOs\n")

        no_GO = []
        for c, p in snps_zip:
            fout.write("%s\t%s\t%s\t" % (count, c, p))
            cur = con.execute("""SELECT chromo, start, end, gene_id, gene_name from GTF
            WHERE chromo=? AND start <= ? AND end >= ?""", (c, p, p))
            data = cur.fetchone()
            if data:
                chromo, start, end, gene_id, gene_name = data
                GOs = "NO_GOs"
                try:
                    GOs = go_dict[gene_id]
                except KeyError:
                    # Nothing done yet with list of SNPs without GOs
                    no_GO.append(count)

                if fst_f:
                    fst = fst_d[c, p]
                    fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_id, gene_name,
                                                                 chromo, start, end, fst, GOs))
                else:
                    fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_id, gene_name,
                                                             chromo, start, end, GOs))
            else:
                if fst_f:
                    fst = fst_d[c, p]
                    fout.write(
                        "NA\tNA\tNA\tNA\tNA\t{0}\t{1}\n".format(fst, ['']))
                else:
                    fout.write("NA\tNA\tNA\tNA\tNA\tNA\tNA\n")

            count += 1

        con.close()


def snpPipe(gtf, vcf, go_terms, out_prefix, gq_thres):
    """Do the complete pipe for SNPs analysis from a gtf files, a vcf
    files and a csv file from ensembl with GO terms annotations
    """
    db = create_gtf_db(gtf)
    snps_zip = import_snps_location(vcf)
    go_dict = import_go_terms(go_terms)
    fetch_snps_location(db, snps_zip, go_dict, out_prefix)
    wvcf.print_genotypes(vcf, gq_thres, out_prefix, db)


def snpPipeFsts(gtf, vcf, go_terms, out_prefix, gq_thres, fst_f):
    """Do the complete pipe for SNPs analysis from a gtf files, a vcf
    files and a csv file from ensembl with GO terms annotations.
    In addition add the fsts values
    """
    db = create_gtf_db(gtf)
    snps_zip = import_snps_location(vcf)
    go_dict = import_go_terms(go_terms)
    fetch_snps_location(db, snps_zip, go_dict, out_prefix, fst_f)
    wvcf.print_genotypes(vcf, gq_thres, out_prefix, db)


@click.command()
@click.argument('gtf')
@click.argument('vcf')
@click.argument('go_terms')
@click.argument('out_prefix')
@click.argument('gq_thres', type=float)
@click.argument('fst_f')
def main(gtf, vcf, go_terms, out_prefix, gq_thres, fst_f):

    snpPipeFsts(gtf, vcf, go_terms, out_prefix, gq_thres, fst_f)


if __name__ == "__main__":

    main()
