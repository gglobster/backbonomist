import os, subprocess, re
from config import blast_prefs, directories as dirs, genomes
from loaders import load_multifasta
from common import ensure_dir
from datetime import datetime
from Bio.Blast.Applications import NcbiblastnCommandline, \
    NcbirpsblastCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIWWW

def make_blastDB(name, infile, db_type):
    """Make BLAST database from FASTA input file."""
    cline = "makeblastdb -in "+ infile +" -dbtype "+ db_type +" -title "+  \
            infile +" -out "+ name +" -parse_seqids"
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def local_blastn_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastn against local database."""
    cline = NcbiblastnCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=prefs['outfmt_pref'])
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_blastp_2file(query_file, dbfile_path, outfile, prefs):
    """Perform blastp against local database."""
    cline = NcbiblastpCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=5) # must output XML!
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def local_rpsblast_2file(query_file, dbfile_path, outfile, prefs):
    """Perform RPS Blast against local database."""
    cline = NcbirpsblastCommandline(query=query_file,
                                  db=dbfile_path,
                                  out=outfile,
                                  evalue=prefs['evalue'],
                                  outfmt=5) # must output XML!
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate() # forces the main script to wait

def remote_blastp_2file(query_string, database, outfile, evalue):
    """Perform blastp against remote database."""
    result_handle = NCBIWWW.qblast('blastp',
                                   database,
                                   query_string,
                                   expect=evalue)
    save_file = open(outfile, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()

def make_genome_DB(genome):
    """Make a Blast DB from a genome FastA file."""
    # load inputs
    fas_dir = dirs['mfas_contigs_dir']
    db_dir = dirs['blast_db_dir']
    g_name = genome['name']
    print " ", g_name, "...",
    # make DB
    make_blastDB(db_dir+g_name, fas_dir+g_name+'_contigs.fas', 'nucl')
    print "DB ready"

def basic_batch_blastn(contig):
    """Send batch jobs to Blast. Muxes to multiple reference DBs."""
    # load inputs
    nick = contig['name']
    in_root = dirs['ref_seg_dir']+nick+"/"
    print " ", nick
    # do blastn
    for ref in contig['refs']:
        input_file = in_root+nick+"_"+ref['type']+".fas"
        out_dir = dirs['blast_out_dir']+nick+"/"+ref['type']+"/"
        ensure_dir(out_dir)
        print "\t", ref['type'], "...",
        for genome in genomes:
            g_name = genome['name']
            db_path = dirs['blast_db_dir']+g_name
            outfile = out_dir+g_name+"_out.txt"
            print g_name,
            local_blastn_2file(input_file, db_path, outfile, blast_prefs)
        print ""
    return "OK"

