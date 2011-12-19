import os, subprocess, re
from config import blast_prefs, fixed_dirs, run_dirs, r_root_dir
from common import ensure_dir
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
    fas_dir = fixed_dirs['mfas_contigs_dir']
    db_dir = fixed_dirs['blast_db_dir']
    ensure_dir([fas_dir, db_dir])
    g_name = genome['name']
    # make DB
    make_blastDB(db_dir+g_name, fas_dir+g_name+'_contigs.fas', 'nucl')

def make_ref_DB(run_id, ref):
    """Make a Blast DB from a reference FastA file."""
    # load inputs
    run_root = r_root_dir+run_id+"/"
    fas_dir = run_root+run_dirs['ref_fas_dir']
    db_dir = fixed_dirs['blast_db_dir']
    ensure_dir([fas_dir, db_dir])
    r_name = ref['name']
    # make DB
    make_blastDB(db_dir+r_name, fas_dir+r_name+'.fas', 'nucl')

def basic_batch_blastn(genomes, run_ref, run_id, timestamp):
    """Send batch jobs to Blast. Muxes to multiple reference DBs."""
    # load inputs
    ref_n = run_ref.name
    run_root = r_root_dir+run_id+"/"
    in_root = run_root+run_dirs['ref_seg_dir']+ref_n+"/"
    print " ", ref_n
    # log
    logstring = "".join(["\n\n# Blast segs to genomes @", timestamp, "\n\n"])
    run_ref.log(logstring)
    # do blastn
    for seg in run_ref.segs:
        input_file = in_root+ref_n+"_"+seg['name']+".fas"
        out_dir = run_root+run_dirs['blast_out_dir']+ref_n+"/"+seg['name']+"/"
        ensure_dir([out_dir])
        print "\t", seg['name'],
        for genome in genomes:
            g_name = genome['name']
            db_path = fixed_dirs['blast_db_dir']+g_name
            outfile = out_dir+g_name+"_out.txt"
            print ".",
            local_blastn_2file(input_file, db_path, outfile, blast_prefs)
        print ""
    run_ref.log("All OK")
    return "OK"

