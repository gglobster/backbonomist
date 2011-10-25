import re
import numpy as np
from os import listdir
from loaders import read_array, load_genbank
from writers import write_genbank
from config import directories as dirs, blast_dtypes, genomes, prox_D, mtype
from common import ensure_dir
from array_tetris import extract_nonzero, clump_rows
from Bio.Blast import NCBIXML

def glompX_blast_out(contig):
    """Collect Blast results and extract match contigs."""
    # load inputs
    nick = contig['name']
    match_root = dirs['match_out_dir']+nick+"/"
    print " ", nick
    # collect
    for ref in contig['refs']:
        print "\t", ref['type'], "...",
        blast_dir = dirs['blast_out_dir']+nick+"/"+ref['type']+"/"
        ensure_dir(blast_dir)
        for genome in genomes:
            g_name = genome['name']
            matches_dir = match_root+g_name+"/"
            ensure_dir(matches_dir)
            blast_infile = blast_dir+g_name+"_out.txt"
            genome_ctg_dir = dirs['gbk_contigs_dir']+g_name+"/"
            rec_array = read_array(blast_infile, blast_dtypes)
            if len(rec_array) > 0:  # if there are hits, take the best one
                contig_id = rec_array[0][1]
                #q_coords = rec_array[0][6], rec_array[0][7]
                #t_coords = rec_array[0][8], rec_array[0][9]
                print contig_id,
                pattern = re.compile(r'('+contig_id+')')
                for item in listdir(genome_ctg_dir):
                    match = re.match(pattern, item)
                    if match:
                        contig = load_genbank(genome_ctg_dir+item)
                        gbk_file = matches_dir+match.group(1)+".gbk"
                        write_genbank(gbk_file, contig)
        print ""

def mauver_load_k0(file, N_count):
    """Parse Mauve coordinates file to extract segment coordinates.

    This loads the coordinates data into a Numpy array. All rows that contain
    a zero value are deleted from the array. Rows are then collapsed if the
    coordinates are close (under a user-specified threshold) and in the same
    orientation. The resulting array is then returned.

    """
    # load file data into numpy array
    raw_array = np.loadtxt(file, skiprows=1, dtype=mtype)
    # set up an array stub to receive non-zero rows (leaves a (1,1) pair)
    stub_array = np.ones(1, dtype=mtype)
    # eliminate rows containing elements of value 0
    nz_array = extract_nonzero(raw_array, stub_array)
    # collapse rows
    cl_array = clump_rows(nz_array, prox_D)
    return cl_array

def collect_cogs(blast_out):
    results = {}
    blast_records = NCBIXML.parse(open(blast_out))
    for record in blast_records:
        rec_key = record.query_id
        results[rec_key+'_def'] = record.query
        if record.alignments: # ignores searches with no hits
            top_hit = record.alignments[0]
            results[rec_key] = top_hit.title
        else:
            results[rec_key] = 'no match'
    return results