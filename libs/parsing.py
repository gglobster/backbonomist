import re, numpy as np
from os import listdir, path
from shutil import copyfile
from loaders import read_array, td_txt_file_load
from config import directories as dirs, p_root_dir, blast_dtypes, genomes, \
    mtype, min_match
from common import ensure_dir
from array_tetris import extract_nonzero, clump_rows
from Bio.Blast import NCBIXML

def glompX_blast_out(ref_ctg, run_id):
    """Collect Blast results and extract match contigs."""
    # load inputs
    nick = ref_ctg['name']
    run_root = p_root_dir+run_id+"/"
    match_root = run_root+dirs['match_out_dir']+nick+"/"
    print " ", nick
    # collect
    for ref in ref_ctg['refs']:
        print "\t", ref['name'], "...",
        blast_dir = run_root+dirs['blast_out_dir']+nick+"/"+ref['name']+"/"
        ensure_dir([blast_dir])
        for genome in genomes:
            g_name = genome['name']
            print g_name,
            matches_dir = match_root+g_name+"/"
            ensure_dir([matches_dir])
            blast_infile = blast_dir+g_name+"_out.txt"
            genome_ctg_dir = dirs['fas_contigs_dir']+g_name+"/"
            rec_array = read_array(blast_infile, blast_dtypes)
            if len(rec_array) > 0:  # take all hits above min_match length
                print "+",
                for line in rec_array:
                    q_start, q_stop = line[6], line[7]
                    if abs(q_stop-q_start) > min_match:
                        contig_id = line[1]
                        pattern = re.compile(r'('+contig_id+')\.fas')
                        for item in listdir(genome_ctg_dir):
                            match = re.match(pattern, item)
                            if match:
                                fas_file = matches_dir+match.group(1)+".fas"
                                if not path.exists(fas_file):
                                    copyfile(genome_ctg_dir+item, fas_file)
            else:
                print "-",
        print ""

def mauver_load2_k0(file, threshold):
    """Parse Mauve coordinates file to extract segment coordinates.

    This loads the coordinates data into a Numpy array. All rows that contain
    a zero value are deleted from the array. Rows are then collapsed if the
    coordinates are close (under a user-specified threshold) and in the same
    orientation. The resulting array is then returned.

    Important note: this function only works for pairwise alignments!

    """
    # load file data into numpy array
    raw_array = np.loadtxt(file, skiprows=1, dtype=mtype)
    # set up an array stub to receive non-zero rows (leaves a (1,1) pair)
    stub_array = np.ones(1, dtype=mtype)
    # eliminate rows containing elements of value 0
    try:
        nz_array = extract_nonzero(raw_array, stub_array)
    except TypeError:
        nz_array = np.append(stub_array, raw_array)
    # collapse rows
    cl_array = clump_rows(nz_array, threshold)
    return cl_array

def collect_cogs(blast_out):
    """Collect hits from Blast XML (not just for COGs anymore)."""
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

def parse_clustal_idstars(filename):
    """Parse ClustalW output file to estimate identity percentage."""
    #from analysis.text_manipulation import td_txt_file_load
    raw_lines = td_txt_file_load(filename, 3)
    single_line = ''.join(raw_lines)
    idnstar = re.compile(r"\*")
    idnalls = re.findall(idnstar, single_line)
    idntot = len(idnalls)
    return idntot # total number of identical nucleotide positions
