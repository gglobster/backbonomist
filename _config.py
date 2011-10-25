import numpy

# input data list
g_names = ['']
plasmids = [{'name': '', 'refs': ({'coords': (0,0), 'type': ''})}]

# concat contigs separator
separator = ""

# project root directory
g_root_dir = 'data/genomes/'
p_root_dir = 'data/plasmids/'

directories = {
'ori_g_dir': g_root_dir+'original/',
'gbk_contigs_dir': g_root_dir+'genbank_ctgs/',
'mfas_contigs_dir': g_root_dir+'mfasta_ctgs/',
'blast_db_dir': g_root_dir+'blast_db/',
'ref_ctg_dir': p_root_dir+'ref_contigs/',
'ref_seg_dir': p_root_dir+'ref_segments/',
'blast_out_dir': p_root_dir+'blast_out/',
'match_out_dir': p_root_dir+'matches/',
}

# Blast parameters
blast_prefs = {'evalue': 0.001,
               'outfmt_pref': 6}
# Blast results arrays datatypes
blast_dtypes = numpy.dtype([('query', 'S16'),
                           ('dbhit', 'S32'),
                           ('idp', 'float'),
                           ('mlen', 'uint8'),
                           ('mms', 'uint8'),
                           ('gaps', 'uint8'),
                           ('q_start', 'uint16'),
                           ('q_end', 'uint16'),
                           ('r_start', 'uint16'),
                           ('r_end', 'uint16'),
                           ('evalue', 'S5'),
                           ('bitscore', 'float')])
