import numpy

# TODO: setup script to generate config file

# input data

project_id = 'BCSL'

genomes = [{'name': 'NZ_ACNE0', 'input': 'mfas', 'file': 'NZ_ACNE0.fas'},
           {'name': 'NZ_ACNF0', 'input': 'mfas', 'file': 'NZ_ACNF0.fas'},
           {'name': 'NZ_ACNI0', 'input': 'mfas', 'file': 'NZ_ACNI0.fas'},
           {'name': 'NZ_ACNJ0', 'input': 'mfas', 'file': 'NZ_ACNJ0.fas'},
           {'name': 'NZ_ACNK0', 'input': 'mfas', 'file': 'NZ_ACNK0.fas'}]

backbones = [{'name': 'pXO1', 'file': 'pXO1.gbk',
              'refs': ({'coords': (260,4650),'type': 'rep'},
                       {'coords': (74185,79173), 'type': 'hel'})}]
                        # TODO: make coords loader

# function categories and legend (keys MUST be lowercase)

fct_flags = {'mge': ('transposase', 'bla'),
             'tra': ('type iv', 'topoisomerase'),
             'ctl': ('transcriptional regulator', 'bla'),
             'unk': ('uncharacterized', 'bla'),
             'def': ('no match', 'bla')} # default

fct_colors = {'mge': ('#66CC00', 'MGE'),
              'tra': ('#6666FF', 'Transfer'),
              'ctl': ('#FFCC00', 'Control'),
              'oth': ('#99FFFF', 'Other'),
              'unk': ('#CCCCCC', 'Uncharacterized'),
              'def': ('#FFFFFF', 'No match')}

# concat contigs separator
separator = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

# project root directory
g_root_dir = 'data/'+project_id+'/genomes/'
b_root_dir = 'data/'+project_id+'/backbones/'

directories = {
'ori_g_dir': g_root_dir+'original/',
'gbk_contigs_dir': g_root_dir+'genbank_ctgs/',
'mfas_contigs_dir': g_root_dir+'mfasta_ctgs/',
'blast_db_dir': g_root_dir+'blast_db/',
'annot_trn_dir': g_root_dir+'annot_trn/',
'ref_ctg_dir': b_root_dir+'ref_contigs/',
'ref_seg_dir': b_root_dir+'ref_segments/',
'blast_out_dir': b_root_dir+'blast_out/',
'match_out_dir': b_root_dir+'matches/',
'scaffolds_dir': b_root_dir+'scaffolds/',
'mauve_out_dir': b_root_dir+'mauve_out/',
'scaff_annot_dir': b_root_dir+'annotation/',
'maps_dir': b_root_dir+'maps/'
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

# Mauve executable
mauve_exec = 'progressiveMauve'
# datatypes
mtype = [('A', 'i4'), ('B', 'i4'), ('C', 'i4'), ('D', 'i4')]
# Proximity threshold for clumping
prox_D = 2000
